#!/bin/env python

def find_elements(thelist,theregexlist):
    import re
    from ROOT import RooArgList
    elems = RooArgList()
    for expr in theregexlist:
        rex = re.compile(expr)
        for obj in thelist:
            matched = rex.match(obj.GetName())
            if matched:
                elems.add(obj)
    return elems

def linearize(func,p,e):
    from ROOT import RooFormulaVar
    c = func.getVal()
    v = p.getVal()
    p.setVal(v+e)
    cplus = func.getVal()
    p.setVal(v-e)    
    cminus = func.getVal()
    p.setVal(v)        
    name = func.GetName()+"_lin_"+p.GetName()
    formula = "{:s}>0 ? 1 + {:s}*{:g} : 1 - {:s}*{:g}".format(p.GetName(),p.GetName(),cplus/c-1,p.GetName(),cminus/c-1)
    return RooFormulaVar(name,name,formula,p)


def get_parametrization(polyfunc):
    parametrization = {}
    terms = polyfunc.terms()
    unity = None
    for term in terms:
        if not unity:
            unity = polyfunc.getCoefficient(term)
            continue
        coef = polyfunc.getCoefficient(term)
        for v in polyfunc.variables():
            e = polyfunc.getExponent(term,v)
            if e.getVal() == 1:
                parametrization[v.GetName()] = coef.getVal()
    return parametrization

def main(args):
    from ROOT import TFile
    infile = TFile.Open(args.input,"READ")
    if not infile or not infile.IsOpen():
        print("unable to open file "*args.input)
        exit(1)
    workspace = infile.Get(args.workspace)
    if not workspace:
        print("file "+args.input+" does not have workspace "+args.workspace)
        exit(1)
    modelConfig = workspace.obj(args.model)
    if modelConfig:
        pdf = modelConfig.GetPdf()
    else:
        pdf = workspace.pdf(args.pdf)
    if not pdf:
        print("workspace does not have a ModelConfig named '"+args.model+"' or a pdf named '"+args.pdf+"', please specify!")
        exit(1)

    xsecs = find_elements(workspace.allFunctions(),args.XS)
    pois = find_elements(workspace.allVars(),args.POIs)    

    parametrization = {}

    from ROOT import RooPolyFunc,RooArgSet
    allnps = RooArgSet()
    allobs = RooArgSet()
    if modelConfig:
        allnps.add(modelConfig.GetNuisanceParameters())
        allobs.add(modelConfig.GetObservables())        
    else:
        for np in pdf.getParameters(pois):
            if np.isConstant():
                allobs.add(np)
            else:
                allnps.add(np)
    cmd = "EDIT::" + pdf.GetName()+"_unparametrized" + "(" + pdf.GetName()
    for xs in xsecs:
        poilist = RooArgSet()
        for p in pois:
            if xs.dependsOn(p):
                poilist.add(p)
        if not len(poilist):
            continue
        param = RooPolyFunc.taylorExpand(str(xs.GetName()+"_parametrization"),str(xs.GetTitle()),xs,poilist)
        local_parametrization =  get_parametrization(param)
        local_parametrization["SM"] = xs.getVal()
        parametrization[xs.GetName()+"_raw"] = local_parametrization
        nps = xs.getParameters(poilist)
        allnps.add(nps)
        unparam = []
        if len(nps) > 0:
            for np in nps:
                lin = linearize(xs,np,1)
                workspace.Import(lin)
                unparam.append(lin.GetName())
        else:
            unparam = None
        if unparam:
            replacement = workspace.factory("RooProduct::"+xs.GetName()+"_modified({" + "{:s}_raw[{:g}],".format(xs.GetName(),xs.getVal()) + ",".join(unparam)+"})")
        else:
            replacement = workspace.factory("{:s}_raw[{:g}]".format(xs.GetName(),xs.getVal()))
        cmd += ","+xs.GetName()+"="+replacement.GetName()
    cmd += ")"

    newpdf = workspace.factory(cmd)
    
    if args.write_yml:
        yml = {"name":pdf.GetName(),"parameterisation":parametrization,"coefficient terms":[p.GetName() for p in pois]+["SM"],"observables":[xs.GetName()+"_raw" for xs in xsecs]}
        import yaml
        with open(args.write_yml,"wt") as outfile:
            documents = yaml.dump(yml, outfile)
            
    if args.write_root:
        from ROOT import RooWorkspace,RooStats
        ws = RooWorkspace(workspace.GetName())
        ws.Import(newpdf)
        mc = RooStats.ModelConfig("ModelConfig")
        mc.SetWS(ws)
        mc.SetPdf(newpdf)
        mc.SetNuisanceParameters(allnps)
        mc.SetObservables(allobs)        
        mc.SetParametersOfInterest(pois)
        ws.Import(mc)
        for data in workspace.allData():
            ws.Import(data)
        ws.writeToFile(args.write_root)

        

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="remove (and extract) a parametrization from a workspace")
    parser.add_argument("--input","-i",help="name of the input file",required=True)
    parser.add_argument("--workspace","-w",help="name of the workspace in the input file",required=True)
    parser.add_argument("--model",default="ModelConfig",help="name of the ModelConfig object")
    parser.add_argument("--pdf",default="simPdf",help="name of the top-level pdf in the workspace")
    parser.add_argument("--POIs",nargs="+",required=True)
    parser.add_argument("--write-yml",type=str,default=None)
    parser.add_argument("--write-root",type=str,default=None)        
    parser.add_argument("--XS",nargs="+",required=True)
    parser.add_argument("--order",type=int,default=1,help="order to which to expand")
    args = parser.parse_args()
    main(args)
