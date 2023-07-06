#!/bin/env python

def find_elements(thelist,theregexlist,theantiregexlist=[]):
    import re
    from ROOT import RooArgList
    elems = RooArgList()
    rexlist = [re.compile(e) for e in theregexlist]
    antirexlist = [re.compile(e) for e in theantiregexlist]
    for obj in thelist:
        matched = False
        for rex in rexlist:
            if rex.match(obj.GetName()):
                matched = True
        for rex in antirexlist:
            if rex.match(obj.GetName()):
                matched = False
        if matched:
            elems.add(obj)
    return elems

def sign(v):
    if v>0: return 1
    elif v<0: return -1
    return 0

def linearize(func,p,e,symmetrize=None):
    from ROOT import RooFormulaVar
    c = func.getVal()
    v = p.getVal()
    p.setVal(v+e)
    cplus = func.getVal()
    p.setVal(v-e)    
    cminus = func.getVal()
    p.setVal(v)
    dplus = cplus/c-1
    dminus = cminus/c-1
    if symmetrize == "max":
        if abs(dplus) >= abs(dminus):
            dminus = -dplus
        else:
            dplus = -dminus
    elif symmetrize == "avg":
        avg = 0.5*(abs(dplus)+abs(dminus))
        if abs(dplus) >= abs(dminus):
            dplus = sign(dplus) * avg
            dminus = -dplus
        else:
            dminus = sign(dminus) * avg
            dplus = -dminus
    name = func.GetName()+"_lin_"+p.GetName()
    if abs(dplus) == abs(dminus):
        formula = "1 + {:s}*{:g}".format(p.GetName(),dplus)
    else:
        formula = "{:s}>0 ? 1 + {:s}*{:g} : 1 - {:s}*{:g}".format(p.GetName(),p.GetName(),dplus,p.GetName(),dminus)
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
    pois = find_elements(workspace.allVars(),args.POIs,args.noPOIS)    

    parametrization = {}

    from ROOT import RooPolyFunc,RooArgSet,RooAbsPdf
    allnps = RooArgSet()
    allobs = RooArgSet()
    allconstraints = []
    if modelConfig:
        allnps.add(modelConfig.GetNuisanceParameters())
        allobs.add(modelConfig.GetObservables())        
    else:
        for np in pdf.getParameters(pois):
            if np.isConstant():
                allobs.add(np)
            else:
                allnps.add(np)
                for c in np.clients():
                    if c.InheritsFrom(RooAbsPdf.Class()):
                        allconstraints.append(c)
    cmd = "EDIT::" + pdf.GetName()+"_unparametrized" + "(" + pdf.GetName()
    newpois = []
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
                lin = linearize(xs,np,1,args.symmetrize)
                workspace.Import(lin)
                unparam.append(lin.GetName())
        else:
            unparam = None
        if unparam:
            replacement = workspace.factory("RooProduct::"+xs.GetName()+"_modified({" + "{:s}_raw[{:g}],".format(xs.GetName(),xs.getVal()) + ",".join(unparam)+"})")
        else:
            replacement = workspace.factory("{:s}_raw[{:g}]".format(xs.GetName(),xs.getVal()))
        cmd += ","+xs.GetName()+"="+replacement.GetName()
        newpois.append(xs.GetName()+"_raw")
    cmd += ")"

    newpdf = workspace.factory(cmd)
    localobs = RooArgSet()
    globalobs = RooArgSet()    
    for obs in allobs:
        glob = False
        for c in allconstraints:
            if c.dependsOn(obs):
                glob = True
        if glob:
            globalobs.add(obs)
        else:
            localobs.add(obs)
        

    if args.extend and not newpdf.canBeExtended():
        from ROOT import RooExtendPdf
        name = newpdf.GetName()+"_ext"
        extpdf = RooExtendPdf(name,name,newpdf,1.)
        newpdf = extpdf

    if args.write_yml:
        yml = {"name":pdf.GetName(),"parameterisation":parametrization,"coefficient terms":[p.GetName() for p in pois]+["SM"],"observables":newpois}
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
        mc.SetObservables(localobs)
        mc.SetGlobalObservables(globalobs)        
        mc.SetParametersOfInterest(",".join(newpois))
        for poi in newpois:
            ws.var(poi).setConstant(False)
        ws.Import(mc)
        for data in workspace.allData():
            ws.Import(data)
        
#        if not ws.data("asimovData"):
#            import ROOT
##            changeMsgLvl = ROOT.RooHelpers.LocalChangeMsgLevel(ROOT.RooFit.ERROR)
#            from RooFitUtils.util import createAsimov
#            ROOT.RooMsgService.instance().addStream(ROOT.RooFit.DEBUG, Topic = ROOT.RooFit.Generation)
#            createAsimov(ws,mc,"asimovData")
            
        allVars = RooArgSet()
        allVars.add(mc.GetParametersOfInterest())
        allVars.add(mc.GetNuisanceParameters())        
        ws.saveSnapshot("AllVars_Nominal",allVars)

        print("writing "+args.write_root)
        ws.writeToFile(args.write_root)

        

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="remove (and extract) a parametrization from a workspace")
    parser.add_argument("--input","-i",help="name of the input file",required=True)
    parser.add_argument("--workspace","-w",help="name of the workspace in the input file",required=True)
    parser.add_argument("--model",default="ModelConfig",help="name of the ModelConfig object")
    parser.add_argument("--pdf",default="simPdf",help="name of the top-level pdf in the workspace")
    parser.add_argument("--extend",action="store_true",help="extend the pdf if need be")
    parser.add_argument("--POIs",nargs="+",required=True)
    parser.add_argument("--no-POIs",nargs="+",dest="noPOIS",required=False)    
    parser.add_argument("--symmetrize",choices=["none","max","avg"],default="none")
    parser.add_argument("--write-yml",type=str,default=None)
    parser.add_argument("--write-root",type=str,default=None)        
    parser.add_argument("--XS",nargs="+",required=True)
    parser.add_argument("--order",type=int,default=1,help="order to which to expand")
    args = parser.parse_args()
    main(args)
