#!/bin/env python

def guessObservables(ws):
    import ROOT
    observables = set()
    for func in ws.allFunctions():
        if func.InheritsFrom(ROOT.RooHistFunc.Class()):
            for obs in func.dataHist().get():
                observables.add(obs.GetName())
    for pdf in ws.allPdfs():
        if pdf.InheritsFrom(ROOT.RooHistPdf.Class()):
            for obs in pdf.dataHist().get():
                observables.add(obs.GetName())
        elif not "Constraint" in pdf.GetName():
            args = ROOT.RooArgSet()
            if pdf.getAnalyticalIntegral(ws.allVars(),args) > 0:
                for obs in args:
                    observables.add(obs.GetName())
            
    for idx in ws.allCats():
        observables.add(idx.GetName())
    return list(observables)
        

def main(args):
    import ROOT

    etcDir = str(ROOT.TROOT.GetEtcDir())

    for arg in args.import_defs:
        ROOT.RooFit.JSONIO.loadFactoryExpressions(arg)

    
    from os.path import join as pjoin
    ws = ROOT.RooWorkspace(args.name)
    tool = ROOT.RooJSONFactoryWSTool(ws)

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
    
    try:
        importJSON = tool.importJSON.__overload__("const string&")
        importJSON(args.infile)
    except Exception as theError:
        print(theError.what())

    if args.simPdf:
        toplevel = {}
        index = ROOT.RooCategory("indexCat","indexCat")
        for pdf in ws.allPdfs():
            if len(pdf.clients()) == 0:
                toplevel[pdf.GetName()] = pdf
        simPdf = ROOT.RooSimultaneous("simPdf","simPdf",toplevel,index)
        ws.Import(simPdf,ROOT.RooFit.RecycleConflictNodes(True))

    if args.modelConfig:
        toplevel = []
        for pdf in ws.allPdfs():
            if len(pdf.clients()) == 0:
                toplevel.append(pdf)
        if len(toplevel) == 1:
            origMC = ROOT.RooStats.ModelConfig("ModelConfig","ModelConfig")
            origMC.SetWS(ws)
            ws.Import(origMC)
            mc = ws.obj("ModelConfig")
            mc.SetPdf(toplevel[0])
            mc.SetParametersOfInterest(",".join(args.pois))
            mc.SetObservables(",".join(guessObservables(ws)))
            for obs in mc.GetObservables():
                if obs.InheritsFrom(ROOT.RooRealVar.Class()):
                    obs.setConstant(True)
        else:
            print("cannot create ModelConfig, no unique top-level p.d.f. found. Please create one with --simPdf")

    if args.asimov:
        toplevel = []
        for mc in ws.allGenericObjects():
            if mc.InheritsFrom(ROOT.RooStats.ModelConfig.Class()):
                data = ROOT.RooStats.AsymptoticCalculator.GenerateAsimovData(mc.GetPdf(),mc.GetObservables())
                data.SetName("asimovData_"+pdf.GetName())
                ws.Import(data,ROOT.RooFit.RecycleConflictNodes(True))
                
    if args.generate:
        toplevel = []
        for mc in ws.allGenericObjects():
            if mc.InheritsFrom(ROOT.RooStats.ModelConfig.Class()):
                data = mc.GetPdf().generate(mc.GetObservables())
                data.SetName("toyData_"+pdf.GetName())
                ws.Import(data,ROOT.RooFit.RecycleConflictNodes(True))                
                
        
        
    ws.writeToFile(args.outfile)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="turn JSON into a workspace")
    parser.add_argument("-i",dest="infile",help="the JSON file")
    parser.add_argument("-o",dest="outfile",help="the ROOT workspace")
    parser.add_argument("--load-importers",nargs="+",type=str,help="load some importers",dest="import_defs",default=[])
    parser.add_argument("--pois",type=str,nargs="+",default="mu")    
    parser.add_argument("--asimov",action="store_true")
    parser.add_argument("--generate",action="store_true")        
    parser.add_argument("--modelConfig",action="store_true")
    parser.add_argument("--simPdf",action="store_true")
    parser.add_argument("--name",help="name of the workspace",default="workspace")

    main(parser.parse_args())
