#!/bin/env python

def main(args):
    import ROOT

    for lib in args.libraries:
        ROOT.gSystem.Load(lib)
    
    infile = ROOT.TFile.Open(args.infile,"READ")

    if args.recover_file:
        infile.Recover()
        infile.ReadStreamerInfo()
    
    ws = None
    for k in infile.GetListOfKeys():
        if k.GetClassName()=="RooWorkspace":
            ws = k.ReadObj()
            break
    if not ws:
        print("no workspace found in file!")
        exit(1)

    if args.default_defs:
        etcDir = str(ROOT.TROOT.GetEtcDir())
        from os.path import join as pjoin
        ROOT.RooJSONFactoryWSTool.loadExportKeys(pjoin(etcDir,"RooFitHS3_wsexportkeys.json"))
        ROOT.RooJSONFactoryWSTool.loadFactoryExpressions(pjoin(etcDir,"RooFitHS3_wsfactoryexpressions.json"))
    for arg in args.export_defs:
        ROOT.RooJSONFactoryWSTool.loadExportKeys(arg)

    if args.print:
        ws.Print()

    if args.patchMC:
        mc = ws.obj(args.patchMC)
        mc.SetWS(ws)
    
    tool = ROOT.RooJSONFactoryWSTool(ws)

    tool.exportJSON(args.outfile)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="turn a workspace into JSON")
    parser.add_argument("infile",help="the ROOT workspace")
    parser.add_argument("outfile",help="the JSON file")
    parser.add_argument("--load-libraries",nargs="+",type=str,help="load some libraries",dest="libraries",default=[])
    parser.add_argument("--recover",action="store_true",dest="recover_file",default=False,help="recover corrupted file")    
    parser.add_argument("--no-default-definitions",action="store_false",dest="default_defs",help="avoid loading the default definitions")
    parser.add_argument("--default-definitions",action="store_true",dest="default_defs",default=True,help="load the default definitions")
    parser.add_argument("--print",action="store_true",dest="print",default=False,help="print the workspace once loaded")
    parser.add_argument("--patch-model-config",dest="patchMC",default=None,help="Model Config to be patched")    
    parser.add_argument("--load-definitions",dest="export_defs",nargs="+",default=[],help="load a specific set of definitions")

    main(parser.parse_args())    
