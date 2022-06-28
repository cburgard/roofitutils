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
    for arg in args.import_defs:
        ROOT.RooJSONFactoryWSTool.loadFactoryExpressions(arg)        
    tool = ROOT.RooJSONFactoryWSTool(ws)
    tool.exportJSON(args.outfile)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="turn a workspace into JSON")
    parser.add_argument("infile",help="the ROOT workspace")
    parser.add_argument("outfile",help="the JSON file")
    parser.add_argument("--load-libraries",nargs="+",type=str,help="load some libraries",dest="libraries",default=[])
    parser.add_argument("--recover",action="store_true",dest="recover_file",default=False,help="recover corrupted file")    
    parser.add_argument("--no-default-definitions",action="store_false",dest="default_defs")
    parser.add_argument("--default-definitions",action="store_true",dest="default_defs",default=True)
    parser.add_argument("--load-export-definitions",dest="export_defs",nargs="+",default=[])
    parser.add_argument("--load-import-definitions",dest="import_defs",nargs="+",default=[])    

    main(parser.parse_args())    
