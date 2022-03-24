#!/bin/env python

def main(args):
    import ROOT
    infile = ROOT.TFile.Open(args.infile,"READ")
    ws = None
    for k in infile.GetListOfKeys():
        if k.GetClassName()=="RooWorkspace":
            ws = k.ReadObj()
            break
    if not ws:
        print("no workspace found in file!")
        exit(1)
    
    etcDir = str(ROOT.TROOT.GetEtcDir())
    from os.path import join as pjoin
    ROOT.RooJSONFactoryWSTool.loadExportKeys(pjoin(etcDir,"RooFitHS3_wsexportkeys.json"))
    ROOT.RooJSONFactoryWSTool.loadFactoryExpressions(pjoin(etcDir,"RooFitHS3_wsfactoryexpressions.json"))
    tool = ROOT.RooJSONFactoryWSTool(ws)
    tool.exportJSON(args.outfile)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="turn a workspace into JSON")
    parser.add_argument("infile",help="the ROOT workspace")
    parser.add_argument("outfile",help="the JSON file")

    main(parser.parse_args())    
