#!/bin/env python

def main(args):
    import ROOT

    etcDir = str(ROOT.TROOT.GetEtcDir())
    from os.path import join as pjoin
    ROOT.RooJSONFactoryWSTool.loadExportKeys(pjoin(etcDir,"RooFitHS3_wsexportkeys.json"))
    ROOT.RooJSONFactoryWSTool.loadFactoryExpressions(pjoin(etcDir,"RooFitHS3_wsfactoryexpressions.json"))
    ws = ROOT.RooWorkspace(args.name)
    tool = ROOT.RooJSONFactoryWSTool(ws)
    tool.importJSON(args.infile)
    ws.writeToFile(args.outfile)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="turn JSON into a workspace")
    parser.add_argument("infile",help="the JSON file")
    parser.add_argument("outfile",help="the ROOT workspace")
    parser.add_argument("--name",help="name of the workspace",default="workspace")

    main(parser.parse_args())
