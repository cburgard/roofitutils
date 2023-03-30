#!/bin/env python

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
        tool.importJSON(args.infile)
    except Exception as theError:
        print(theError)
    ws.writeToFile(args.outfile)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="turn JSON into a workspace")
    parser.add_argument("-i",dest="infile",help="the JSON file")
    parser.add_argument("-o",dest="outfile",help="the ROOT workspace")
    parser.add_argument("--load-importers",nargs="+",type=str,help="load some importers",dest="import_defs",default=[])        
    parser.add_argument("--name",help="name of the workspace",default="workspace")

    main(parser.parse_args())
