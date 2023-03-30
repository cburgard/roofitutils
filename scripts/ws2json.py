#!/bin/env python

def main(args):
    import ROOT

    for lib in args.libraries:
        ROOT.gSystem.Load(lib)

    ws = None

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)    
        
    for infilename in args.infile:
        thisws = None
        infile = ROOT.TFile.Open(infilename,"READ")

        if args.recover_file:
            infile.Recover()
            infile.ReadStreamerInfo()
    
        for k in infile.GetListOfKeys():
            if k.GetClassName()=="RooWorkspace":
                thisws = k.ReadObj()
                break
            if ROOT.TClass.GetClass(k.GetClassName()).InheritsFrom(ROOT.RooAbsArg.Class()):
                ws.Import(k.ReadObj())

        if not thisws:
            continue

        if ws:
            for pdf in thisws.allPdfs():
                ws.Import(pdf,ROOT.RooFit.RecycleConflictNodes())
            for func in thisws.allFunctions():
                ws.Import(func,ROOT.RooFit.RecycleConflictNodes())
        else:
            ws = thisws.Clone()
            if args.patchMC:
                mc = ws.obj(args.patchMC)
                mc.SetWS(ws)
                if mc.GetPdf():
                    if not mc.GetPdf().getStringAttribute("combined_data_name"):
                        mc.GetPdf().setStringAttribute("combined_data_name","obsData")

    if not ws:
        print("no workspace found in file!")
        exit(1)

    for arg in args.export_defs:
        ROOT.RooFit.JSONIO.loadExportKeys(arg)

    if args.print:
        ws.Print()

    tool = ROOT.RooJSONFactoryWSTool(ws)

    tool.exportJSON(args.outfile)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="turn a workspace into JSON")
    parser.add_argument("-i",dest="infile",help="the ROOT workspace(s)",nargs="+")
    parser.add_argument("-o",dest="outfile",help="the JSON file")
    parser.add_argument("--load-exporters",nargs="+",type=str,help="load some exporters",dest="export_defs",default=[])    
    parser.add_argument("--load-libraries",nargs="+",type=str,help="load some libraries",dest="libraries",default=[])
    parser.add_argument("--recover",action="store_true",dest="recover_file",default=False,help="recover corrupted file")
    parser.add_argument("--print",action="store_true",dest="print",default=False,help="print the workspace once loaded")
    parser.add_argument("--patch-model-config",dest="patchMC",default=None,help="Model Config to be patched")    

    main(parser.parse_args())    
