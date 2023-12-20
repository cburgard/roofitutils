#!/bin/env python

_memory = []
def protect(obj):
    _memory.append(obj)

def main(args):
    import ROOT

    for lib in args.libraries:
        ROOT.gSystem.Load(lib)

    try:
        ROOT.RooFit.JSONIO.setupKeys()
    except:
        print("unable to find default keys, skipping")
        pass

    ws = None
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)    
        
    for infilename in args.infile:
        thisws = None
        infile = ROOT.TFile.Open(infilename,"READ")
        protect(infile)

        if args.recover_file:
            infile.Recover()
            infile.ReadStreamerInfo()

        if ws and args.verbose:
            ws.Print()
            
        for k in infile.GetListOfKeys():
            if k.GetClassName()=="RooWorkspace":
                thisws = k.ReadObj()
                break
            if ROOT.TClass.GetClass(k.GetClassName()).InheritsFrom(ROOT.RooAbsArg.Class()):
                if not ws:
                    ws = ROOT.RooWorkspace("workspace")
                obj = k.ReadObj()
                ws.Import(obj.Clone())

        if not thisws:
            continue

        if ws:
            for pdf in thisws.allPdfs():
                ws.Import(pdf,ROOT.RooFit.RecycleConflictNodes())
            for func in thisws.allFunctions():
                ws.Import(func,ROOT.RooFit.RecycleConflictNodes())
        else:
            protect(ws)
            ws = thisws

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

    if args.force_rename:
        for var in ws.allVars():
            if not var.InheritsFrom(ROOT.RooConstVar().Class()) and "." in var.GetName():
                var.SetName(var.GetName().replace(".","_"))
        for func in ws.allFunctions():
            if "." in func.GetName():
                func.SetName(func.GetName().replace(".","_"))
        for pdf in ws.allPdfs():
            if "." in pdf.GetName():
                pdf.SetName(pdf.GetName().replace(".","_"))
        for data in ws.allData():
            if "." in data.GetName():
                data.SetName(data.GetName().replace(".","_"))
            for var in data.get():
                if "." in var.GetName():
                    var.SetName(var.GetName().replace(".","_"))
        for data in ws.allEmbeddedData():
            if "." in data.GetName():
                data.SetName(data.GetName().replace(".","_"))
            for var in data.get():
                if "." in var.GetName():
                    var.SetName(var.GetName().replace(".","_"))                                    
                    
    if args.print:
        ws.Print()

    extConstraints = ROOT.RooArgSet()
    for c in args.extConstraint:
        if not ws.obj(c):
            print("unable to find object "+c)
        else:
            extConstraints.add(ws.obj(c))
    if extConstraints.size() > 0:
        for obj in ws.allGenericObjects():
            if obj.InheritsFrom(ROOT.RooStats.ModelConfig.Class()):
                obj.SetExternalConstraints(extConstraints)
        
    tool = ROOT.RooJSONFactoryWSTool(ws)
    if not args.use_lists:
        tool.useListsInsteadOfDicts = args.use_lists

    tool.exportJSON(args.outfile)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="turn a workspace into JSON")
    parser.add_argument("-i",dest="infile",help="the ROOT workspace(s)",nargs="+")
    parser.add_argument("-v","--verbose",action="store_true",help="be more verbose")
    parser.add_argument("-o",dest="outfile",help="the JSON file")
    parser.add_argument("--load-exporters",nargs="+",type=str,help="load some exporters",dest="export_defs",default=[])    
    parser.add_argument("--load-libraries",nargs="+",type=str,help="load some libraries",dest="libraries",default=[])
    parser.add_argument("--recover",action="store_true",dest="recover_file",default=False,help="recover corrupted file")
    parser.add_argument("--force-rename",action="store_true",dest="force_rename",default=False,help="change forbidden names")    
    parser.add_argument("--print",action="store_true",dest="print",default=False,help="print the workspace once loaded")
    parser.add_argument("--patch-model-config",dest="patchMC",default=None,help="Model Config to be patched")
    parser.add_argument("--register-external-constraint",dest="extConstraint",nargs="+",default=[])
    parser.add_argument("--use-lists",dest="use_lists",action="store_true",default=True)
    parser.add_argument("--use-dicts",dest="use_lists",action="store_false")    

    main(parser.parse_args())    
