#!/usr/bin/env python
def main(args):
    import gc
    gc.disable()

    import ROOT
    infile = ROOT.TFile.Open(args.input)
    if not infile or not infile.IsOpen():
        print("unable to open file "+args.input)
        exit(0)

    ws = infile.Get(args.workspace)
    if not ws:
        print("unable to obtain workspace "+args.workspace)
        exit(0)

    mc = ws.obj(args.mc)
    if not mc:
        print("unable to obtain model config "+args.mc)
        exit(0)

    if ws.data(args.name):
        print("a dataset with the name "+args.name+" is already present in the workspace!")
        exit(0)

    allVars = ws.allVars().snapshot()

    for assignment in args.parameters:
        p,v = assignment.split("=")
        par = allVars.find(p)
        if not par:
            print("unable to find parameter "+p)
            exit(0)
        val = float(v)
        if val > par.getMax():
            par.setMax(val)
        if val < par.getMin():
            par.setMin(val)
        par.setVal(val)

    from RooFitUtils.util import createAsimov
    createAsimov(ws,mc,args.name)

    ws.writeToFile(args.output)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("create an asimov dataset for a workspace",formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("input", help="Input workspace to be used")
    parser.add_argument("--output", default="workspace-output.root", help="Output workspace")
    parser.add_argument("--workspace", default="combined", help="Workspace name")
    parser.add_argument("--name", default="asimovData", help="Dataset name")
    parser.add_argument("--parameters", nargs="+", help="values of parameters to be used", metavar="par=x", default=[])
    parser.add_argument("--mc", default="ModelConfig", help="ModelConfig name")

    args = parser.parse_args()
    main(args)
