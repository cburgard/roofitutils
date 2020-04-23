
def main(args):
    import ROOT
    infile = ROOT.TFile.Open(args.input)
    ws = infile.Get(args.name)
    mc = ws.obj(args.mc)
    pdf = mc.GetPdf()
    pdfname = pdf.GetName()
    pdf.SetName(pdfname+"_old")

    oldpoi = ws.var(args.old)
    if not oldpoi:
        print("unable to find parameter '"+args.old+"', parmeters are:")
        ws.allVars().Print("v")
        exit(0)
    splits = {}
    for k in args.splits:
        elems = k.split(":")
        splits[elems[0]]=elems[1]
    for k,v in splits.items():
        ws.factory("{:s}[{:g},{:g},{:g}]".format(v,oldpoi.getVal(),oldpoi.getMin(),oldpoi.getMax()))

    toreplace = []
    for f in ws.allFunctions():
        if not f.dependsOn(oldpoi): continue
        for g in toreplace:
            if f.dependsOn(g):
                continue
            if g.dependsOn(f):
                toreplace.remove(g)
        toreplace.append(f)

    pdfreplace = {}
    for f in toreplace:
        fname = f.GetName()
        f.SetName(fname+"_old")
        pdfreplace[f.GetName()]=fname
        newpar = None
        for k,v in splits.items():
            if k in fname:
                newpar = v
        if not newpar:
            raise RuntimeError("function {:s} does not fit any of the categories!".format(fname))
        s = "EDIT::{:s}({:s},{:s}={:s})".format(fname,f.GetName(),oldpoi.GetName(),newpar)
        print("replacing {:s} by {:s} in {:s}".format(oldpoi.GetName(),newpar,fname))
        ws.factory(s)

    ws.factory("EDIT::{:s}({:s},".format(pdfname,pdf.GetName())+",".join([k+"="+v for k,v in pdfreplace.items()])+")")
    newpdf = ws.pdf(pdfname)
    mc.SetPdf(newpdf)
    if newpdf.dependsOn(oldpoi):
        print("new Pdf still depends on old parameter - not all occurences were replaced!")
    else:
        print("all occurences were replaced!")

    import RooFitUtils
    newws = ROOT.RooFitUtils.makeCleanWorkspace(ws,"combined")
    newws.writeToFile(args.output)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("split a parameter into several independent ones based on region and/or sample",formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("input", help="Input workspace to be used")
    parser.add_argument("--output", default="workspace-output.root", help="Output workspace")
    parser.add_argument("--name", default="combined", help="Workspace name")
    parser.add_argument("--mc", default="ModelConfig", help="ModelConfig name")
    parser.add_argument("--old", required=True, help="old parameter name")
    parser.add_argument("--splits", nargs="+", help="list of splits to be made", metavar="key:newparameter")

    args = parser.parse_args()
    main(args)
