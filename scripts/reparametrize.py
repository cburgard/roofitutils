#!/usr/bin/env python

import re
varchars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"
varex    = re.compile("["+varchars+"]+")
assignex = re.compile("(?P<lhs>["+varchars+"]+)\\s*=\\s*(?P<rhs>.+)")
paramex  = re.compile("(?P<name>["+varchars+"]+)\[[0-9,\.]+\]")
exprex   = re.compile("expr::(?P<name>["+varchars+"]+)\s*\('(?P<expr>.*)',(?P<args>.+)\)")
prodex   = re.compile("prod::(?P<name>["+varchars+"]+)\s*\((?P<args>.+)\)")

def parse(expression):
    assignment = assignex.match(expression)
    product    = prodex.match(expression)
    expr       = exprex.match(expression)
    parameter  = paramex.match(expression)
    if parameter:
        return (parameter.group("name"),[])
    elif assignment:
        return (assignment.group("lhs"),sorted(set(varex.findall(assignment.group("rhs")))))
    elif product:
        return (product.group("name"),sorted(set(varex.findall(product.group("args")))))
    elif expr:
        return (expr.group("name"),sorted(set(varex.findall(expr.group("args")))))
    raise RuntimeError("unable to parse expression '"+expression+"'")

def namelist(coll):
    itr = coll.createIterator()
    var = itr.Next()
    retval = []
    while var :
        retval.append(var.GetName())
        var = itr.Next()
    return retval

def loadRooFitUtils():
    # retrieve the root core dir environment variable
    from os import getenv
    rcdir = getenv ("ROOTCOREDIR")
    if rcdir:
        from ROOT import gROOT
        if gROOT.ProcessLine(".x $ROOTCOREDIR/scripts/load_packages.C"):
            raise ImportError("unable to load RootCore!")
    else:
        from ROOT import gSystem
        if gSystem.Load("libRooFitUtils"):
            raise ImportError("unable to load standalone libRooFitUtils.so!")
    # force the loading
    import ROOT
    ROOT.RooFitUtils.ExtendedMinimizer
    return rcdir

def parseItems(expressions):
    insertions = []
    assignments = []
    for expression in expressions:
        if not expression.strip() or expression.startswith("#"):
            continue
        if "=" in expression:
            assignments.append(expression.strip())
        else:
            insertions.append(expression.strip())
    return insertions,assignments

def parseTXT(m):
    with open(m,"rt") as modelfile:
        expressions = [ e for e in modelfile ]
    return parseItems(expressions)

def parseXML(modelfile):
    from BeautifulSoup import BeautifulStoneSoup
    handler = open(modelfile).read()
    soup = BeautifulStoneSoup(handler)
    expressions = [ item["name"] for item in soup.findAll() if item.name == "item" ]
    return parseItems(expressions)

def main(args):
    import ROOT
    loadRooFitUtils()

    wsfile = ROOT.TFile.Open(args.inFileName)
    if not wsfile or not wsfile.IsOpen():
        raise IOError("unable to open file '{:s}'".format(wsfile))

    workspace = wsfile.Get(args.workspace)
    if not workspace:
        workspaces = [o.GetName() for o in wsfile.GetListOfKeys() if o.GetClassName() == "RooWorkspace"]
        if len(workspaces) == 1:
            workspace = wsfile.Get(workspaces[0])
            raise RuntimeError("unable to find workspace '{:s}' in file '{:s}', using RooWorkspace {:s} instead".format(args.workspace,args.inFileName,workspace.GetName()))
        elif len(workspaces) > 0:
            raise IOError("unable to find workspace '{:s}' in file '{:s}', available objects of type RooWorkspace are: {:s}".format(args.workspace,args.inFileName,",".join(workspaces)))
        else:
            raise IOError("unable to any workspace in file '{:s}'".format(args.inFileName))

    mc = workspace.obj(args.modelconfig)
    if not mc:
        raise IOError("unable to find model config '{:s}' in workspace '{:s}'".format(args.modelconfig,workspace.GetName()))

    orig_pois = namelist(mc.GetParametersOfInterest())
    orig_nps  = namelist(mc.GetNuisanceParameters())
    orig_globs= namelist(mc.GetGlobalObservables())

    for poi in orig_pois:
        workspace.var(poi).setConstant(True)

    pdf = mc.GetPdf()
    pdfName = pdf.GetName()

    params = []
    funcs = {}
    new_pois = []
    new_nps = []
    new_globs = []

    for m in args.modelFiles:
        if m.endswith(".xml"):
            insertions,assignments = parseXML(m)
        else:
            insertions,assignments = parseTXT(m)

        for expression in insertions:
            print("inserting new object into workspace using expression '{:s}'".format(expression))
            newobj = workspace.factory(expression);
            if not newobj:
                raise RuntimeError("error evaluating expression '{:s}'".format(expression))
            if newobj.InheritsFrom(ROOT.RooRealVar.Class()):
                new_pois.append(newobj.GetName())

        editstr = "EDIT::" + pdfName + "(" + pdfName
        for expression in assignments:
            editstr = editstr + ",\n" + expression.replace(" ","")
        editstr = editstr + ")"

        # replacement
        print("editing PDF with command '{:s}'".format(editstr))
        workspace.factory(editstr)

        newpdf = workspace.obj(pdfName)

        if not newpdf:
            raise RuntimeError("something went wrong when editing the Pdf")

    print("finalizing workspace")
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.ObjectHandling)
    newws = ROOT.RooFitUtils.makeCleanWorkspace(workspace)
    newmc = newws.obj(mc.GetName())
    newpoinames = sorted(set(new_pois))
    for poi in newpoinames:
        p = newws.var(poi)
        if p:
            p.setConstant(False)
    newmc.SetParametersOfInterest(",".join(newpoinames))
    newmc.SetNuisanceParameters(",".join(sorted(set(new_nps+orig_nps))))
    newmc.SetGlobalObservables(",".join(sorted(set(new_globs+orig_globs))))
    print("writing output to "+args.outFileName)
    newws.writeToFile(args.outFileName,True)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("reparametrize a workspace given one (or more) models")
    parser.add_argument( "--input"         , type=str,     dest="inFileName"                 , help="Input file.", required=True, metavar="path/to/workspace.root")
    parser.add_argument( "--models"        , type=str,     dest="modelFiles"                  , help="Model file.", required=True, metavar="path/to/model.txt", nargs="+")
    parser.add_argument( "--output"        , type=str,     dest="outFileName"                , help="Output file.", required=True, metavar="output.root")
    parser.add_argument( "--workspace"     , type=str,     dest="workspace"                     , help="name of the workspace." , default="combined" )
    parser.add_argument( "--modelconfig"   , type=str,     dest="modelconfig"                     , help="name of the model config." , default="ModelConfig" )

    args = parser.parse_args()
    main(args)
