#!/bin/env python

def mod(a,b):
    return a%b

def linspace(vmin,vmax,npoints):
    step = (vmax - vmin)/(npoints)
    vals = [ float(vmin + i*step) for i in range(0,npoints+1) ]
    print(vals)
    return vals

def rootcore():
    # retrieve the root core dir environment variable
    from os import getenv
    from ROOT import gROOT
    rcdir = getenv ("ROOTCOREDIR")
    if rcdir:
        gROOT.ProcessLine(".x $ROOTCOREDIR/scripts/load_packages.C")
    else:
        raise ImportError("unable to load RootCore!")
    return rcdir

def printtime(seconds):
    s = "{0:.3f}s".format(mod(seconds,60))
    if seconds < 60:
        return s
    m = "{0:d}m".format(int(mod(seconds/60,60)))
    if seconds < 3600:
        return m+" "+s
    h = "{0:d}h".format(int(seconds/3600))
    return h+" "+m+":"+s

def union(listoflists):
    s = set()
    for l in listoflists:
        s.update(l)
    return list(s)

def vec(l,t):
    import ROOT
    v = ROOT.vector(t)()
    for e in l:
        v.push_back(e)
    return v

def main(args):
    from time import time
    import ROOT
    rootcore()

    # setup verbosity
    ROOT.Log.SetReportingLevel(ROOT.LOG.FromString(args.loglevel))
    if args.loglevel == "DEBUG":
        ROOT.Math.MinimizerOptions.SetDefaultPrintLevel(1)
    else:
        ROOT.Math.MinimizerOptions.SetDefaultPrintLevel(-1)

    if ROOT.Math.MinimizerOptions.DefaultPrintLevel() < 0: ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

    # Configuration of minimizer
    ROOT.Math.MinimizerOptions.SetDefaultMinimizer(args.minimizerType, args.minimizerAlgo)
    ROOT.Math.MinimizerOptions.SetDefaultStrategy(args.defaultStrategy)

    ROOT.RooFitUtils.RooStarMomentMorphFix = args.fixCache
    ROOT.RooFitUtils.RooMultiPdfFix = args.fixMulti

    # Load the model
    model = ROOT.ExtendedModel("model", args.inFileName, args.wsName,
                               args.modelConfigName, args.dataName, args.snapshot,
                               args.binnedLikelihood, "pdf_")

    ws = model.GetWorkspace()
    mc = model.GetModelConfig()
    pdf = model.GetPdf()
    data = model.GetData()
    nuis = model.GetNuisanceParameters()
    globs = model.GetGlobalObservables()
    pois = model.GetParametersOfInterest()
    obs = model.GetObservables()

    if args.fixAllNP:          model.fixNuisanceParameters()
    if args.setInitialError:   model.setInitialErrors()
    if args.fixParameters:     model.fixNuisanceParameters(",".join(args.fixParameters))
    model.fixParametersOfInterest()

    model.profileParameters(",".join(args.profile))

    # Collect POIs
    poiset = ROOT.RooArgSet()
    poilist =  []
    for poi in args.pois:
        p = model.configureParameter(poi)
        if not p:
            raise(RuntimeError("unable to find parameter '{0:s}'".format(poi)))
        p.removeRange()
        p.setError(0.2)
        poiset.add(p)
        poilist.append(p)

    if args.makeParameterSnapshots:
        # Save the snapshots of nominal parameters
        print("Saving nominal snapshots.")
        ws.saveSnapshot("nominalGlobs", mc.GetGlobalObservables())
        ws.saveSnapshot("nominalNuis", mc.GetNuisanceParameters())
        ws.saveSnapshot("nominalPois", mc.GetParametersOfInterest())

    argelems = [ROOT.RooFit.Minimizer(args.minimizerType, args.minimizerAlgo), ROOT.RooFit.Strategy(args.defaultStrategy), 
                ROOT.ExtendedMinimizer.Eps(args.eps), ROOT.RooFit.Constrain(nuis), ROOT.RooFit.GlobalObservables(globs),
                ROOT.RooFit.NumCPU(args.numCPU, args.numThreads), ROOT.RooFit.Offset(args.offsetting), ROOT.RooFit.Optimize(args.constOpt),
                ROOT.RooFit.Precision(args.precision)]
    arglist = ROOT.RooLinkedList()
    for arg in argelems: arglist.Add(arg)
    minimizer = ROOT.ExtendedMinimizer("minimizer", pdf, data, ws,arglist)
    
    if args.fit:
        start = time()
        if not args.dummy:
            if args.findSigma:
                minimizer.minimize(ROOT.ExtendedMinimizer.Scan(poiset), ROOT.RooFit.Hesse(), ROOT.RooFit.Save())
            else:
                minimizer.minimize(ROOT.RooFit.Hesse(), ROOT.RooFit.Save())
        
        end = time()
        print("Fitting time: " + printtime(end-start))
        minNll = minimizer.GetMinNll()
        print("NLL after minimisation: "+str(minNll))
    else:
        print("no fit requested")

    parnames = None
    coords = None
    if args.scan:
        print(args.scan)
        parname = args.scan[0]
        parrange = linspace(float(args.scan[2]),float(args.scan[3]),int(args.scan[1]))
        parnames = vec([parname],"string")
        coords = vec([ vec([val],"double") for val in parrange],"vector<double>")

    if args.points != None:
        with open(args.points) as infile:
            points = [ { n.strip():float(v) for (n,v) in [x.split("=") for x in line.split(",") ]} for line in infile if len(line)>0]
        parnames = vec(sorted(union([p.keys() for p in points])),"string")
        coords = vec( [ vec( [ point[p] for p in parnames ] , "double") for point in points ], "vector<double>")

    if parnames and coords and not args.dummy:
        minimizer.scan(parnames,coords)
    else:
        print("no scan requested")

    result = minimizer.getResult()

    if result:
        if args.outFileName:
            with open(args.outFileName,'w') as out:
                if result.min.nll:
                    out.write("Minimization: minNll = {0:g}\n".format(result.min.nll))
                    for p in result.parameters:
                        out.write("{0:s} = {1:g} - {2:g} + {3:g}\n".format(p.name,p.value,abs(p.errLo),abs(p.errHi)))
                for scan in result.scans:
                    out.write((" ".join(scan.parNames)) + " nll\n")
                    for i in range(0,len(scan.nllValues)):
                        out.write((" ".join([ str(scan.parValues[i][j]) for j in range(0,len(scan.parNames)) ]))+" "+str(scan.nllValues[i])+"\n")
            print("wrote output to "+args.outFileName)
        else:
            print("no output requested")
    else:
        print("received invalid result")


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("run a fit")
    parser.add_argument( "--input"         , type=str,     dest="inFileName"                 , help="File to run over.", required=True, metavar="path/to/workspace.root")
    parser.add_argument( "--output"        , type=str,     dest="outFileName"                , help="Output file.", required=False, metavar="out.txt")
    parser.add_argument( "--poi"           , type=str,     dest="pois"                       , help="POIs to measure.", metavar="POI", nargs="+", default=[])
    parser.add_argument( "--scan"          , nargs=4,      dest="scan"                       , help="POI range to scan the Nll.", metavar=("POI","n","min","max"), default=None)
    parser.add_argument( "--points"        , type=str,     dest="points"                     , help="Points to scan the Nll at.", metavar="points.txt", default=None)
    parser.add_argument( "--snapshot"      , type=str,     dest="snapshot"                   , help="Initial snapshot.", default="nominalNuis" )
    parser.add_argument( "--makeSnapshots" , type=bool,    dest="makeParameterSnapshots"     , help="Make parameter snapshots.", default=True )
    parser.add_argument('--fit',                           dest='fit', action='store_true'   , help="Actually run the fit.", default=True )
    parser.add_argument('--no-fit',                        dest='fit', action='store_false'  , help="Do not run the fit.", default=True )
    parser.add_argument('--findSigma',                     dest='findSigma', action='store_true' , help="Search for crossings to identify the 1-sigma-band.", default=True )
    parser.add_argument('--no-findSigma',                  dest='findSigma', action='store_false', help="Do not Search for crossings.", default=True )
    parser.add_argument('--dummy',                         dest='dummy', action='store_true' , help="Perform a dummy run.", default=False )
    parser.set_defaults(fit=True)
    parser.add_argument( "--folder"        , type=str,     dest="folder"                     , help="Output folder.", default="test" )
    parser.add_argument( "--profile"       , type=str,     dest="profile"                    , help="Parameters to profile.", nargs="+", metavar="NP", default=[] )
    parser.add_argument( "--fix"           , type=str,     dest="fixParameters"              , help="Parameters to fix.", nargs="+", metavar="NP", default=[])
    parser.add_argument( "--workspace"     , type=str,     dest="wsName"                     , help="WS to grab." , default="combined" )
    parser.add_argument( "--modelconfig"   , type=str,     dest="modelConfigName"            , help="MC to load.", default="ModelConfig" )
    parser.add_argument( "--data"          , type=str,     dest="dataName"                   , help="Data to use.", default="combData" )
    parser.add_argument( "--minimizerType" , type=str,     dest="minimizerType"              , help="Minimizer type.", default="Minuit2" )
    parser.add_argument( "--minimizerAlgo" , type=str,     dest="minimizerAlgo"              , help="Minimizer algorithm.", default="Migrad" )
    parser.add_argument( "--strategy"      , type=int,     dest="defaultStrategy"            , help="Default strategy.", default=0 )
    parser.add_argument( "--numCPU"        , type=int,     dest="numCPU"                     , help="Number of CPUs.", default=1 )
    parser.add_argument( "--numThreads"    , type=int,     dest="numThreads"                 , help="Number of CPUs.", default=3 )
    parser.add_argument( "--binned"        , type=bool,    dest="binnedLikelihood"           , help="Binned likelihood.", default=True )
    parser.add_argument( "--starfix"       , type=bool,    dest="fixCache"                   , help="Fix StarMomentMorph cache.", default=True )
    parser.add_argument( "--multifix"      , type=bool,    dest="fixMulti"                   , help="Fix MultiPdf level 2.", default=True )
    parser.add_argument( "--precision"     , type=float,   dest="precision"                  , help="Precision for scan.", default=0.001 )
    parser.add_argument( "--eps"           , type=float,   dest="eps"                        , help="Convergence criterium.", default=1.0 )
    parser.add_argument( "--eigen"         , type=bool,    dest="eigendecomposition"         , help="Eigenvalues and vectors.", default=False )
    parser.add_argument( "--offset"        , type=bool,    dest="offsetting"                 , help="Offset likelihood.", default=True )
    parser.add_argument( "--initError"     , type=bool,    dest="setInitialError"            , help="Pre-set the initial error.", default=False )
    parser.add_argument( "--optimize"      , type=int,     dest="constOpt"                   , help="Optimize constant terms." , default=2)
    parser.add_argument( "--loglevel"      , type=str,     dest="loglevel"                   , help="Verbosity.", choices=["DEBUG","INFO","WARNING","ERROR"], default="DEBUG" )
    parser.add_argument( "--fixAllNP"      , type=bool,    dest="fixAllNP"                   , help="Fix all NP.", default=False )
    args = parser.parse_args()
    main(args)