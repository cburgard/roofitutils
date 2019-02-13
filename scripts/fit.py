#!/bin/env python 

import sys

def setup(args):
    # general setup and loading of modules
    import ROOT

    from RooFitUtils.util import loadRooFitUtils
    # load libraries
    loadRooFitUtils()
    
    # setup verbosity
    ROOT.RooFitUtils.Log.SetReportingLevel(ROOT.RooFitUtils.Log.FromString(args.loglevel))
    if args.loglevel == "DEBUG":
        ROOT.Math.MinimizerOptions.SetDefaultPrintLevel(1)
    else:
        ROOT.Math.MinimizerOptions.SetDefaultPrintLevel(-1)

    if ROOT.Math.MinimizerOptions.DefaultPrintLevel() < 0: ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

    # Configuration of minimizer
    ROOT.Math.MinimizerOptions.SetDefaultMinimizer(args.minimizerType, args.minimizerAlgo)
    ROOT.Math.MinimizerOptions.SetDefaultStrategy(args.defaultStrategy)

    # patches for HCombRoot
    ROOT.RooFitUtils.RooStarMomentMorphFix = args.fixCache
    ROOT.RooFitUtils.RooMultiPdfFix = args.fixMulti


def buildModel(args):
    # Load the model
    import ROOT
    model = ROOT.RooFitUtils.ExtendedModel("model", args.inFileName, args.wsName,
                                           args.modelConfigName, args.dataName, args.snapshot,
                                           args.binnedLikelihood, "pdf_")
    
    if args.fixAllNP:          model.fixNuisanceParameters()
#   if args.breakdownErrors:   model.breakdownErrors()
    if args.setInitialError:   model.setInitialErrors()
    if args.fixParameters:     model.fixParameters(",".join(args.fixParameters))

    model.fixParametersOfInterest()
    model.profileParameters(",".join(args.profile))
    
    return model

def buildMinimizer(args,model):
    import ROOT
    
    ws = model.GetWorkspace()
    mc = model.GetModelConfig()
    allparams = ROOT.RooArgSet()
    nuis = model.GetNuisanceParameters()
    pois = model.GetParametersOfInterest()
    globs = model.GetGlobalObservables()
    ROOT.RooFitUtils.addArgSet(allparams, nuis)
    ROOT.RooFitUtils.addArgSet(allparams, pois)
    ROOT.RooFitUtils.addArgSet(allparams, globs)
    obs = model.GetObservables()

    if args.makeParameterSnapshots:
       # Save the snapshots of nominal parameters
        print("Saving nominal snapshots.")
        ws.saveSnapshot("nominalGlobs", globs)
        ws.saveSnapshot("nominalNuis", nuis)
        ws.saveSnapshot("nominalPois", pois)
        ws.saveSnapshot("nominalAll", allparams)


    # Collect POIs
    pois = model.GetParametersOfInterest()
    poiset = ROOT.RooArgSet()
    if args.pois:
        poinames = args.pois
    else:
        poinames = [ p.GetName() for p in makelist(pois) ]
    for poi in poinames:
        p = model.configureParameter(poi)
        if not p:
            raise(RuntimeError("unable to find parameter '{0:s}'".format(poi)))
        p.setConstant(False)
        poiset.add(p)

    argelems = [ROOT.RooFit.Minimizer(args.minimizerType, args.minimizerAlgo), 
                ROOT.RooFit.Strategy(args.defaultStrategy), 
                ROOT.RooFitUtils.ExtendedMinimizer.Eps(args.eps), 
                ROOT.RooFitUtils.ExtendedMinimizer.ReuseMinimizer(args.reuseMinimizer), 
                ROOT.RooFitUtils.ExtendedMinimizer.ReuseNLL(args.reuseNll), 
                ROOT.RooFit.Constrain(nuis), 
                ROOT.RooFit.GlobalObservables(globs),
                ROOT.RooFit.NumCPU(args.numCPU, args.numThreads), 
                ROOT.RooFit.Offset(args.offsetting), 
                ROOT.RooFit.Optimize(args.constOpt),
                ROOT.RooFit.Precision(args.precision),
                ROOT.RooFit.Hesse(args.hesse),
                ROOT.RooFit.Save()]
    if args.findSigma:
        argelems.append(ROOT.RooFitUtils.ExtendedMinimizer.Scan(poiset)) 

    from RooFitUtils.util import nodel
    nodel(poiset)
    nodel(argelems)
    arglist = ROOT.RooLinkedList()
    for arg in argelems: arglist.Add(arg)

    minimizer = ROOT.RooFitUtils.ExtendedMinimizer("minimizer", model,arglist)
    return minimizer


def fit(args,model,minimizer):
    from time import time
    from RooFitUtils.util import parsePoint,timestamp,linspace,vec,mkdir
    from RooFitUtils.io import writeResult
    import ROOT

    # Collect POIs
    pois = model.GetParametersOfInterest()
    if args.pois:
        poinames = args.pois
    else:
        poinames = [ p.GetName() for p in makelist(pois) ]
    for poi in poinames:
        p = model.configureParameter(poi)
        if not p:
            raise(RuntimeError("unable to find parameter '{0:s}'".format(poi)))
        p.setConstant(False)

    if args.fit:
        start = time()
        if not args.dummy:
            minimizer.minimize()
        
        end = time()
        print("Fitting time: " + timestamp(end-start))
        minNll = minimizer.GetMinNll()
        print("NLL after minimisation: "+str(minNll))

        if args.makeParameterSnapshots:
            ws = model.GetWorkspace()
            mc = model.GetModelConfig()
            allparams = ROOT.RooArgSet()
            nuis = model.GetNuisanceParameters()
            allparams.add(nuis)
            globs = model.GetGlobalObservables()
            allparams.add(globs)
            pois = model.GetParametersOfInterest()
            allparams.add(pois)
            obs = model.GetObservables()
            # Save the snapshots of nominal parameters
            print("Saving minimum snapshots.")
            ws.saveSnapshot("minimumGlobs", globs)
            ws.saveSnapshot("minimumNuis", nuis)
            ws.saveSnapshot("minimumPois", pois)
            ws.saveSnapshot("minimumAll", allparams)
    else:
        print("no fit requested")

    parnames = None
    coords = None
    if args.scan:
        val = args.scan.split()
        parname = val[0]
        parrange = linspace(float(val[2]),float(val[3]),int(val[1]))
        parnames = vec([parname],"string")
        coords = vec([ vec([val],"double") for val in parrange],"vector<double>")

    if args.points != None:
        with open(args.points) as infile:
            points = [ parsePoint(line) for line in infile if len(line)>0 ]
            parnames = vec(sorted(union([p.keys() for p in points])),"string")
            coords = vec( [ vec( [ point[p] for p in parnames ] , "double") for point in points ], "vector<double>")

    if args.point != None:
        point = parsePoint(args.point)
        parnames = vec(sorted(point.keys()),"string")
        coords = vec( [ vec( [ point[p] for p in parnames ] , "double") ], "vector<double>")

    if parnames and coords and not args.dummy:
        minimizer.scan(parnames,coords)
    else:
        print("no scan requested")

    result = minimizer.getResult()

    if result:
        if args.outFileName:
            import os
            outpath,outfile  = os.path.split(args.outFileName)
            mkdir(outpath)
            with open(args.outFileName,'w') as out:
                writeResult(out,result,args.hesse)
            print("wrote output to "+args.outFileName)
        else:
            print("no output requested")
    else:
        print("received invalid result")

    if args.outWsName:
        ws = model.GetWorkspace()        
        ws.writeToFile(args.outWsName)


def createScanJobs(args,arglist,pointsPerJob):
    from RooFitUtils.util import stringify,makepoint,reconstructCall,generateCoordsDict,mkdir
    options = reconstructCall(args,arglist,["scan","findSigma","writeSubmit","refineScan"])
    import sys
    name = sys.argv[0]
    coords = generateCoordsDict(args.scan)
    idx = 0
    import os
    outpath,outfile = os.path.split(args.writeSubmit)
    pointspath = outpath+"/coords_0.txt"
    mkdir(outpath)
    from RooFitUtils.util import clearfile
    clearfile(args.writeSubmit)
    clearfile(pointspath)
 
    idx = 0
    ipoints = 0
    for coord in coords:
        ipoints = ipoints + 1
        if  ipoints % pointsPerJob == 0:  
	    pointspath =outpath+"/coords" +"_"+str(idx)+".txt"
	    clearfile(pointspath)
            with open(args.writeSubmit,"a") as jobs:
	        options[" --no-findSigma --points"]=pointspath
                if args.outFileName:
                    options["--output"]=args.outFileName+".part"+str(idx)
                cmd = " ".join([k+" "+stringify(v) for k,v in options.items()])
		if not os.path.exists(args.outFileName+".part"+str(idx)):
                    jobs.write(name+" "+cmd+"\n")
	    idx = idx + 1
        with open(pointspath,"a") as coordlist:
          point = makepoint(coord)
          coordlist.write(point+"\n")
   
    print("wrote "+args.writeSubmit)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("run a fit")
    arglist = []
    arglist.append(parser.add_argument( "--input"         , type=str,     dest="inFileName"                 , help="File to run over.", required=True, metavar="path/to/workspace.root"))
    arglist.append(parser.add_argument( "--output"        , type=str,     dest="outFileName"                , help="Output file.", required=False, metavar="out.txt"))
    arglist.append(parser.add_argument( "--poi"           , type=str,     dest="pois"                       , help="POIs to measure.", metavar="POI", nargs="+", default=[]))
    arglist.append(parser.add_argument( "--scan"          , type=str,     dest="scan"                       , help="POI ranges to scan the Nll.", metavar=("POI","N","min","max"), default=None,nargs=4,action="append"))
    arglist.append(parser.add_argument( "--points"        , type=str,     dest="points"                     , help="Points to scan the Nll at.", metavar="points.txt", default=None))
    arglist.append(parser.add_argument( "--singlepoint"   , type=str,     dest="point"                      , help="A single point to scan the Nll at.", metavar="POI_A=1,POI_B=0", default=None))
    arglist.append(parser.add_argument( "--snapshot"      , type=str,     dest="snapshot"                   , help="Initial snapshot.", default="nominalNuis" ))
    arglist.append(parser.add_argument( "--make-snapshots", action="store_true",    dest="makeParameterSnapshots"     , help="Make parameter snapshots." ))
    arglist.append(parser.add_argument('--fit',                           dest='fit', action='store_true'   , help="Actually run the fit.", default=True ))
    arglist.append(parser.add_argument('--no-fit',                        dest='fit', action='store_false'  , help="Do not run the fit.", default=True ))
    arglist.append(parser.add_argument('--findSigma',                     dest='findSigma', action='store_true' , help="Search for crossings to identify the 1-sigma-band.", default=True ))
    arglist.append(parser.add_argument('--no-findSigma',                  dest='findSigma', action='store_false', help="Do not Search for crossings.", default=True ))
    arglist.append(parser.add_argument('--dummy',                         dest='dummy', action='store_true' , help="Perform a dummy run.", default=False ))
    parser.set_defaults(fit=True)
    arglist.append(parser.add_argument( "--folder"        , type=str,     dest="folder"                     , help="Output folder.", default="test" ))
    arglist.append(parser.add_argument( "--profile"       , type=str,     dest="profile"                    , help="Parameters to profile.", nargs="+", metavar="NP", default=[] ))
    arglist.append(parser.add_argument( "--fix"           , type=str,     dest="fixParameters"              , help="Parameters to fix.", nargs="+", metavar="NP", default=[]))
    arglist.append(parser.add_argument( "--workspace"     , type=str,     dest="wsName"                     , help="WS to grab." , default="combined" ))
    arglist.append(parser.add_argument( "--write-workspace", type=str,    dest="outWsName"                  , help="Filename of the output workspace." , default=None ))
    arglist.append(parser.add_argument( "--modelconfig"   , type=str,     dest="modelConfigName"            , help="MC to load.", default="ModelConfig" ))
    arglist.append(parser.add_argument( "--data"          , type=str,     dest="dataName"                   , help="Data to use.", default="combData" ))
    arglist.append(parser.add_argument( "--minimizerType" , type=str,     dest="minimizerType"              , help="Minimizer type.", default="Minuit2" ))
    arglist.append(parser.add_argument( "--minimizerAlgo" , type=str,     dest="minimizerAlgo"              , help="Minimizer algorithm.", default="Migrad" ))
    arglist.append(parser.add_argument( "--hesse"         , action='store_true',    dest="hesse"       , help="enable HESSE", default=False ))
    arglist.append(parser.add_argument( "--no-hesse"      , action='store_false',   dest="hesse"	     , help="disable HESSE", default=False ))
    arglist.append(parser.add_argument( "--strategy"      , type=int,     dest="defaultStrategy"            , help="Default strategy.", default=0 ))
    arglist.append(parser.add_argument( "--numCPU"        , type=int,     dest="numCPU"                     , help="Number of CPUs.", default=1 ))
    arglist.append(parser.add_argument( "--numThreads"    , type=int,     dest="numThreads"                 , help="Number of CPU Threads.", default=3 ))
    arglist.append(parser.add_argument( "--writeSubmit"   , type=str,     dest="writeSubmit"                , help="Instead of fitting, write a job definition file.", metavar="jobs.txt" ))
    arglist.append(parser.add_argument( "--binned"        , action='store_true',    dest="binnedLikelihood"           , help="Binned likelihood.", default=True ))
    arglist.append(parser.add_argument( "--unbinned"      , action='store_false',   dest="binnedLikelihood"           , help="Unbinned likelihood.", default=False ))
    arglist.append(parser.add_argument( "--starfix"       , action='store_true',    dest="fixCache"                   , help="Fix StarMomentMorph cache.", default=True ))
    arglist.append(parser.add_argument( "--no-starfix"    , action='store_false',   dest="fixCache"                   , help="Do not fix StarMomentMorph cache.", default=False ))
    arglist.append(parser.add_argument( "--multifix"      , action='store_true',    dest="fixMulti"                   , help="Fix MultiPdf level 2.", default=True ))
    arglist.append(parser.add_argument( "--no-multifix"   , action='store_false',   dest="fixMulti"                   , help="Do not fix MultiPdf level 2.", default=False ))
    arglist.append(parser.add_argument( "--precision"     , type=float,   dest="precision"                  , help="Precision for scan.", default=0.001 ))
    arglist.append(parser.add_argument( "--eps"           , type=float,   dest="eps"                        , help="Convergence criterium.", default=1.0 ))
    arglist.append(parser.add_argument( "--eigen"         , action='store_true',   dest="eigendecomposition"         , help="Eigenvalues and vectors.", default=False ))
    arglist.append(parser.add_argument( "--offset"        , action='store_true',   dest="offsetting"                 , help="Offset likelihood.", default=True ))
    arglist.append(parser.add_argument( "--no-offset"     , action='store_false',  dest="offsetting"                 , help="Do not offset likelihood.", default=False ))
    arglist.append(parser.add_argument( "--reuse-minimizer"        , action='store_true',   dest="reuseMinimizer"                 , help="Allow to reuse the minimizer.", default=False ))
    arglist.append(parser.add_argument( "--no-reuse-minimizer"     , action='store_false',  dest="reuseMinimizer"                 , help="Do not allow to reuse the minimizer.", default=True ))
    arglist.append(parser.add_argument( "--reuse-nll"        , action='store_true',   dest="reuseNll"                 , help="Allow to reuse the nll.", default=False ))
    arglist.append(parser.add_argument( "--no-reuse-nll"     , action='store_false',  dest="reuseNll"                 , help="Do not allow to reuse the nll.", default=True ))
    arglist.append(parser.add_argument( "--initError"     , type=bool,    dest="setInitialError"            , help="Pre-set the initial error.", default=False ))
    arglist.append(parser.add_argument( "--optimize"      , type=int,     dest="constOpt"                   , help="Optimize constant terms." , default=2))
    arglist.append(parser.add_argument( "--loglevel"      , type=str,     dest="loglevel"                   , help="Verbosity.", choices=["DEBUG","INFO","WARNING","ERROR"], default="DEBUG" ))
    arglist.append(parser.add_argument( "--logsave"       , type=bool,    dest="logsave"                    , help="saving output as log" , default=False ))
    arglist.append(parser.add_argument( "--fixAllNP"      , action='store_true',    dest="fixAllNP"                   , help="Fix all NP.", default=False ))

    args = parser.parse_args()

    if args.writeSubmit and args.scan:
        createScanJobs(args,arglist,40)
        exit(0)

    from sys import flags
    if not flags.interactive:
        if args.logsave:
            log_file = open(args.outFileName+".log","w")
            sys.stdout = log_file

    setup(args)
    model = buildModel(args)
    minimizer = buildMinimizer(args,model)
 
    if not flags.interactive:
        fit(args,model,minimizer)

        if args.logsave:
            sys.stdout = sys.__stdout__ 
            log_file.close()
    else:
        print("prepared fit:")
        print("  ExtendedModel model")
        print("  ExtendedMinimizer minimizer")
        print("call 'fit(args,model,minimizer)' to run!")
