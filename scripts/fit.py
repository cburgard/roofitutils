#!/usr/bin/env python

import sys

def main(args):
    if args["writeSubmit"]:
        if args["scan"] or args["refineScan"]:
            from RooFitUtils.fit import createScanJobs
            createScanJobs(args,arglist)
        if args["breakdown"]:
            from RooFitUtils.fit import createBreakdownJobs
            createBreakdownJobs(args,arglist)
        if args["impacts"]:
            from RooFitUtils.fit import createImpactJobs            
            createImpactJobs(args,arglist)                        
        exit(0)
    if args["breakdown"] or args["impacts"]:
        print("breakdown/impact calculation only supported in batch mode, please use --writeSubmit")
        exit(0)

    from sys import flags
    if not flags.interactive:
        logfile = args.get("logsave","log.txt")
        if logfile:
            from RooFitUtils.util import redirect_c_output
            fd = redirect_c_output(logfile,False)

    from RooFitUtils.fit import setup,buildModel,buildMinimizer,fit
    setup(args)
    model = buildModel(args)
    minimizer = buildMinimizer(args,model)

    if not flags.interactive:
        fit(args,model,minimizer)
        if logfile:
            from RooFitUtils.util import restore_c_output            
            restore_c_output(fd)
    else:
        print("prepared fit:")
        print("  ExtendedModel model")
        print("  ExtendedMinimizer minimizer")
        print("call 'fit(args,model,minimizer)' to run!")

    
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("run a fit")
    arglist = []
    
    basic = parser.add_argument_group("Basic Usage")
    arglist.append(basic.add_argument( "-i","--input"    , type=str,     dest="inFileName"                 , help="File to run over.", required=True, metavar="path/to/workspace.root"))
    arglist.append(basic.add_argument( "-o","--output"   , type=str,     dest="outFileName"                , help="Output file.", required=False, metavar="out.txt",default=None))
    arglist.append(basic.add_argument( "-d", "--data"    , type=str,     dest="dataName"                   , help="Data to use.", default="combData" ))
    arglist.append(basic.add_argument( "--poi"           , type=str,     dest="pois"                       , help="POIs to measure.", metavar="POI", nargs="+", default=[]))
    arglist.append(basic.add_argument( "--snapshot"      , type=str,     dest="snapshot"                   , help="Initial snapshot.", default="nominalNuis" ))
    arglist.append(basic.add_argument( "--make-snapshots", action="store_true",    dest="makeParameterSnapshots"     , help="Make parameter snapshots." ))
    arglist.append(basic.add_argument('--fit',                           dest='fit', action='store_true'   , help="Actually run the fit.", default=True ))
    arglist.append(basic.add_argument('--no-fit',                        dest='fit', action='store_false'  , help="Do not run the fit.", default=False ))
    arglist.append(basic.add_argument('--minos',                     dest='minos', action='store_true' , help="Run Minos to identify the 1-sigma-band.", default=True ))
    arglist.append(basic.add_argument('--no-minos',                  dest='minos', action='store_false', help="Do not run minos.", default=False ))    
    arglist.append(basic.add_argument('--dummy',                         dest='dummy', action='store_true' , help="Perform a dummy run.", default=False ))
    arglist.append(basic.add_argument( "--folder"        , type=str,     dest="folder"                     , help="Output folder.", default="test" ))
    arglist.append(basic.add_argument( "--profile"       , type=str,     dest="profile"                    , help="Parameters to profile.", nargs="+", metavar="NP", default=[] ))
    arglist.append(basic.add_argument( "--fix"           , type=str,     dest="fixParameters"              , help="Parameters to fix.", nargs="+", metavar="NP", default=[]))
    arglist.append(basic.add_argument( "--float"         , type=str,     dest="floatParameters"            , help="Parameters to float.", nargs="+", metavar="NP", default=[]))
    arglist.append(basic.add_argument( "--workspace"     , type=str,     dest="wsName"                     , help="WS to grab." , default="combWS" ))
    arglist.append(basic.add_argument( "--write-workspace", "--writeWorkspace", type=str,    dest="outWsName"                  , help="Filename of the output workspace." , default=None ))
    arglist.append(basic.add_argument( "--write-result", "--writeResult", action='store_true',   dest="writeResult",help="option to write the RooFitResult and HESSE matrix", default=False ))
    arglist.append(basic.add_argument( "--modelconfig", "--modelConfig"   , type=str,     dest="modelConfigName"            , help="MC to load.", default="ModelConfig" ))
    arglist.append(basic.add_argument( "--log", "--logsave"       , type=str,    dest="logsave"                    , help="saving output as log" , default=None ))
    arglist.append(basic.add_argument( "--print", action='store_true',   dest="printResult",help="option to print the fit result", default=True ))
    arglist.append(basic.add_argument( "--no-print", action='store_false',   dest="printResult",help="option to print the fit result"))        
    arglist.append(basic.add_argument( "--print-hesse", "--printHesse", action='store_true',   dest="printHesse",help="option to print Hesse matrix", default=False ))
    
    roofit = parser.add_argument_group("RooFit Options")    
    arglist.append(roofit.add_argument( "--minimizerType" , type=str,     dest="minimizerType"              , help="Minimizer type.", default="Minuit2" ))
    arglist.append(roofit.add_argument( "--minimizerAlgo" , type=str,     dest="minimizerAlgo"              , help="Minimizer algorithm.", default="Migrad" ))
    arglist.append(roofit.add_argument( "--print-level", "--printLevel" , type=int,     dest="printLevel"              , help="print level of minimizer", default=1))
    arglist.append(roofit.add_argument( "--hesse"         , action='store_true',    dest="hesse"       , help="enable HESSE", default=False ))
    arglist.append(roofit.add_argument( "--no-hesse"      , action='store_false',   dest="hesse"         , help="disable HESSE", default=False ))
    arglist.append(roofit.add_argument( "--strategy"      , type=int,     dest="defaultStrategy"            , help="Default strategy.", default=1 ))
    arglist.append(roofit.add_argument( "--num-cpu", "--numCPU"        , type=int,     dest="numCPU"                     , help="Number of CPUs.", default=1 ))
    arglist.append(roofit.add_argument( "--mp-strategy", "--mpStrategy"    , type=int,     dest="mpStrategy"                 , help="Multi-Processing strategy.", default=3 ))
    arglist.append(roofit.add_argument( "--offset"        , action='store_true',   dest="offsetting"                 , help="Offset likelihood.", default=True ))
    arglist.append(roofit.add_argument( "--eps"           , type=float,   dest="eps"                        , help="Convergence criterium.", default=0.05 ))
    arglist.append(roofit.add_argument( "--no-offset"     , action='store_false',  dest="offsetting"                 , help="Do not offset likelihood.", default=False ))
    arglist.append(roofit.add_argument( "--reuse-minimizer"        , action='store_true',   dest="reuseMinimizer"                 , help="Allow to reuse the minimizer.", default=False ))
    arglist.append(roofit.add_argument( "--no-reuse-minimizer"     , action='store_false',  dest="reuseMinimizer"                 , help="Do not allow to reuse the minimizer.", default=True ))
    arglist.append(roofit.add_argument( "--reuse-nll"        , action='store_true',   dest="reuseNll"                 , help="Allow to reuse the nll.", default=False ))
    arglist.append(roofit.add_argument( "--no-reuse-nll"     , action='store_false',  dest="reuseNll"                 , help="Do not allow to reuse the nll.", default=True ))
    arglist.append(roofit.add_argument( "--initial-error", "--initError"     , type=bool,    dest="setInitialError"            , help="Pre-set the initial error.", default=False ))
    arglist.append(roofit.add_argument( "--optimize"      , type=int,     dest="constOpt"                   , help="Optimize constant terms." , default=2))
    arglist.append(roofit.add_argument( "--log-level", "--logLevel", "--loglevel"      , type=str,     dest="loglevel"                   , help="Verbosity.", choices=["DEBUG","INFO","WARNING","ERROR"], default="ERROR" ))
    arglist.append(roofit.add_argument( "--fixAllNP"      , action='store_true',    dest="fixAllNP"                   , help="Fix all NP.", default=False ))
    arglist.append(roofit.add_argument( "--correlation-matrix", "--correlationMatrix", action='store_true',   dest="correlationMatrix",help="option to save correlation matrix", default=False ))
    arglist.append(roofit.add_argument( "--binned"        , action='store_true',    dest="binnedLikelihood"           , help="Binned likelihood.", default=True ))
    arglist.append(roofit.add_argument( "--unbinned"      , action='store_false',   dest="binnedLikelihood"           , help="Unbinned likelihood.", default=False ))

    multi = parser.add_argument_group("Scans, Rankings, and Batch-Interaction")
    arglist.append(multi.add_argument( "--scan"          , type=str,     dest="scan"                       , help="POI ranges to scan the Nll.", metavar=("POI","N","min","max"), default=None,nargs=4,action="append"))
    arglist.append(multi.add_argument( "--breakdown"     , type=str,     dest="breakdown"                  , help="nuisance parameter groups to perform breakdown on.", metavar="REGEX", default=None,nargs="+",action="append"))
    arglist.append(multi.add_argument( "--impacts"       , type=str,     dest="impacts"                    , help="nuisance parameters to calculate impacts of.", metavar="REGEX", default=None,nargs="+"))
    arglist.append(multi.add_argument( "--refine-scan"   , type=str,     dest="refineScan"                 , help="Previous scan results to refine.", default=None,nargs="+"))
    arglist.append(multi.add_argument( "--refine-scan-thresholds", type=float,     dest="refineScanThresholds", help="Likelihood thresholds to use to refine previous scan.", default=None,nargs="+"))
    arglist.append(multi.add_argument( "--points"        , type=str,     dest="points"                     , help="Points to scan the Nll at.", metavar="points.txt", default=None))
    arglist.append(multi.add_argument( "--singlepoint"   , type=str,     dest="point"                      , help="A single point to scan the Nll at.", metavar="POI_A=1,POI_B=0", default=None))
    arglist.append(multi.add_argument( "--write-submit", "--writeSubmit"   , type=str,     dest="writeSubmit"                , help="Instead of fitting, write a job definition file.", metavar="jobs.txt" ))
    arglist.append(multi.add_argument( "--submit-command", "--submitCommand" , type=str,     dest="submitCommand"              , help="Submission command to be used.", default=""))    
    arglist.append(multi.add_argument( "--job-size", "--jobSize"   , type=int,     dest="writeSubmitPoints"                , help="How many points to use per job when writing out jobs for scans.", metavar="N", default=2))
    
    
    advanced = parser.add_argument_group("Advanced Features")
    arglist.append(advanced.add_argument( "--penalty"       , type=str,     dest="penalty"                    , help="Penalty terms", metavar=("Penalty","vars"), default=None,nargs=2,action="append"))
    arglist.append(advanced.add_argument( "--penaltyfile"   , type=str,     dest="penaltyfile"                , help="Penalty terms", metavar=("Penalty"), default=None))

    arglist.append(advanced.add_argument("--find-sigma", '--findSigma',                     dest='findSigma', action='store_true' , help="Search for crossings to identify the 1-sigma-band.", default=False ))
    arglist.append(advanced.add_argument("--no-find-sigma", '--no-findSigma',                  dest='findSigma', action='store_false', help="Do not Search for crossings.", default=True ))
    arglist.append(advanced.add_argument( "--starfix"       , action='store_true',    dest="fixCache"                   , help="Fix StarMomentMorph cache.", default=True ))
    arglist.append(advanced.add_argument( "--no-starfix"    , action='store_false',   dest="fixCache"                   , help="Do not fix StarMomentMorph cache.", default=False ))
    arglist.append(advanced.add_argument( "--multifix"      , action='store_true',    dest="fixMulti"                   , help="Fix MultiPdf level 2.", default=True ))
    arglist.append(advanced.add_argument( "--no-multifix"   , action='store_false',   dest="fixMulti"                   , help="Do not fix MultiPdf level 2.", default=False ))
    arglist.append(advanced.add_argument( "--eigen"         , action='store_true',   dest="eigendecomposition"         , help="Eigenvalues and vectors.", default=False ))


    args = parser.parse_args()

    main(vars(args))
