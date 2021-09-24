#!/usr/bin/env python

def runPruning(args):
  import ROOT, sys
  ROOT.PyConfig.IgnoreCommandLineOptions = True
  from RooFitUtils.util import loadRooFitUtils,makepctstring,getjobdims
 # load libraries
  loadRooFitUtils()

  msrmnt = ROOT.RooFitUtils.Measurement("meas")

  from RooFitUtils.io import getFitResult
  pois = args.pois

  fitresult = getFitResult(*args.fitResult)   
  chesse = fitresult.covarianceMatrix().Invert()
  dim = chesse.GetNcols()

  if args.prunefromalist != None:
    with open(args.prunefromalist,"r") as f:
      lines = [line.strip("\n") for line in f.readlines()]
    pruneNPs_list = ",".join(lines)
    msrmnt.MakeConstSnapshot(args.inws, fitresult, pruneNPs_list, args.outws, args.snapname)
    return

  pcts = []
  for pct in args.pct:
    pcts.append(makepctstring(pois,pct))
  
  order = ROOT.std.set("pair<double,std::string>")()
  for rankfile in args.orderfiles[0]:
    with open(rankfile,"r") as f:
      lines = [line.strip("\n") for line in f.readlines()]
    for line in lines:
      pair = line.split(" ")
      x = ROOT.pair("double,std::string")(float(pair[0]),pair[1])
      order.insert(x)
  pruneNPs = list(msrmnt.PruneNuisanceParameters(order, chesse, fitresult, ",".join(pois), ",".join(pcts), args.NPfilter))
  print("Setting {} NPs out of {} to constant ".format(len(pruneNPs),len(order)))

  if args.outws:
    msrmnt.MakeConstSnapshot(args.inws, fitresult, ",".join(pruneNPs), args.outws, args.snapname)
  if args.outlist:
    with open(args.outlist,"wt") as outfile:
      for np in pruneNPs:
        outfile.write(np+"\n")
    

if __name__ == "__main__":
   from argparse import ArgumentParser
   parser = ArgumentParser(description="prune nuisance parameters")
   parser.add_argument( "-i","--input"         , type=str,  dest="inws"        ,  help="input workspace name.", required=True, default=None)
   parser.add_argument( "-o","--output"        , type=str,  dest="outws"       , help="output workspace name.",default=None)
   parser.add_argument( "--write-list",type=str,  dest="outlist"       , help="output list of pruned parameters.",default=None)
   parser.add_argument( "--pois"          , nargs="+", type=str,    dest="pois"        , help="POIs to measure.", required=True)
   parser.add_argument( "--percentages"   , nargs="+", type=float,  dest="pct"         , help="percentage change in poi variances",default=[1])
   parser.add_argument( "--snapshot-name"      , type=str,  dest="snapname"    , help="name of the pruned snapshot", required=True)
   parser.add_argument( "--NPfilter"      , type=str,  dest="NPfilter"    , help="NPs for prune check", default=".*")
   parser.add_argument( "--fitResult"     , type=str,  dest="fitResult"   , nargs="+", help="path to the file containing the fit result")
   parser.add_argument( "--order"         ,action='append',nargs="+"      , dest="orderfiles", help="files with the NPs and their ranks")
   parser.add_argument( "--prunefromalist",type=str,  dest="prunefromalist"       , help="txt file specifying NPs to be pruned",default=None)
   args = parser.parse_args()

   from sys import flags
   from ROOT import gSystem
   runPruning(args)
