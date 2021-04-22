#!/usr/bin/env python

def runPruning(args):
  import ROOT, sys
  ROOT.PyConfig.IgnoreCommandLineOptions = True
  from RooFitUtils.util import loadRooFitUtils,makepctstring,getjobdims,makePOIstring,retriveObj
 # load libraries
  loadRooFitUtils()

  msrmnt = ROOT.RooFitUtils.Measurement("meas")
  if not len(args.pois):
     pois = makePOIstring(args.inws)
  else: 
   pois = args.pois

  file0 = ROOT.TFile.Open(args.fitResult[0][0])
  fitresult = file0.Get(args.fitResult[0][1])

  chesse = fitresult.covarianceMatrix().Invert()
  dim = chesse.GetNcols()
  
  pcts = ""
  if "," in args.pct:
    if len(args.pct.split(",") == len(pois.split(","): pcts = args.pct
    else: raise RuntimeError("Please provide same number of POIs and thresholds")
  else: pcts = makepctstring(pois,args.pct)
  
  order = ROOT.std.set("pair<double,std::string>")()
  for rankfile in args.orderfiles[0]:
    with open(rankfile,"r") as f:
      lines = [line.strip("\n") for line in f.readlines()]
    for line in lines:
      pair = line.split(" ")
      x = ROOT.pair("double,std::string")(float(pair[0]),pair[1])
      order.insert(x)
  
  pruneNPs = msrmnt.PruneNuisanceParameters(order, chesse, fitresult, pois, pcts, args.NPfilter)
  NPs = ""
  print("Setting {} NPs out of {} to constant ".format(len(pruneNPs),len(order)))
  for x in pruneNPs : NPs = NPs + x + ","
  NPs = NPs[0:len(NPs)-1]

  print("INFO: making pruned snapshot: "+args.snapname)
  msrmnt.MakeConstSnapshot(args.inws, fitresult, NPs, args.outws, args.snapname)
  print("INFO: saved pruned snapshots in "+ args.outws)

if __name__ == "__main__":
   from argparse import ArgumentParser
   parser = ArgumentParser("run pruning")
   arglist = []
   arglist.append(parser.add_argument( "--input"         , type=str,  dest="inws"        ,  help="input workspace name.", required=True, default=None))
   arglist.append(parser.add_argument( "--output"        , type=str,  dest="outws"       , help="output workspace name.", default=None))
   arglist.append(parser.add_argument( "--pois"          , type=str,  dest="pois"        , help="POIs to measure.", default=""))
   arglist.append(parser.add_argument( "--snapshot-name"      , type=str,  dest="snapname"    , help="name of the pruned snapshot", required=True))
   arglist.append(parser.add_argument( "--NPfilter"      , type=str,  dest="NPfilter"    , help="NPs for prune check", default=".*"))
   arglist.append(parser.add_argument( "--fitResult"     , type=str,  dest="fitResult"   , action='append',nargs="+", help="path to fit result"))
   arglist.append(parser.add_argument( "--percentages"   , type=str,  dest="pct"         , help="percentage change in poi variances",default="1"))
   arglist.append(parser.add_argument( "--order"         ,action='append',nargs="+"      , dest="orderfiles", help="files with the NPs and their ranks"))
   arglist.append(parser.add_argument( "--logsave"       , action='store_true',  dest="logsave" , help="save a log"))
   args = parser.parse_args()

   from sys import flags

   from ROOT import gSystem
   runPruning(args)

# python scripts/prune_NPs.py --input test/WS-Comb-5XS_80ifb_postFit.root --fitResult test/WS-Comb-5XS_80ifb_fitresult.root fitresult_minimizer_combData --pois r_ggF,r_VBF,r_WH,r_ZH,r_ttH --snapshot-name prune_combData_1pct --order test/order_*.txt --output test/WS-Comb-5XS_80ifb_postFit_prune.root
