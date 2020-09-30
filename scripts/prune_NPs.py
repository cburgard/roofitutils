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

  fitresult = retriveObj(args.fitResult)
  chesse = fitresult.covarianceMatrix().Invert()

  dim = chesse.GetNcols()
  
  pcts = {1}
  pctstrings = {}
  if args.pct: percentages.append(args.pct)
  else: 
    for i in pcts : pctstrings[str(i)] = makepctstring(pois,i) 

  
  order = ROOT.std.set("pair<double,std::string>")()
  for rankfile in args.orderfiles[0]:
    with open(rankfile,"r") as f:
      lines = [line.strip("\n") for line in f.readlines()]
    for line in lines:
      pair = line.split(" ")
      x = ROOT.pair("double,std::string")(float(pair[0]),pair[1])
      order.insert(x)
  
  pruneNPs, count = {}, 0
  for pct in pcts:
   pruneNPs[str(pct)] = msrmnt.PruneNuisanceParameters(order, chesse, fitresult, pois, pctstrings[str(pct)], args.NPfilter)
   NPs = ""
   for x in pruneNPs[str(pct)] : NPs = NPs + x + ","
   NPs = NPs[0:len(NPs)-1]
   if len(args.snapname):
     snapname = args.snapname 
   else:
     snapname = "prune_"+str(pct)+"pct"

   print("INFO: making pruned snapshot: "+snapname)
   if count == 0: msrmnt.MakeConstSnapshot(args.inws, fitresult, NPs, args.outws, snapname)
   else : msrmnt.MakeConstSnapshot(args.outws, fitresult, NPs, args.outws, snapname)
   count = count + 1
  print("INFO: saved pruned snapshots in "+ args.outws)

if __name__ == "__main__":
   from argparse import ArgumentParser
   parser = ArgumentParser("run pruning")
   arglist = []
   arglist.append(parser.add_argument( "--input"         , type=str,  dest="inws"        ,  help="input workspace name.", required=True, default=None))
   arglist.append(parser.add_argument( "--outputWS"      , type=str,  dest="outws"       , help="output workspace name.", default=None))
   arglist.append(parser.add_argument( "--pois"          , type=str,  dest="pois"        , help="POIs to measure.", default=""))
   arglist.append(parser.add_argument( "--snapshot"      , type=str,  dest="snapname"    , help="name of the pruned snapshot", default=""))
   arglist.append(parser.add_argument( "--NPfilter"      , type=str,  dest="NPfilter"    , help="NPs for prune check", default=".*"))
   arglist.append(parser.add_argument( "--fitResult"     , type=str,  dest="fitResult"   , help="path to fit result"))
   arglist.append(parser.add_argument( "--percentages"   , type=str,  dest="pct"         , help="percentage change in poi variances"))
   arglist.append(parser.add_argument( "--order"         ,action='append',nargs="+"      , dest="orderfiles", help="files with the NPs and their ranks"))
   arglist.append(parser.add_argument( "--logsave"       , action='store_true',  dest="logsave" , help="save a log"))
   args = parser.parse_args()

   from sys import flags

   from ROOT import gSystem
   runPruning(args)
#python scripts/prune_NP.py --input prune_test/WS-Comb-5XS_80ifb_postFit.root --pois r_ggF,r_VBF,r_WH,r_ZH,r_ttH --data combData --hesse prune_test/5XS_hesse.root --fitResult prune_test/5XS_fitresult.root --order prune_test/order_NPs*.txt --outputWS prune_test/WS-Comb-5XS_80ifb_postFit_prune.root --logsave


