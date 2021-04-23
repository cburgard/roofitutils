#!/usr/bin/env python

def runOrdering(args):
  import ROOT,os
  ROOT.PyConfig.IgnoreCommandLineOptions = True
  from RooFitUtils.util import loadRooFitUtils,getjobdims,getnofNPs,getfitresult
  # load libraries
  loadRooFitUtils()
  pois, NPfilter = "", ".*"

  pois = args.pois
  if args.NPfilter : NPfilter = args.NPfilter

  from RooFitUtils.io import getFitResult
  fitresult = getFitResult(*args.fitResult)
  chesse = fitresult.covarianceMatrix().Invert()
  dim = chesse.GetNcols()

  if args.writeSubmit:
    matindices = getjobdims(dim,args.jobtime) 
    with open(args.writeSubmit,"w") as f:
      for matind in matindices:
        outfilename = args.outFileName.replace(".txt","")+"_"+str(matind[0])+"_"+str(matind[1])+".txt"
        f.write("order_NPs.py --output {} --pois {} --fitResult {} {} --pos-lo {} --pos-hi {} & \n".format(outfilename,pois,args.fitResult[0][0],args.fitResult[0][1],matind[0],matind[1])) 
    print("INFO : wrote joblines in "+args.writeSubmit)
   
  else:
    msrmnt = ROOT.RooFitUtils.Measurement("meas")
    if args.nhi == 0: nlo, nhi = 0, getnofNPs(chesse,pois)
    else: nlo, nhi = args.nlo, args.nhi
    orderNP = ROOT.std.set("pair<double,string>>")()
    orderNP = msrmnt.OrderNuisanceParameters(chesse,fitresult, pois, NPfilter, nlo, nhi)

    with open(args.outFileName, "w") as f:
      for pair in orderNP:
        f.write(str(pair[0]) + " "+pair[1]+"\n")
    print("Wrote ranks of "+ str(nhi-nlo+1)+ " to "+ args.outFileName)

if __name__ == "__main__":
   from argparse import ArgumentParser
   parser = ArgumentParser("run pruning")
   arglist = []
   arglist.append(parser.add_argument( "--output"        , type=str,  dest="outFileName"  , help="Output text file. (outname.txt)", required=True, default=None))
   arglist.append(parser.add_argument( "--pois"          , type=str,  dest="pois"         , help="POIs to measure.", required=True))
   arglist.append(parser.add_argument( "--NPfilter"      , type=str,  dest="NPfilter"     , help="NPs for prune check", default=".*"))
   arglist.append(parser.add_argument( "--pos-lo"           , type=int,  dest="nlo"          , help="start position in NP list",default=0))
   arglist.append(parser.add_argument( "--pos-hi"           , type=int,  dest="nhi"          , help="end pposition in NP list",default=0))
   arglist.append(parser.add_argument( "--writeSubmit"   , type=str,  dest="writeSubmit"  , help="create a txt file for splitting the process",default=""))
   arglist.append(parser.add_argument( "--jobTime"       , type=float,dest="jobtime"      , help="ballpark-time you want each job to take in mins",default=10))
   args = parser.parse_args()
   runOrdering(args)
