#!/usr/bin/env python

def runOrdering(args):
  import ROOT,os
  ROOT.PyConfig.IgnoreCommandLineOptions = True
  from RooFitUtils.util import loadRooFitUtils,getjobdims,makePOIstring, retriveObj, getnofNPs
  msrmnt = ROOT.RooFitUtils.Measurement("meas")
 # load libraries
  loadRooFitUtils()
  pois, NPfilter = "", ".*"

  if args.pois :
     pois = args.pois

  if args.NPfilter :
     NPfilter = args.NPfilter

  else:
     pois = makePOIstring(args.inFileName)

 # fitresult = retriveObj(args.fitResult)
  print args.fitResult
  file0 = ROOT.TFile.Open(args.fitResult[0][0])
  fitresult = file0.Get(args.fitResult[0][1])

  chesse = fitresult.covarianceMatrix().Invert()

  dim = chesse.GetNcols()

  if args.writeSubmit:
    matindices = getjobdims(dim,args.jobtime) 
    with open(args.writeSubmit,"w") as f:
      for matind in matindices:
        outfilename = args.outFileName.replace(".txt","")+"_"+str(matind[0])+"_"+str(matind[1])+".txt"
        f.write("python /project/atlas/users/rahulb/project/hcomb/caf4hcomb/RooFitUtils/scripts/order_NPs.py --output "+outfilename+" --pois "+pois+" --fitResult "+args.fitResult[0][0] + " "+ args.fitResult[0][1] +" --nlo "+str(matind[0])+" --nhi "+str(matind[1])+" & \n")  
    print("INFO : wrote joblines in "+args.writeSubmit)
   
  else:
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
   arglist.append(parser.add_argument( "--pois"          , type=str,  dest="pois"         , help="POIs to measure."))
   arglist.append(parser.add_argument( "--NPfilter"      , type=str,  dest="NPfilter"     , help="NPs for prune check", default=".*"))
   arglist.append(parser.add_argument( "--fitResult"     , type=str,  dest="fitResult"    , action='append',nargs="+", help="path to fit result"))
   arglist.append(parser.add_argument( "--nlo"           , type=int,  dest="nlo"          , help="start position in NP list",default=0))
   arglist.append(parser.add_argument( "--nhi"           , type=int,  dest="nhi"          , help="end pposition in NP list",default=0))
   arglist.append(parser.add_argument( "--writeSubmit"   , type=str,  dest="writeSubmit"  , help="create a txt file for splitting the process",default=""))
   arglist.append(parser.add_argument( "--jobTime"       , type=float,dest="jobtime"      , help="ballpark-time you want each job to take in mins",default=10))
   args = parser.parse_args()
   runOrdering(args)
