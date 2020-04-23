#!/bin/env python

from RooFitUtils.pgfplotter import writecorrmatrix
from RooFitUtils.io import collectcorrelationmatrix

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser("plot a correlation matrix")
    parser.add_argument('-i','--input',action='append',nargs="+",metavar=('drawoptions','file.txt'),help="text files with the input information",required=True)
    parser.add_argument("--atlas",type=str,help="ATLAS plot label, will enable ATLAS style if used",required=False,default=None)
    parser.add_argument("--labels",type=str,help="label(s) of the parameter axis",nargs="*",default=["\\mu"])
    parser.add_argument("--output",type=str,help="output file name",default="scan.tex",required=True)
    parser.add_argument("--parFilter",type=str,help="filter for parameters in correlation matrix",default="",required=False)
    args = parser.parse_args()

    parnames = []
    results  = []
    for inset in args.input:
        label = inset[0]
        files = inset[1:][0]
        collectcorrelationmatrix(parnames,results,files,args.parFilter,label)

   # print(parnames)
   # print(results)

  #  print(results)
    #print(scans)
#    points = {}
#    for inset in args.points:
#        label = inset[0]
#        files = inset[1:]
#        collectpoints(points,files,label)

#    if len(scans) == 0 and len(results) == 0:
#        print("no files found matching")

#    pois = [tuple(p.split(",")) for p in args.poi]

#    scans1d = {k: v for k, v in scans.items() if len(k) == 1 and (len(pois)==0 or k in pois)}
#    scans2d = {k: v for k, v in scans.items() if len(k) == 2 and (len(pois)==0 or k in pois)}

    with open(args.output,"w") as outfile:
    writecorrmatrix(args.atlas,parnames,results,args.output,ymax=None)

#        if len(scans1d) > 0:
#            writescans1d(args.atlas,args.labels[0],scans1d,args.output,args.ymax)
#        if len(scans2d) > 0:
#            writescans2d(args,scans2d,points)
