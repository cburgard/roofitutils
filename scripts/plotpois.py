#!/bin/env python

from RooFitUtils.pgfplotter import writepois
from RooFitUtils.io import collectresults,collectpoints

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser("plot a likelihood scan")
    parser.add_argument('-i','--input',action='append',nargs="+",metavar=('drawoptions','file.txt'),help="text files with the input information",required=True)
    parser.add_argument("--points",action='append',nargs="+",metavar=('drawoptions','file.txt'),help="text files with some additional points",required=False,default=[])
    parser.add_argument("--atlas",type=str,help="ATLAS plot label, will enable ATLAS style if used",required=False,default=None)
    parser.add_argument("--labels",type=str,help="label(s) of the parameter axis",nargs="*",default=["\\mu"])
    parser.add_argument("--pois",type=str,help="POIs to select",nargs="*",default=[])
    parser.add_argument("--output",type=str,help="output file name",default="scan.tex",required=True)
    parser.add_argument("--ymax",type=float,help="y axis maximum",default=None)
    parser.add_argument("--flip-axes",action="store_true",dest="flipAxes",default=False)
    parser.add_argument("--smooth",action="store_true",default=False)
    args = parser.parse_args()

    scans = {}
    results = {}
    for inset in args.input:
        label = inset[0]
        files = inset[1:]
        collectresults(scans,results,files,label)

    print(args.pois)
    poivals = { k: v for k, v in results.items() if k in args.pois}

    print(poivals)

    with open(args.output,"w") as outfile:
       writepois(args.atlas,args.labels[0],poivals,args.output,args.ymax)
#
