#!/bin/env python

from RooFitUtils.pgfplotter import writescans1d,writescans2d
from RooFitUtils.io import collectresults

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser("plot a likelihood scan")
    parser.add_argument('-i','--input',action='append',nargs="+",metavar=('label','file.txt'),help="text files with the input information")
    parser.add_argument("--atlas",type=str,help="ATLAS plot label, will enable ATLAS style if used",required=False,default=None)
    parser.add_argument("--labels",type=str,help="label(s) of the parameter axis",nargs="*",default=["\\mu"])
    parser.add_argument("--poi",type=str,help="POIs to select",nargs="*",default=[])
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
        
    if len(scans) == 0 and len(results) == 0:
        print("no files found matching")
    
    pois = [tuple(p.split(",")) for p in args.poi]

    scans1d = {k: v for k, v in scans.items() if len(k) == 1 and (len(pois)==0 or k in pois)}
    scans2d = {k: v for k, v in scans.items() if len(k) == 2 and (len(pois)==0 or k in pois)}

    with open(args.output,"w") as outfile:
        if len(scans1d) > 0:
            writescans1d(args.atlas,args.labels[0],scans1d,args.output,args.ymax)
        if len(scans2d) > 0:
            writescans2d(args,scans2d)


