#!/usr/bin/env python

from RooFitUtils.pgfplotter import writepulls
from RooFitUtils.io import collectresults,collectpoints

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser("plot a likelihood scan")
    parser.add_argument('--range',nargs=2,default=[-5,5],help="range to plot",type=float)    
    parser.add_argument('-i','--input',action='append',nargs="+",metavar=('drawoptions','file.txt'),help="text files with the input information",required=True)
    parser.add_argument("--atlas",type=str,help="ATLAS plot label, will enable ATLAS style if used",required=False,default=None)
    parser.add_argument('-o',"--output",type=str,help="output file name",default="pulls.tex",required=True)
    args = parser.parse_args()

    scans = {}
    results = {}
    for inset in args.input:
        label = inset[0]
        files = inset[1:]
        collectresults(scans,results,files,label)

    with open(args.output,"wt") as outfile:
        writepulls(args,results,outfile)


