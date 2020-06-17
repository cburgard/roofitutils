#!/usr/bin/env python

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


    with open(args.output,"w") as outfile:
    writecorrmatrix(args.atlas,parnames,results,args.output,ymax=None)
