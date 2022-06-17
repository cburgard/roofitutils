#!/usr/bin/env python

from RooFitUtils.pgfplotter import writepois
from RooFitUtils.io import collectresults,collectpoints
from RooFitUtils.util import superset

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser("plot the results of a multi-POI fit")
    parser.add_argument('-i','--input',action='append',nargs="+",metavar=('drawoptions','file.txt'),help="text files with the input information",required=True)
    parser.add_argument("-o", "--output",type=str,help="output file name",default="scan.tex",required=True)
    parser.add_argument("--atlas",type=str,help="ATLAS plot label, will enable ATLAS style if used",required=False,default=None)
    parser.add_argument("--labels",type=str,help="additional labels to be typeset",nargs="+",default=[])
    parser.add_argument("--pois",type=str,help="POIs to select",nargs="*",default=[])
    parser.add_argument("--y-spread",type=float,help="spread of POIs along the y axis",default=1)
    parser.add_argument("--x-range",type=float,help="range to show along the x axis",nargs=2,default=[None,None])
    parser.add_argument("--mark-SM",type=float,help="value to mark as SM",default=0)
    parser.add_argument('--values',action='store_true',help="print the values")
    args = parser.parse_args()

    results = {}
    for inset in args.input:
        label = inset[0]
        files = inset[1:]
        collectresults(results,files,label)

    poivals = results["MLE"]
    if args.pois:
        pois = args.pois
    else:
        pois = superset([s.keys() for s in poivals.values()])
    
    with open(args.output,"w") as outfile:
       writepois(args.atlas,pois,poivals,args.output,range=args.x_range,smval=args.mark_SM,spread=args.y_spread,printvalues=args.values,plotlabels=args.labels)
