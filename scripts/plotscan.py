#!/usr/bin/env python

from RooFitUtils.pgfplotter import writescans1d,writescans2d
from RooFitUtils.io import collectresults,collectpoints
from RooFitUtils.util import getPercent

def getPercentages(args,ndim):
    if args.show_percent == True or (args.show_percent == None and ndim == 2 and not args.show_sigma):
        return [0.01 * p for p in args.percent_levels ]
    elif args.show_sigma == True or (args.show_sigma == None and ndim == 1):
        return [getPercent(n) for n in args.sigma_levels]

def parseInput(inset):
    label = inset[0]
    files = inset[1:]
    if len(files) == 0:
        print("malformed argument '"+" ".join(inset)+"', needs to have format 'drawoptions file.txt [file.txt [...]]")
        exit(0)
    return label,files

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser("plot a likelihood scan")
    parser.add_argument('-i','--input',action='append',nargs="+",metavar=('drawoptions','file.txt'),help="files with input information",required=True)
    parser.add_argument('--merge-input',action='append',nargs="+",metavar=('drawoptions','file.txt'),help="more files with input information",default=[])
    parser.add_argument("--points",action='append',nargs="+",metavar=('drawoptions','file.txt'),help="text files with some additional points",required=False,default=[])
    parser.add_argument("--drawpoints",action="store_true",default=False,help="draw scan points on 1d curve")
    parser.add_argument("--atlas",type=str,help="ATLAS plot label, will enable ATLAS style if used",required=False,default=None)
    parser.add_argument("--labels",type=str,help="label(s) of the parameter axis",nargs="*",default=None)
    parser.add_argument("--poi",type=str,help="POIs to select",nargs="*",default=[])
    parser.add_argument("-o","--output",type=str,help="output file name",default="scan.tex",required=True)
    parser.add_argument("--ymax",type=float,help="y axis maximum",default=None)
    parser.add_argument("--flip-axes",action="store_true",dest="flipAxes",default=False,help="flip X and Y axes")
    parser.add_argument("--smooth",action="store_true",default=False,help="apply smoothing")
    parser.add_argument("--contour-algorithm",choices=['ROOT', 'skimage'],default="ROOT",dest="contourAlg",help="contour finding algorithm to be used")
    parser.add_argument("--npoints",type=int,default=100,help="granularity of the interpolation grid")
    parser.add_argument("--sigma-levels",type=int,nargs="+",default=[1,2])
    parser.add_argument("--percent-levels",type=int,nargs="+",default=[68,95])
    parser.add_argument( "--show-sigma", action='store_true', dest="show_sigma", help="Show sigma levels.", default=None)
    parser.add_argument( "--no-show-sigma", action='store_false', dest="show_sigma", help="Do not show sigma levels.")
    parser.add_argument( "--show-percent", action='store_true', dest="show_percent", help="Show percent levels.", default=None)
    parser.add_argument( "--no-show-percent", action='store_false', dest="show_percent", help="Do not show percent levels.")

    args = parser.parse_args()

    results = {}
    mergeresults = {}    

    for inset in args.input:
        label,files = parseInput(inset)
        collectresults(results,files,label)

    for inset in args.merge_input:
        label,files = parseInput(inset)
        collectresults(mergeresults,files,label)

    points = {}
    for inset in args.points:
        label,files = parseInput(inset)
        collectpoints(points,files,label)

    scans = results["scans"]
    fits = results["MLE"]
    
    if len(scans) == 0 and len(fits) == 0:
        print("no files found matching any of the given input paths")
        exit(1)

    pois = [tuple(p.split(",")) for p in args.poi]

    scans1d       = {k: v for k, v in scans.items() if len(k) == 1 and (len(pois)==0 or k in pois)}
    scans2d       = {k: v for k, v in scans.items() if len(k) == 2 and (len(pois)==0 or k in pois)}

    with open(args.output,"w") as outfile:
        if len(scans1d) > 0:
            if not args.labels:
                labels = list(scans1d.keys())[0]
            else:
                labels = args.labels
            writescans1d(args.atlas,labels[0],scans1d,args.output,getPercentages(args,1),args.drawpoints,args.ymax,otherscans1d=[mergeresults.get("scans",{})])
        elif len(scans2d) > 0:
            if not args.labels:
                labels = list(scans2d.keys())[0]
            else:
                labels = args.labels            
            writescans2d(args.atlas,labels,scans2d,args.output,points,args.npoints,getPercentages(args,2),contourAlg=args.contourAlg,smooth=args.smooth,flipAxes=args.flipAxes,otherscans2d=[mergeresults.get("scans",{})])
        else:
            for p in pois:
                print("no scans found for pois '"+",".join(p)+"'")


