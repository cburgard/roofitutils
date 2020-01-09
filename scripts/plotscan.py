#!/bin/env python

from RooFitUtils.pgfplotter import writescans1d,writemergescans2d,writescans2d
from RooFitUtils.io import collectresults,collectpoints

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
    parser.add_argument('-i','--input',action='append',nargs="+",metavar=('drawoptions','file.txt'),help="text files with input information",required=True)
    parser.add_argument('--mergeinput',action='append',nargs="+",metavar=('drawoptions','file.txt'),help="text files with input information")
    parser.add_argument("--points",action='append',nargs="+",metavar=('drawoptions','file.txt'),help="text files with some additional points",required=False,default=[])
    parser.add_argument("--atlas",type=str,help="ATLAS plot label, will enable ATLAS style if used",required=False,default=None)
    parser.add_argument("--labels",type=str,help="label(s) of the parameter axis",nargs="*",default=["\\mu"])
    parser.add_argument("--poi",type=str,help="POIs to select",nargs="*",default=[])
    parser.add_argument("--output",type=str,help="output file name",default="scan.tex",required=True)
    parser.add_argument("--ymax",type=float,help="y axis maximum",default=None)
    parser.add_argument("--flip-axes",action="store_true",dest="flipAxes",default=False,help="flip X and Y axes")
    parser.add_argument("--smooth",action="store_true",default=False,help="apply smoothing")
    parser.add_argument("--contour-algorithm",choices=['ROOT', 'skimage'],default="ROOT",dest="contourAlg",help="contour finding algorithm to be used")
    parser.add_argument("--npoints",type=int,default=100,help="granularity of the interpolation grid")

    args = parser.parse_args()

    scans = {}
    results = {}

    scans_merge = {}
    results_merge = {}  

    mergeinp = ""
    for inset in args.input:
        label,files = parseInput(inset)
        collectresults(scans,results,files,label)

    if args.mergeinput != None:
        for inset in args.mergeinput:
            label,files = parseInput(inset)
            collectresults(scans_merge,results_merge,files_merge,label_merge)
    
    points = {}
    for inset in args.points:
        label,files = parseInput(inset)
        collectpoints(points,files,label)
        
    if len(scans) == 0 and len(results) == 0:
        print("no files found matching any of the given input paths")
    
    pois = [tuple(p.split(",")) for p in args.poi]

    scans1d       = {k: v for k, v in scans.items() if len(k) == 1 and (len(pois)==0 or k in pois)}
    scans2d       = {k: v for k, v in scans.items() if len(k) == 2 and (len(pois)==0 or k in pois)}
    scans2d_merge = {k: v for k, v in scans_merge.items() if len(k) == 2 and (len(pois)==0 or k in pois)}

    with open(args.output,"w") as outfile:
        if len(scans1d) > 0:
            writescans1d(args.atlas,args.labels[0],scans1d,args.output,args.ymax)
        elif len(scans2d) > 0:
            writescans2d(args,scans2d,points,args.npoints)
        elif len(scans2d_merge) > 0 and mergeinp != "" :
            writemergescans2d(args,scans2d,scans2d_merge,points,args.npoints)
        else:
            for p in pois:
                print("no scans found for pois '"+",".join(p)+"'")

