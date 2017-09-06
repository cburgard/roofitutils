#!/bin/env python

import os
inf = float("inf")
nan = float("nan")
from math import isnan

def writehead(stream):
    stream.write("\\documentclass{standalone}\n")
    stream.write("\\usepackage{pgfplots,tikz}\n")
    stream.write("\\usetikzlibrary{calc}\n")
    stream.write("\\begin{document}\n")

def writefoot(stream):
    stream.write("\\end{document}\n")

def getvals(d,nllmin):
    xvals = sorted(d.keys())
    yvals = [ max(d[k] - nllmin,0) for k in xvals ]
    return zip(xvals,yvals)
        

def findcrossings(points,nllval):
    from scipy.interpolate import PchipInterpolator as interpolate
    xvals = [ x for x,y in points ]
    yvals = [ y for x,y in points ]
    i0 = 0
    ymin = inf
    for i in range(0,len(xvals)):
        if yvals[i] < ymin:
            ymin = yvals[i]
            i0 = i
    x0 = xvals[i0]
    xl,xr = min(xvals),max(xvals)
    r = abs(xr - xl)
    xll = xl - 0.1*r
    xrr = xr + 0.1*r
    interp = interpolate(xvals, yvals, extrapolate=True)
    from scipy.optimize import ridder as solve
    up = nan
    down = nan
    if   interp(xl ) > nllval:      down = solve(lambda x : interp(x) - nllval,xl ,x0)
    elif interp(xll) > nllval:      down = solve(lambda x : interp(x) - nllval,xll,x0)
    if   interp(xr ) > nllval:      up   = solve(lambda x : interp(x) - nllval,x0,xr )
    elif interp(xrr) > nllval:      up   = solve(lambda x : interp(x) - nllval,x0,xrr)
    r = (x0,x0-down,up-x0)
    return r

def collectscan(paths):
    scan = {}
    for path in paths:
        with open(path,'r') as infile:
            lines = [ line for line in infile ]
            for line in lines[1:]:
                vals = line.split()
                scan[float(vals[0])] = float(vals[1])
    return scan

def writescan(par,allpoints,outfilename,ymax=None):
    nllmin = min(allpoints.values())
    nllmax = max(allpoints.values())
    with open(outfilename,"w") as outfile:
        writehead(outfile)
        domain = "domain={0:f}:{1:f}".format(min(allpoints.keys()),max(allpoints.keys()))
        outfile.write("\\begin{tikzpicture}\n")
        outfile.write("\\begin{axis}[\n")
        outfile.write("    ymin=0,\n")
        if ymax:
            outfile.write("    ymax="+str(ymax)+",\n")
        outfile.write("    "+domain+",\n")
        outfile.write("    xlabel=${0:s}$, ylabel=$\\Delta \\log L$\n".format(par))
        outfile.write("]\n")
        outfile.write("\\addplot[color=black,mark=none,smooth] coordinates {\n")
        points =  getvals(allpoints,nllmin)
        for x,y in points:  outfile.write("    ({0:f},{1:f})\n".format(x,y))
        outfile.write("};\n")
        outfile.write("\\addplot[draw=none,mark=x] coordinates {\n")
        for x,y in points:  outfile.write("    ({0:f},{1:f})\n".format(x,y))
        outfile.write("};\n")
        outfile.write("\\addplot[red] {0.5};\n")

        cv1,down1,up1 = findcrossings(points,0.5)
        if not isnan(down1):
            outfile.write("\\draw[green] (axis cs:"+str(cv1-down1)+",0) -- (axis cs:"+str(cv1-down1)+",0.5);\n")
        if not isnan(up1):
            outfile.write("\\draw[green] (axis cs:"+str(cv1+up1)+",0) -- (axis cs:"+str(cv1+up1)+",0.5);\n")
        cv2,down2,up2 = findcrossings(points,2)
        if not isnan(down2):
            outfile.write("\\draw[blue] (axis cs:"+str(cv2-down2)+",0) -- (axis cs:"+str(cv2-down2)+",2.);\n")
        if not isnan(up2):
            outfile.write("\\draw[blue] (axis cs:"+str(cv2+up2)+",0) -- (axis cs:"+str(cv2+up2)+",2.);\n")
        print("POI = {:f}, 1sigma = +{:f} -{:f}, 2sigma = +{:f} -{:f}".format(cv1,up1,down1,up2,down2))

        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser("plot a likelihood scan")
    parser.add_argument("input",type=str,help="text files with the input information",nargs="+")
    parser.add_argument("--label",type=str,help="label of the x-axis",default="\\mu")
    parser.add_argument("--output",type=str,help="output file name",default="scan.tex")
    parser.add_argument("--ymax",type=float,help="y axis maximum",default=None)
    args = parser.parse_args()

    scan = collectscan(args.input)
    writescan(args.label,scan,args.output,args.ymax)