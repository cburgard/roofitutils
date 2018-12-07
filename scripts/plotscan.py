#!/bin/env python

import os
inf = float("inf")
nan = float("nan")
from math import isnan

def writehead(stream):
    stream.write("\\documentclass{standalone}\n")
    stream.write("\\usepackage{pgfplots,tikz}\n")
    stream.write("\\usetikzlibrary{calc}\n")
    stream.write("\\usetikzlibrary{shapes.misc}\n")
    stream.write("\\tikzset{cross/.style={cross out, draw=black, minimum size=2*(#1-\pgflinewidth), inner sep=0pt, outer sep=0pt},cross/.default={3pt}}\n")
    stream.write("\\begin{document}\n")

def writefoot(stream):
    stream.write("\\end{document}\n")

def getvals(d,nllmin):
    xvals = sorted(d.keys())
    yvals = [ max(d[k] - nllmin,0) for k in xvals ]
    return zip(xvals,yvals)
        
def findminimum(points):
    from scipy.interpolate import PchipInterpolator as interpolate
    from scipy.optimize import minimize
    from numpy import array
    xvals = sorted(points.keys())
    yvals = [ points[x] for x in xvals ]
    interp = interpolate(xvals, yvals, extrapolate=True)
    minimum = minimize(lambda v:interp(v[0]), array([min(xvals)]))
    return minimum.fun

def smooth(graph):
    newgraph = []
    from math import isnan
    lastx,lasty = graph[-1]
    for x,y in graph:
        if not isnan(lastx) and not isnan(lasty):
            newx = 0.5*(x+lastx)
            newy = 0.5*(y+lasty)
            newgraph.append((newx,newy))
        lastx,lasty = x,y
    return newgraph

def findcontours(points,values):
    from numpy import array
    from math import isnan
    keys = sorted(points.keys())

    xvals = [ p[0] for p in keys ]
    yvals = [ p[1] for p in keys ]
    zvals = [ points[p] for p in keys ]

    npoints = 100
    
    from scipy.interpolate import griddata
    from numpy import mgrid
    grid_x, grid_y = mgrid[min(xvals):max(xvals):complex(npoints), min(yvals):max(yvals):complex(npoints)]

    grid_z = griddata(array(keys),array(zvals),(grid_x, grid_y), method='linear')

    nllmin = inf
    for i in range(0,npoints):
        for j in range(0,npoints):
            zval = grid_z[i][j]
            if zval < nllmin:
                nllmin = zval
                xval = grid_x[i][j]
                yval = grid_y[i][j]
    minimum = (xval,yval)
                
    from skimage import measure
    allcontours = []
    allimgcontours = []
    for v in values:
        imgcontours = measure.find_contours(grid_z, v+nllmin)
        allimgcontours.append(imgcontours)
        contours = []
        for c in imgcontours:
            realc = []
            for i,j in c:
                if isnan(i) or isnan(j): continue
                realc.append((i/npoints * (max(xvals)-min(xvals)) + min(xvals),j/npoints * (max(yvals)-min(yvals)) + min(yvals)))
            if len(realc) > 0:
                realc.append(realc[0])
                contours.append(smooth(realc))
        allcontours.append(contours)

    # Display the image and plot all contours found
#    import matplotlib.pyplot as plt
#    fig, ax = plt.subplots()
#    ax.imshow(grid_z, interpolation='nearest', cmap=plt.cm.gray)
#    for contours in allimgcontours:
#        for n, contour in enumerate(contours):
#            ax.plot(contour[:, 1], contour[:, 0], linewidth=2)
#        
#    ax.axis('image')
#    ax.set_xticks([])
#    ax.set_yticks([])
#    plt.show() 

    return allcontours,minimum
    
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

def collectresults(scans,results,files,label):
    import re
    parpat = re.compile("([a-zA-Z0-9_.]+)[ ]*=[ ]*([0-9.naife-]+)[ ]*-[ ]*([0-9.naife-]+)[ ]*\+[ ]*([0-9.naife-]+)[ ]*")
    nllpat = re.compile("Minimization:[ ]*minNll[ ]*=[ ]*([0-9.naife-]+)")
    import glob
    filenames = []
    for expression in files:
        filenames.extend(glob.glob(expression))
    if len(filenames) == 0:
        print("no results found in "+expression)
        exit(0)
    for filename in filenames:
        if os.path.isfile(filename):
            with open(filename,'r') as infile:
                lines = [ line for line in infile ]
                for lineno in range(0,len(lines)):
                    line = lines[lineno]
                    try:
                        parts = line.split()
                        nllmatch = nllpat.match(line)
                        match = parpat.match(line)
                        if nllmatch:
                            minnll = float(nllmatch.group(1))
                        elif match:
                            pname,cv,ed,eu = match.group(1).strip(),match.group(2),match.group(3),match.group(4)
                            result = (float(cv),float(ed),float(eu))
                            if not pname in results.keys():
                                results[pname] = {}
                            results[pname][label] = result
                        elif parts[-2].strip() == "nll":
                            npars = len(parts)-2
                            scanps = parts[0:npars]
                            key = tuple(scanps)
                            if key not in scans.keys():
                                scans[key] = {}
                            if label not in scans[key].keys():
                                scans[key][label] = {}
                        else:
                            pvals   = tuple([float(parts[i]) for i in range(0,len(parts)-2)])
                            nllval = float(parts[-2])
                            scans[key][label][pvals] = nllval
                    except:
                        if nllmatch:
                            print("unable to parse line '"+line.strip()+"' in file '"+filename+"', attempted to parse as nll value")
                        elif match:
                            print("unable to parse line '"+line.strip()+"' in file '"+filename+"', attempted to parse as parameter value")
                        else:      
                            print("unable to parse line '"+line.strip()+"' in file '"+filename+"'")
                        continue
                            
                for p in scans.keys():
                    if p in results.keys():
                        scans[p][label][results[p][label][0]]=minnll;

def writeATLAS(label,outfile):
    outfile.write("\\node at (rel axis cs:0,1) [anchor=south west,text width=3cm]{\\textbf{ATLAS} "+label+"};\n")

def writescans1d(atlas,par,allscans,outfilename,ymax=None):
    with open(outfilename,"w") as outfile:
        writehead(outfile)
        domain = "domain={0:f}:{1:f}".format(min([ v[0] for scan in allscans.values() for v in scan.keys() ]),max([ v[0] for scan in allscans.values() for v in scan.keys() ]))
        outfile.write("\\begin{tikzpicture}\n")
        if atlas: writeATLAS(atlas,outfile)
        outfile.write("\\begin{axis}[\n")
        outfile.write("    ymin=0,\n")
        if ymax:
            outfile.write("    ymax="+str(ymax)+",\n")
        outfile.write("    "+domain+",\n")
        outfile.write("    legend pos=north east,legend style={draw=none},\n")
        outfile.write("    xlabel=${0:s}$, ylabel=$\\Delta \\log L$\n".format(par))
        outfile.write("]\n")
        for pnamelist,scan in allscans.items():
            writescan1d(pnamelist[0],par,{k[0]:v for k,v in scan.items()},outfile,ymax)
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)

def parsedict(s):
    d = {}
    for kv in s.split(","):
        k,v = kv.split("=")
        d[k] = v
    return d

def writescans2d(args,scans2d):
    with open(args.output,"w") as outfile:    
        writehead(outfile)
        if args.atlas:
            outfile.write("\\renewcommand\\sfdefault{phv}\n")
            outfile.write("\\renewcommand\\rmdefault{phv}\n")
            outfile.write("\\renewcommand\\ttdefault{pcr}\n")            
        outfile.write("\\begin{tikzpicture}[\n")
        if args.atlas:
            outfile.write("  font={\\fontfamily{qhv}\\selectfont}\n")
        outfile.write("]\n")
        outfile.write("\\begin{axis}[clip=false,\n")
        if args.atlas:
            outfile.write("legend pos=outer north east,legend style={anchor=south east,draw=none},\n")            
            outfile.write("xticklabel={\\pgfmathprintnumber[assume math mode=true]{\\tick}},\n")
            outfile.write("yticklabel={\\pgfmathprintnumber[assume math mode=true]{\\tick}},\n")
        if len(args.labels) == 2:
            outfile.write("    xlabel="+args.labels[0]+",\n")
            outfile.write("    ylabel="+args.labels[1]+",\n")
        outfile.write("]\n")
        if args.atlas: writeATLAS(args.atlas,outfile)
    	# 1 sigma (=68.26895% CL):  2.296
    	# 2 sigma (=95.44997% CL):  6.180
        contours = {0.5*2.296:"solid",0.5*6.180:"dashed"}
        for pnamelist,scan in scans2d.items():
            for dset,points in scan.items():
                writescan2d(points,outfile,contours,parsedict(dset))
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)    
        
def writescan2d(allpoints,outfile,contourdefs,style):
    thresholds = sorted(contourdefs.keys())
    contours,minimum = findcontours(allpoints,thresholds)
    outfile.write("\\draw (axis cs:{:f},{:f}) node[cross,".format(minimum[0],minimum[1])+",color="+style.get("color","black")+"] {};\n")
    first = True
    for v,conts in zip(thresholds,contours):
        icont = 0
        for c in conts:
            outfile.write("% contour {:d} of {:f}\n".format(icont,v))
            outfile.write("\\addplot[color="+style.get("color","black")+","+contourdefs[v]+",mark=none,smooth")
            if not first: outfile.write(",forget plot")
            outfile.write("] coordinates {\n")
            for x,y in c:
                outfile.write("    ({:f},{:f})\n".format(x,y))
            outfile.write("} ;\n")
            if first and "title" in style.keys():
                outfile.write("\\addlegendentry{"+style["title"]+"};\n")
            first=False
            icont = icont+1

def writescan1d(parname,parlabel,allpoints,outfile,ymax=None):
    outfile.write("\\addplot[color=black,mark=none,smooth] coordinates {\n")
    othermin = min(allpoints.values())
    nllmin = findminimum(allpoints)
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
    outfile.write("\\addlegendentry{{${:s} = {:.3f}^{{+{:.3f}}}_{{{:.3f}}}$}}".format(parlabel,cv1,up1,down1))
    print("{:s} = {:f}, 1sigma = +{:f} -{:f}, 2sigma = +{:f} -{:f}".format(parname,cv1,up1,down1,up2,down2))

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser("plot a likelihood scan")
    parser.add_argument('-i','--input',action='append',nargs="+",metavar=('label','file.txt'),help="text files with the input information")
    parser.add_argument("--atlas",type=str,help="ATLAS plot label, will enable ATLAS style if used",required=False,default=None)
    parser.add_argument("--labels",type=str,help="label(s) of the parameter axis",nargs="*",default=["\\mu"])
    parser.add_argument("--poi",type=str,help="POIs to select",nargs="*",default=[])
    parser.add_argument("--output",type=str,help="output file name",default="scan.tex",required=True)
    parser.add_argument("--ymax",type=float,help="y axis maximum",default=None)
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
            writescans1d(args.atlas,args.label[0],scans1d,args.output,args.ymax)
        if len(scans2d) > 0:
            writescans2d(args,scans2d)


