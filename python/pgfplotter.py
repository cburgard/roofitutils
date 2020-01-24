#!/bin/evn python

from RooFitUtils.util import getThreshold
thresholdColors = ["blue","green","yellow","orange","red"]
thresholdStyles = ["solid","dashed","loosely dashed","dotted","loosely dotted"]

def writehead(stream):
    stream.write("\\documentclass{standalone}\n")
    stream.write("\\usepackage{scalerel}\n")
    stream.write("\\usepackage{pgfplots,tikz}\n")
    stream.write("\\usetikzlibrary{calc}\n")
    stream.write("\\usepackage[scaled=1]{helvet}\n")
    stream.write("\\usepackage[helvet]{sfmath}\n")
    stream.write("\\usepackage{amsmath,latexsym}\n")
    stream.write("\\usetikzlibrary{shapes.misc}\n")
    stream.write("\\tikzset{cross/.style={cross out, draw=black, minimum size=2*(#1-\pgflinewidth), inner sep=0pt, outer sep=0pt},cross/.default={3pt}}\n")
    stream.write("\\catcode`_=\\active\n")
    stream.write("\\newcommand_[1]{\\ensuremath{\\sb{\\scriptscriptstyle #1}}}\n")
    stream.write("\\begin{document}\n")
    stream.write("\\definecolor{myyellow}{rgb}{0.96,0.742,0.29}\n")
    stream.write("\\definecolor{myblue}{rgb}{0.1,0.32,0.738}\n")

def writefoot(stream):
    stream.write("\\end{document}\n")

def writeATLAS(label,outfile):
    outfile.write("\\node at (rel axis cs:0,1) [anchor=south west,text width=3cm]{\\textbf{ATLAS} "+label+"};\n")

def writepois(atlas,par,allscans,outfilename,ymax=None):
   """write a pois/NP ranking styled plot to a pgfplots tex file"""
   with open(outfilename,"w") as outfile:
       writehead(outfile)
       allvals = [v[0] for pois in allscans.values() for v in pois.values()]

def writecorrelations(stream,allcorrs):
    for x in range(0,len(allcorrs)):
        stream.write(allcorrs[x]+"\n")

def concat(strlist):
    string = ""
    for x in range(0,len(strlist)):
      string = string +","+ strlist[x]
    return string[1:len(string)]

def writecorrmatrix(atlas,parslist,allcorrs,outfilename,ymax=None):
    """write a correlation matrix to a pgfplots tex file"""
    paraxisy = parslist
    paraxisx = parslist[::-1]
    with open(outfilename,"w") as outfile:
        writehead(outfile)
        outfile.write("\\begin{tikzpicture}\n")
        outfile.write("\\begin{axis}[\n")
        outfile.write("    colormap={bluewhiteyellow}{color=(myyellow) color=(white) color=(myblue)},\n")
        outfile.write("    clip = false,\n")
        outfile.write("    colorbar,\n")
        outfile.write("    colormap name={bluewhiteyellow},\n")
        outfile.write("    x=1em,\n")
        outfile.write("    y=1em,\n")
        outfile.write("    xtick=data,\n")
        outfile.write("    ytick=data,\n")
        outfile.write("    ymin="+ paraxisy[0] + ",\n")
        outfile.write("    ymax="+ paraxisy[len(parslist)-1] +",\n")     
        outfile.write("    xmin="+ paraxisx[0] +",\n")
        outfile.write("    xmax="+ paraxisx[len(parslist)-1] +",\n")
        outfile.write("    enlarge x limits={abs=0.5em},\n")
        outfile.write("    enlarge y limits={abs=0.5em},\n")
        outfile.write("    point meta min=-1,\n")
        outfile.write("    point meta max=+1,\n")
        outfile.write("    grid=both,\n")
        outfile.write("    major grid style={draw=none},\n")
        outfile.write("    minor tick num=1,\n")
        outfile.write("    symbolic x coords={"+ concat(paraxisx) + "},\n")
        outfile.write("    symbolic y coords={"+ concat(paraxisy) + "},\n")
        outfile.write("    axis on top,")
        outfile.write("       x tick label style={rotate=90},\n")
        outfile.write("    tick style={draw=none}\n ]\n")
        outfile.write("\\addplot [matrix plot*,point meta=explicit,mesh/cols="+str(len(parslist))+"] table [meta=correlations] {\n")
        outfile.write("x  y  correlations\n")
        writecorrelations(outfile,allcorrs) 
        outfile.write("};\n")
        outfile.write("\\node (atlas) [above right, font={\\fontfamily{phv}\\fontseries{b}\selectfont}] at (rel axis cs:0,1) {ATLAS};\n")
        outfile.write("\\node [anchor=west] at (atlas.east) {Internal};\n")
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)

def writepois(atlas,allpois,outfilename,ymax=None):
    """write a POI plot to a pgfplots tex file"""
    from RooFitUtils.io import texprep
    with open(outfilename,"w") as outfile:
        writehead(outfile)
        outfile.write("\\begin{tikzpicture}\n")
        outfile.write("\\begin{axis}[\n")
        outfile.write("    width = 0.8\\textwidth,\n")
        outfile.write("    height = 1\\textwidth, \n")
        outfile.write("    xlabel = Best fit value \n")
        outfile.write("    clip = false,\n")
        outfile.write("    ymin=-1,\n")
        outfile.write("    ymax= "+str(4.75 + len(allpois)*1.25)+ ",\n")     
        outfile.write("    xmin= -1.75,\n")
        outfile.write("    xmax=  2,\n")
        outfile.write("    minor tick num=4,\n")
        outfile.write("    ytick style={draw=none},\n")
        outfile.write("    yticklabels=\empty\n")
        outfile.write("]")
        outfile.write("\\addplot+ [color=black, sharp plot,only marks,error bars/.cd,x dir=both, x explicit]\n coordinates{\n")
        count = 0
        for x in allpois:
            outfile.write("("+str(allpois[x][0])+","+ str(0.25+1.25*count)+") += ("+str(allpois[x][1]) +",0) -= ("+str(allpois[x][2])+",0) \n")
            count = count + 1
        outfile.write("};\n")
        count = 0
        for x in allpois:
            outfile.write("\\node at ( axis cs:-1.75,"+ str(1.00+1.25*count)+") [anchor = north west]{"+texprep(x)+"};\n")
            count = count + 1
        outfile.write("\\node (atlas) [above right, font={\\fontfamily{phv}\\fontseries{b}\selectfont}] at (rel axis cs:0,1) {ATLAS};\n")
        outfile.write("\\node [anchor=west] at (atlas.east) {Internal};\n")
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)
        print("wrote "+outfilename)

def writecorrmatrix(atlas,parslist,allcorrs,outfilename,ymax=None):
    """write a correlation matrix to a pgfplots tex file"""
    paraxisy = parslist
    paraxisx = parslist[::-1]
    with open(outfilename,"w") as outfile:
        writehead(outfile)
        outfile.write("\\begin{tikzpicture}\n")
        outfile.write("\\begin{axis}[\n")
        outfile.write("    colormap={bluewhiteyellow}{color=(myyellow) color=(white) color=(myblue)},\n")
        outfile.write("    clip = false,\n")
        outfile.write("    colorbar,\n")
        outfile.write("    colormap name={bluewhiteyellow},\n")
        outfile.write("    x=1em,\n")
        outfile.write("    y=1em,\n")
        outfile.write("    xtick=data,\n")
        outfile.write("    ytick=data,\n")
        outfile.write("    ymin="+ paraxisy[0] + ",\n")
        outfile.write("    ymax="+ paraxisy[len(parslist)-1] +",\n")     
        outfile.write("    xmin="+ paraxisx[0] +",\n")
        outfile.write("    xmax="+ paraxisx[len(parslist)-1] +",\n")
        outfile.write("    enlarge x limits={abs=0.5em},\n")
        outfile.write("    enlarge y limits={abs=0.5em},\n")
        outfile.write("    point meta min=-1,\n")
        outfile.write("    point meta max=+1,\n")
        outfile.write("    grid=both,\n")
        outfile.write("    major grid style={draw=none},\n")
        outfile.write("    minor tick num=1,\n")
        outfile.write("    symbolic x coords={"+ concat(paraxisx) + "},\n")
        outfile.write("    symbolic y coords={"+ concat(paraxisy) + "},\n")
        outfile.write("    axis on top,")
        outfile.write("       x tick label style={rotate=90},\n")
        outfile.write("    tick style={draw=none}\n ]\n")
        outfile.write("\\addplot [matrix plot*,point meta=explicit,mesh/cols="+str(len(parslist))+"] table [meta=correlations] {\n")
        outfile.write("x  y  correlations\n")
        writecorrelations(outfile,allcorrs) 
        outfile.write("};\n")
        outfile.write("\\node (atlas) [above right, font={\\fontfamily{phv}\\fontseries{b}\selectfont}] at (rel axis cs:0,1) {ATLAS};\n")
        outfile.write("\\node [anchor=west] at (atlas.east) {Internal};\n")
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)
        print("wrote "+outfilename)

def writescans1d(atlas,par,allscans,outfilename,percent_thresholds=None,drawpoints=False,ymax=None):
    """write a bunch of 1d scans to a pgfplots tex file"""
    with open(outfilename,"w") as outfile:
        writehead(outfile)
        allvals = [ v[0] for curves in allscans.values() for scan in curves.values() for v in scan.keys()]
        domain = "domain={0:f}:{1:f}".format(min(allvals),max(allvals))
        outfile.write("\\begin{tikzpicture}\n")
        if atlas: writeATLAS(atlas,outfile)
        outfile.write("\\begin{axis}[\n")
        outfile.write("    ymin=0,\n")
        if ymax:
            outfile.write("    ymax="+str(ymax)+",\n")
        outfile.write("    "+domain+",\n")
        outfile.write("    legend pos=north east,legend style={draw=none},\n")
        outfile.write("    xlabel=${0:s}$, ylabel=$-2\\ \\ln \\Lambda$,\n".format(par))
        outfile.write("    every axis x label/.style={at={(axis description cs:1.0,-0.1)},anchor=north east},\n")
	outfile.write("    every axis y label/.style={at={(axis description cs:-0.1,1.0)},rotate=90,anchor=south east},\n")
        outfile.write("    xmin={0:f},xmax={1:f}\n".format(min(allvals),max(allvals)))
        outfile.write("]\n")
        for pnamelist,curve in allscans.items():
            for options,scan in curve.items():
                print("writing scan for "+pnamelist[0])
                writescan1d(pnamelist[0],par,{k[0]:v for k,v in scan.items()},outfile,percent_thresholds,drawpoints,ymax)
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)
        print("wrote "+outfilename)

def writescan1d(parname,parlabel,allpoints,outfile,percent_thresholds,drawpoints=False,ymax=None):
    """write a single 1d sncan to a pgfplots tex file"""
    from math import isnan
    from RooFitUtils.interpolate import findminimum,findcrossings
    from RooFitUtils.util import graphrescale
    """obtaining the - 2 log Lambda(x), where Lambda = LL(x)/LL(x0)  """
    nllmin = findminimum(allpoints)
    points = graphrescale(allpoints,nllmin,2)
    outfile.write("\\addplot[color=black,very thick,mark=none,smooth] coordinates {\n")
    for x,y in points:  outfile.write("    ({0:f},{1:f})\n".format(x,y))
    outfile.write("};\n")
    if drawpoints:
      outfile.write("\\addplot[draw=none,mark=x] coordinates {\n")
      for x,y in points:  outfile.write("    ({0:f},{1:f})\n".format(x,y))
      outfile.write("};\n")
    outfile.write("\\addplot[gray,densly dashed,thick] {1};\n")
    outfile.write("\\addplot[gray,densly dashed,thick] {4};\n") 

    if percent_thresholds:
        xmin,ymax = points[0]
    	xmax = xmin
        for x,y in points:
	    xmin = min(x,xmin)
   	    xmax = max(x,xmax)
       	    ymax = max(y,ymax)

        for i in range(0,len(percent_thresholds)):
  	    t = 1*getThreshold(percent_thresholds[i],1)
	    if t < ymax:
	        cv,down,up = findcrossings(points,t)
	        if cv-down > xmin and not isnan(down):
		    outfile.write("\\draw["+thresholdColors[i]+"] (axis cs:"+str(cv-down)+",0) -- (axis cs:"+str(cv-down)+","+str(t)+");\n")
	        if cv+up < xmax and not isnan(up):
		    outfile.write("\\draw["+thresholdColors[i]+"] (axis cs:"+str(cv+up  )+",0) -- (axis cs:"+str(cv+up  )+","+str(t)+");\n")
	        if i == 0:
	  	    s = "{:s} = {:f}".format(parname,cv)
		    outfile.write("\\addlegendentry{{${:s} = {:.1f}^{{+{:.3f}}}_{{-{:.3f}}}$}}".format(parlabel,cv,abs(up),abs(down)))
	        s = s + ", {:.3f}% CL = +{:f} -{:f}".format(100*percent_thresholds[i],abs(up),abs(down))
        print(s)

def writescans2d(args,scans2d,extrapoints,npoints,percent_thresholds):
    """write a bunch of 2d scans to a pgfplots tex file"""
    from RooFitUtils.util import parsedict
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
        outfile.write("\\begin{axis}[clip=false,minor tick num=4,\n")
        if args.atlas:
            outfile.write("legend pos=outer north east,legend style={anchor=south east,draw=none},\n")            
            outfile.write("xticklabel={\\pgfmathprintnumber[assume math mode=true]{\\tick}},\n")
            outfile.write("yticklabel={\\pgfmathprintnumber[assume math mode=true]{\\tick}},\n")
        if len(args.labels) == 2:
            if args.flipAxes:
                outfile.write("    ylabel=$"+args.labels[0]+"$,\n")
                outfile.write("    xlabel=$"+args.labels[1]+"$,\n")
            else:
                outfile.write("    xlabel=$"+args.labels[0]+"$,\n")
                outfile.write("    ylabel=$"+args.labels[1]+"$,\n")
        outfile.write("]\n")
        if args.atlas: writeATLAS(args.atlas,outfile)

        for pnamelist,scan in scans2d.items():
            for drawopts,points in scan.items():
                writescan2d(args,points,outfile,percent_thresholds,parsedict(drawopts),npoints)
        for drawopts,points in extrapoints.items():
            writepoints2d(args,points,outfile,parsedict(drawopts))
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)    
        print("wrote "+args.output)

def writemergescans2d(args,scans2d,scans2d_merge,extrapoints,npoints,percent_thresholds):
    """write a bunch of 2d scans to a pgfplots tex file"""
    from RooFitUtils.util import parsedict
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
        outfile.write("\\begin{axis}[clip=false,minor tick num=4,\n")
        if args.atlas:
            outfile.write("legend pos=outer north east,legend style={anchor=south east,draw=none},\n")
            outfile.write("xticklabel={\\pgfmathprintnumber[assume math mode=true]{\\tick}},\n")
            outfile.write("yticklabel={\\pgfmathprintnumber[assume math mode=true]{\\tick}},\n")
        if len(args.labels) == 2:
            if args.flipAxes:
                outfile.write("    ylabel="+args.labels[0]+",\n")
                outfile.write("    xlabel="+args.labels[1]+",\n")
            else:
                outfile.write("    xlabel="+args.labels[0]+",\n")
                outfile.write("    ylabel="+args.labels[1]+",\n")
        outfile.write("]\n")
        if args.atlas: writeATLAS(args.atlas,outfile)

        for pnamelist,scan in scans2d.items():
            for drawopts,points in scan.items():
                for pnamelist_merge,scan_merge in scans2d_merge.items():
                    for drawopts,points_merge in scan_merge.items():
                         writemergescan2d(args,points,points_merge,outfile,percent_thresholds,parsedict(drawopts))
        for drawopts,points in extrapoints.items():
            writepoints2d(args,points,outfile,parsedict(drawopts))
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)
        print("wrote "+args.output)


def writepoints2d(args,points,outfile,style):
    outfile.write("\\addplot[mark=x,mark options={scale=.5},only marks,draw="+style.get("color","black")+"] coordinates {\n")
    i = 0
    red = int(style.get("reduce",0))
    from random import randint    
    for point in points:
        i = i + 1
        if red > 0 and randint(0, red) != 0: continue
        keys = sorted(point.keys())
        if args.flipAxes:
            outfile.write("    ({:f},{:f})\n".format(point[keys[1]],point[keys[0]]))
        else:
            outfile.write("    ({:f},{:f})\n".format(point[keys[0]],point[keys[1]]))            
    outfile.write("};\n")

                          
def writescan2d(args,allpoints,outfile,percent_thresholds,style,npoints):
    """write a single 2d scan to a pgfplots tex file"""
    from RooFitUtils.interpolate import findcontours
    thresholds = [ getThreshold(p,2) for p in percent_thresholds ]
    contours,minimum = findcontours(allpoints,thresholds,args.smooth,npoints,args.contourAlg)
    outfile.write("\\draw (axis cs:")
    if args.flipAxes:
        outfile.write("{:f},{:f}".format(minimum[1],minimum[0]))
    else:
        outfile.write("{:f},{:f}".format(minimum[0],minimum[1]))
    outfile.write(") node[cross,color="+style.get("color","black")+"] {};\n")
    
    first = True
    i = 0
    for v,conts in zip(thresholds,contours):
        icont = 0
        for c in conts:
            outfile.write("% contour {:d} of {:f}\n".format(icont,v))
            outfile.write("\\addplot[color="+style.get("color","black")+","+thresholdStyles[i]+",mark=none,smooth")
            if not first: outfile.write(",forget plot")
            outfile.write("] coordinates {\n")
            for x,y in c:
                if args.flipAxes:
                    outfile.write("    ({:f},{:f})\n".format(y,x))                    
                else:
                    outfile.write("    ({:f},{:f})\n".format(x,y))
            outfile.write("} ;\n")
            if first and "title" in style.keys():
                outfile.write("\\addlegendentry{"+style["title"]+"};\n")
            first=False
            icont = icont+1
        i = i+1

def writemergescan2d(args,allpoints1,allpoints2,outfile,percent_thresholds,style):
     """merge two 2d scan to obtain the minimum envelope and write to pgfplots tex file"""
     from RooFitUtils.interpolate import findmergecontours
     thresholds = [ getThreshold(p,2) for p in percent_thresholds ]
     contours,minimum = findmergecontours(allpoints1,allpoints2,thresholds,args.smooth,npoints)
     outfile.write("\\draw (axis cs:")
     if args.flipAxes:
         outfile.write("{:f},{:f}".format(minimum[1],minimum[0]))
     else:
         outfile.write("{:f},{:f}".format(minimum[0],minimum[1]))
     outfile.write(") node[cross,color="+style.get("color","black")+"] {};\n")
 
     first = True
     for v,conts in zip(thresholds,contours):
         icont = 0
         for c in conts:
             outfile.write("% contour {:d} of {:f}\n".format(icont,v))
             outfile.write("\\addplot[color="+style.get("color","black")+","+contourdefs[v]+",mark=none,smooth")
             if not first: outfile.write(",forget plot")
             outfile.write("] coordinates {\n")
             for x,y in c:
                 if args.flipAxes:
                     outfile.write("    ({:f},{:f})\n".format(y,x))
                 else:
                     outfile.write("    ({:f},{:f})\n".format(x,y))
             outfile.write("} ;\n")
             if first and "title" in style.keys():
                 outfile.write("\\addlegendentry{"+style["title"]+"};\n")
             first=False

