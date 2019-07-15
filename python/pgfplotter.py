#!/bin/evn python

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
	print allvals

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
    print parslist
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
        outfile.write("	   x tick label style={rotate=90},\n")
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

def writescans1d(atlas,par,allscans,outfilename,ymax=None):
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
        outfile.write("    xlabel=${0:s}$, ylabel=$\\Delta \\log L$\n".format(par))
        outfile.write("]\n")
        for pnamelist,curve in allscans.items():
            for options,scan in curve.items():
                print("writing scan for "+pnamelist[0])
                writescan1d(pnamelist[0],par,{k[0]:v for k,v in scan.items()},outfile,ymax)
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

def writecorrmatrix(atlas,parslist,allcorrs,outfilename,ymax=None):
    """write a correlation matrix to a pgfplots tex file"""
    paraxisy = parslist
    paraxisx = parslist[::-1]
    print parslist
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
        outfile.write("	   x tick label style={rotate=90},\n")
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

def writescans1d(atlas,par,allscans,outfilename,ymax=None):
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
        outfile.write("    xlabel=${0:s}$, ylabel=$\\Delta \\log L$\n".format(par))
        outfile.write("]\n")
        for pnamelist,curve in allscans.items():
            for options,scan in curve.items():
                print("writing scan for "+pnamelist[0])
                writescan1d(pnamelist[0],par,{k[0]:v for k,v in scan.items()},outfile,ymax)
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)

def writescan1d(parname,parlabel,allpoints,outfile,ymax=None):
    """write a single 1d scan to a pgfplots tex file"""
    from math import isnan
    from RooFitUtils.interpolate import findminimum,findcrossings
    from RooFitUtils.util import graphoffset
    othermin = min(allpoints.values())
    nllmin = findminimum(allpoints)
    points =  graphoffset(allpoints,nllmin)
    
    outfile.write("\\addplot[color=black,mark=none,smooth] coordinates {\n")
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

def writescans2d(args,scans2d,extrapoints):
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
    	# 1 sigma (=68.26895% CL):  2.296
    	# 2 sigma (=95.44997% CL):  6.180
	# 95% CL		 :  2.28
        # 68% CL 		 :  5.99
#        contours = {0.5*2.296:"solid",0.5*6.180:"dashed"}
        contours = {0.5*2.28:"solid",0.5*5.99:"dashed"}
        for pnamelist,scan in scans2d.items():
            for drawopts,points in scan.items():
                writescan2d(args,points,outfile,contours,parsedict(drawopts))
        for drawopts,points in extrapoints.items():
            writepoints2d(args,points,outfile,parsedict(drawopts))
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)    

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
                          
def writescan2d(args,allpoints,outfile,contourdefs,style):
    """write a single 2d scan to a pgfplots tex file"""
    from RooFitUtils.interpolate import findcontours
    thresholds = sorted(contourdefs.keys())
    contours,minimum = findcontours(allpoints,thresholds,args.smooth)
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
            icont = icont+1

