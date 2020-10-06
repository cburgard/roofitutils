
from RooFitUtils.util import getThreshold
thresholdColors = ["blue","green","yellow","orange","red"]
thresholdStyles = ["solid","dashed","loosely dashed","dotted","loosely dotted"]

def writehead(stream):
    stream.write("\\documentclass{standalone}\n")
    stream.write("\\usepackage{scalerel}\n")
    stream.write("\\usepackage{pgfplots,tikz}\n")
    stream.write("\\pgfplotsset{compat=1.15}\n")
    stream.write("\\usetikzlibrary{calc}\n")
    stream.write("\\usepackage[scaled=1]{helvet}\n")
    stream.write("\\usepackage[helvet]{sfmath}\n")
    stream.write("\\usepackage{amsmath,latexsym}\n")
    stream.write("\\usetikzlibrary{shapes.misc}\n")
    stream.write("\\tikzset{cross/.style={cross out, draw=black, minimum size=2*(#1-\pgflinewidth), inner sep=0pt, outer sep=0pt},cross/.default={3pt}}\n")
#    stream.write("\\catcode`_=\\active\n")
#    stream.write("\\newcommand_[1]{\\ensuremath{\\sb{\\scriptscriptstyle #1}}}\n")
    stream.write("\\begin{document}\n")
    stream.write("\\pgfplotsset{scaled x ticks=false}\n")
    stream.write("\\begin{tikzpicture}[font={\\fontfamily{qhv}\\selectfont}]\n")
    stream.write("\\renewcommand\\sfdefault{phv}\n")
    stream.write("\\renewcommand\\rmdefault{phv}\n")
    stream.write("\\renewcommand\\ttdefault{pcr}\n")

    stream.write("\\definecolor{myyellow}{rgb}{0.96,0.742,0.29}\n")
    stream.write("\\definecolor{myblue}{rgb}{0.1,0.32,0.738}\n")

def writefoot(stream):
    stream.write("\\end{document}\n")

def writeATLAS(label,outfile):
    outfile.write("\\node at (rel axis cs:0.025,0.975) [anchor= north west,text width=4cm]{\\textit{\\textbf{ATLAS}} Internal};\n")
    outfile.write("\\node at (rel axis cs:0.025,0.920) [anchor= north west,text width=4cm]{\\scriptsize{$\\sqrt{s}=$13 TeV, 139 fb$^{\\scriptsize{-1}}$}};\n")
    outfile.write("\\node at (rel axis cs:0.025,0.865) [anchor= north west,text width=4cm]{\\scriptsize{$m_{\\scalebox{.9}{$\\scriptstyle H$}}=$ 125.09 GeV, $|y_{\\scalebox{.9}{$\\scriptstyle H$}}|$ $<$ 2.5}};\n")


def concat(strlist):
    string = ""
    for x in range(0,len(strlist)):
        string = string +","+ strlist[x]
    return string[1:len(string)] # remove comma which is the first character

def writepois(atlas,poinames,allpois,outfilename,ymax=None):
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
	for ipoi in poinames:
	    x = ipoi
	    t = allpois[x]['color=black']
            outfile.write("("+str(t[0])+","+ str(0.25+1.25*count)+") += ("+str(t[1]) +",0) -= ("+str(t[2])+",0) \n")
            count = count + 1
        outfile.write("};\n")
        count = 0
	for x in poinames:
            outfile.write("\\node at ( axis cs:-1.75,"+ str(1.00+1.25*count)+") [anchor = north west]{"+texprep(x)+"};\n")
            count = count + 1
        outfile.write("\\node (atlas) [above right, font={\\fontfamily{phv}\\fontseries{b}\selectfont}] at (rel axis cs:0,1) {ATLAS};\n")
        outfile.write("\\node [anchor=west] at (atlas.east) {Internal};\n")
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)
        print("wrote "+outfilename)

def makelabels(listofstrings):
    x = [ i.replace("_","\_") for i in listofstrings]
    return x

def writematrix(atlas,xcoords,ycoords,allvalues,outfilename,minval=None,maxval=None,rotatelabels=90):
    """write a correlation matrix to a pgfplots tex file"""
    if len(ycoords) != len(allvalues):
        print(len(allvalues),len(ycoords))
        raise RuntimeError("incompatible lengths in y")
    if len(xcoords) != len(allvalues[0]):
        raise RuntimeError("incompatible lengths in x")    
    xlabels = [ i.replace("_","\_") for i in xcoords ]
    ylabels = [ i.replace("_","\_") for i in ycoords ]
    with open(outfilename,"w") as outfile:
        writehead(outfile)
        outfile.write("\\begin{tikzpicture}[\n")
        if atlas:
            outfile.write("  font={\\fontfamily{qhv}\\selectfont}\n")        
        outfile.write("]\n")
        outfile.write("\\begin{axis}[\n")
        outfile.write("    colormap={bluewhiteyellow}{color=(myyellow) color=(white) color=(myblue)},\n")
        outfile.write("    clip = false,\n")
        outfile.write("    colorbar,\n")
        outfile.write("    colormap name={bluewhiteyellow},\n")
        outfile.write("    x=3em,\n")
        outfile.write("    y=3em,\n")
        outfile.write("    xtick=data,\n")
        outfile.write("    ytick=data,\n")
        outfile.write("    ymin="+ ycoords[0] + ",\n")
        outfile.write("    ymax="+ ycoords[len(ycoords)-1] +",\n")
        outfile.write("    xmin="+ xcoords[0] +",\n")
        outfile.write("    xmax="+ xcoords[len(xcoords)-1] +",\n")
        outfile.write("    enlarge x limits={abs=1.5em},\n")
        outfile.write("    enlarge y limits={abs=1.5em},\n")
        if minval:
            outfile.write("    point meta min="+str(minval)+",\n")
        if maxval:
            outfile.write("    point meta max="+str(maxval)+",\n")
        outfile.write("    grid=both,\n")
        outfile.write("    major grid style={draw=none},\n")
        outfile.write("    minor tick num=1,\n")
        outfile.write("    symbolic x coords={"+ concat(xcoords) + "},\n")
        outfile.write("    symbolic y coords={"+ concat(ycoords) + "},\n")
        outfile.write("    xticklabels={"+ concat(xlabels) + "},\n") # no typo
        outfile.write("    yticklabels={"+ concat(ylabels) + "},\n") # no typo
        outfile.write("    axis on top,")
        outfile.write("       x tick label style={rotate="+str(rotatelabels)+"},\n")
        outfile.write("    tick style={draw=none}\n ]\n")
        outfile.write("\\addplot [matrix plot*,point meta=explicit,mesh/cols="+str(len(ycoords))+",mesh/rows="+str(len(xcoords))+"] table [meta=correlations] {\n")
        outfile.write("x  y  correlations\n")
        for x in range(0,len(xcoords)):
            for y in range(0,len(ycoords)):
                outfile.write(" "+xcoords[x]+" "+ycoords[y]+" "+str(allvalues[y][x])+"\n")
        outfile.write("};\n")
        for x in range(0,len(xcoords)):
            for y in range(0,len(ycoords)):
                if abs(allvalues[y][x]) > 0:
                    outfile.write("\\node at (axis cs:"+str(xcoords[x])+","+str(ycoords[y])+"){"+str(allvalues[y][x])+"};\n")
        if atlas:
            outfile.write("\\node (atlas) [scale=2,above right] at (rel axis cs:0,1) {\\bfseries\\itshape ATLAS};\n")
            outfile.write("\\node [scale=2,anchor=west] at (atlas.east) {"+atlas+"};\n")
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)
        print("wrote "+outfilename)

def writecorrmatrix(atlas,parslist,allcorrs,outfilename,ymax=None):
    writematrix(atlas,parslist,parslist,allcorrs,outfilename,ymax)
        
def writescans1d(atlas,par,allscans,outfilename,percent_thresholds=None,drawpoints=False,ymax=None):
    """write a bunch of 1d scans to a pgfplots tex file"""
    with open(outfilename,"w") as outfile:
        writehead(outfile)
        allvals = [ v[0] for curves in allscans.values() for scan in curves.values() for v in scan.keys()]
        domain = "domain={0:f}:{1:f}".format(min(allvals),max(allvals))
      #  outfile.write("\\begin{tikzpicture}\n")
        outfile.write("\\begin{axis}[\n")
        outfile.write("    ymin=0,ymax=10.\n")
        if ymax:
            outfile.write("    ymax="+str(ymax)+",\n")
        outfile.write("    "+domain+",\n")
        outfile.write("    legend pos=north east,legend style={draw=none},\n")
        outfile.write("    xlabel=${0:s}$, ylabel=$-2\\ \\ln \\Lambda$,\n".format(par))
        outfile.write("    every axis x label/.style={at={(axis description cs:1.0,-0.1)},anchor=north east},\n")
        outfile.write("    every axis y label/.style={at={(axis description cs:-0.1,1.0)},rotate=90,anchor=south east},\n")
        outfile.write("    xmin={0:f},xmax={1:f}\n".format(min(allvals),max(allvals)))
        outfile.write("]\n")
        if atlas: writeATLAS(atlas,outfile)
        for pnamelist,curve in allscans.items():
            for options,scan in curve.items():
                print("writing scan for "+pnamelist[0])
                writescan1d(pnamelist[0],par,{k[0]:v for k,v in scan.items()},options,outfile,percent_thresholds,drawpoints,ymax)
        outfile.write("\\addplot[gray,densely dashed,thick] {1};\n")
        outfile.write("\\addplot[gray,densely dashed,thick] {4};\n") 
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)
        print("wrote "+outfilename)

def writescan1d(parname,parlabel,allpoints,options,outfile,percent_thresholds,drawpoints=False,ymax=None):
    """write a single 1d sncan to a pgfplots tex file"""
    from math import isnan
    from RooFitUtils.interpolate import findminimum,findcrossings
    from RooFitUtils.util import graphrescale
    """obtaining the - 2 log Lambda(x), where Lambda = LL(x)/LL(x0)  """
    nllmin = findminimum(allpoints)
    points = graphrescale(allpoints,nllmin,2)
    outfile.write("\\addplot["+options+",very thick,mark=none,smooth] coordinates {\n")
    for x,y in points:  outfile.write("    ({0:f},{1:f})\n".format(x,y))
    outfile.write("};\n")
    if drawpoints:
        outfile.write("\\addplot[draw=none,mark=x] coordinates {\n")
        for x,y in points:  outfile.write("    ({0:f},{1:f})\n".format(x,y))
        outfile.write("};\n")
    outfile.write("\\addplot[gray,densely dashed,thick] {1};\n")
    outfile.write("\\addplot[gray,densely dashed,thick] {4};\n")

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
               # if cv-down > xmin and not isnan(down):
               #     outfile.write("\\draw["+thresholdColors[i]+"] (axis cs:"+str(cv-down)+",0) -- (axis cs:"+str(cv-down)+","+str(t)+");\n")
                #if cv+up < xmax and not isnan(up):
                #    outfile.write("\\draw["+thresholdColors[i]+"] (axis cs:"+str(cv+up  )+",0) -- (axis cs:"+str(cv+up  )+","+str(t)+");\n")
                if i == 0:
                    s = "{:s} = {:f}".format(parname,cv)
                    outfile.write("\\addlegendentry{{${:s} = {:.5f}^{{+{:.5f}}}_{{-{:.5f}}}$}}".format(parlabel,cv,abs(up),abs(down)))
                s = s + ", {:.3f}% CL = +{:f} -{:f}".format(100*percent_thresholds[i],abs(up),abs(down))
       # print(s)

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
    thresholds = [ getThreshold(p,2)*0.5 for p in percent_thresholds ]
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
     thresholds = [ getThreshold(p,2)*0.5 for p in percent_thresholds ]
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

def getColorDefStringLaTeX(name,color):
    c = ROOT.gROOT.GetColor(color)
    return getColorDefStringLaTeX(name,c)

def getColorDefStringLaTeX(name,color):
  if not color: return ""
  r,g,b = 0.,0.,0.
  color.GetRGB(r,g,b);
  return "\\definecolor{"+name+"}{rgb}{"+str(r)+","+str(g)+","+str(b)+"}"


def writepulls(args,results,outfile):
    writehead(outfile)
    xunit = 3
    yunit = .5
    outfile.write("\\begin{tikzpicture}[x="+str(xunit)+"cm,y="+str(yunit)+"cm,%\n")
    outfile.write("  lbl/.style={scale=1,anchor=west},%\n")
    outfile.write("  axlbl/.style={scale=0.5,anchor=center},%\n")
    outfile.write("  pull/.style={{|[scale=1]}-{|[scale=1]}},%\n")
    outfile.write("  dot/.style={circle,fill,inner sep=1pt},\n")
    outfile.write("  every node/.append style={font=\\sffamily}\n")
    outfile.write("]\n")
    outfile.write("\\pgfdeclarelayer{background}\\pgfsetlayers{background,main}\n")
    allpars = sorted(results.keys())
    npar = len(allpars)
    for np in range(0,npar):
        text = allpars[np]
        outfile.write("\\node[lbl] at ("+str(args.range[1]+1)+","+str(np-npar)+") {")
        if "{" in text and not "$" in text:
            outfile.write("\\ensuremath{"+text+"}")
        else:
            outfile.write(text.replace("_","\\_"))
        outfile.write("};\n")
    for np in range(0,npar):
        res = results[allpars[np]]
        ires = 0
        for style,(cv,edn,eup) in res.items():
            ires = ires+1
            offset = float(ires)/(len(res)+1) - 0.5
            if cv+abs(eup) > args.range[1] or cv-abs(edn) < args.range[0]:
                print("unable to print parameter "+allpars[np]+", dimension too large: "+str(cv)+" +"+str(eup)+" "+str(edn))
            else:
                outfile.write("  \\draw["+style+"] (" +str(cv-abs(edn)) + "," + str(np-npar+offset) + ") -- (" + str(cv+abs(eup)) + "," + str(np-npar+offset) + ");")
                outfile.write("  \\node[dot,"+style+"] at ("+str(cv)+","+str(np-npar+offset)+ ") {};\n")
    from numpy import arange
    for x in arange(int(args.range[0]),int(args.range[1])+1.1,step=0.25):
        outfile.write("\\draw[black] (" +str(x)+ "," + str(-npar-0.5+0.2) + ") -- (" +str(x)+ "," +str(-npar-0.5)+ ") node [axlbl,below=3pt]{" +str(x)+ "};\n")
#        outfile.write("\\node[axlbl,black,below left=5pt] at (" +str(args.range[0]-0.1)+ "," +str(-npar-0.5)+ ") {$\\frac{\\theta-\\theta_0}{\\Delta\\theta_0}$};\n")
        outfile.write("\\draw[black] (" +str(int(args.range[0])-0.1)+ "," +str(-npar-0.5)+ ") -- (" +str(int(args.range[1])+0.1)+ "," +str(-npar-0.5)+ ");\n")
        for x in range(int(args.range[0]),int(args.range[1])):
            outfile.write("\\draw[dashed,black] (" +str(x)+ "," +str(-npar-0.5)+ ") -- (" +str(x)+ "," +str(0)+ ");\n")
                
    outfile.write("\\begin{pgfonlayer}{background}\n")
    outfile.write("  \\foreach \\i in ")
    if npar>3: outfile.write("{-1,-3,...,"+str(-2*((npar+1)/2))+"}")
    else: outfile.write("{1}")
    outfile.write("{\\fill[fill=white!90!black] let \\p1 = (current bounding box.east), \\p2 = (current bounding box.west) in (\\x1,\\i-0.5) rectangle (\\x2,\\i+0.5); }\n")
    outfile.write("\\end{pgfonlayer}\n")
    outfile.write("\\end{tikzpicture}\n")
    writefoot(outfile)
