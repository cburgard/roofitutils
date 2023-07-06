
from RooFitUtils.util import getThreshold
histocolors = ["blue!50","red","green","gray","orange","violet!20","cyan","violet","brown"]
thresholdColors = ["blue","green","yellow","orange","red"]
thresholdStyles = ["solid,fill,opacity=.5","dashed","loosely dashed","dotted","loosely dotted"]


def writehead(stream,atlas=True,varwidth=None):
    stream.write("% this plot was automatically generated by RooFitUtils/pgfplotter.py\n")
    stream.write("% some files written this way can occasionally become very large\n")
    stream.write("% should you encounter an error like this one while trying to compile this file\n")
    stream.write("%     ! TeX capacity exceeded, sorry [main memory size=12000000].\n")
    stream.write("% try the following:\n")
    stream.write("%  - use lualatex or xelatex\n")
    stream.write("%  - increase your tex main memory size by editing your 'texmf.cnf'\n")
    stream.write("%    which you can find by typing 'kpsewhich -a texmf.cnf' in your terminal\n")
    stream.write("%    and change the entry 'main_memory' to some very large number, e.g. 'main_memory.xelatex = 500000000'\n")
    stream.write("%    finally, run 'fmtutil-sys --all' to update your configuration\n")                    
    stream.write("\\documentclass[margin=1pt,")
    if varwidth: stream.write("varwidth="+varwidth)
    stream.write("]{standalone}\n")
    stream.write("\\usepackage{scalerel}\n")
    stream.write("\\usepackage{xcolor}\n")    
    stream.write("\\usepackage{pgfplots,tikz}\n")
    stream.write("\\usepackage{iftex}\n")
    stream.write("\\usetikzlibrary{calc}\n")
    if atlas:
        stream.write("\\ifPDFTeX\n\\usepackage[scaled=1]{helvet}\\usepackage{sansmath}\n")
        stream.write("\\else\n\\usepackage{fontspec}\\setsansfont{TeX Gyre Heros}\n\\fi\n")
        stream.write("\\usepackage[helvet]{sfmath}\n")
    stream.write("\\usepackage{amsmath,latexsym}\n")
    stream.write("\\usetikzlibrary{shapes,positioning,patterns}\n")
    stream.write("\\tikzset{cross/.style={cross out, draw=black, minimum size=2*(#1-\pgflinewidth), inner sep=0pt, outer sep=0pt},cross/.default={3pt}}\n")
    stream.write("\\pgfplotsset{every axis legend/.append style={draw=none},")
    stream.write("    /pgfplots/area legend/.style={\n")
    stream.write("      /pgfplots/legend image code/.code={\n")
    stream.write("        \\fill[draw] (0pt,0pt) ellipse [x radius=\\pgfplotbarwidth,y radius=0.3em];\n")
    stream.write("        \\draw (1.5pt,1.5pt) -- (-1.5pt,-1.5pt);\n")
    stream.write("        \\draw (1.5pt,-1.5pt) -- (-1.5pt,1.5pt);\n")      
    stream.write("}, },\n")
    stream.write("}\n")
    
    stream.write("\\pgfplotsset{compat=newest}\n")
    stream.write("\\begin{document}\n")
    stream.write("\\pgfplotsset{scaled x ticks=false}\n")
    stream.write("\\definecolor{myyellow}{rgb}{0.96,0.742,0.29}\n")
    stream.write("\\definecolor{myblue}{rgb}{0.1,0.32,0.738}\n")
    stream.write("\\definecolor{myotheryellow}{rgb}{0.98828125,0.5625,0.13671875}\n")
    stream.write("\\definecolor{myotherblue}{rgb}{0.18359375,0.3515625,0.7578125}\n")
    if atlas:
        if atlas == True:
            atlaslabel = "Internal"
        else:
            atlaslabel = atlas
        stream.write("\\providecommand\\ATLASlabel{"+str(atlaslabel)+"}\n")
        stream.write("\\ifpdftex\n")
        stream.write("\\renewcommand\\sfdefault{phv}\n")
        stream.write("\\renewcommand\\rmdefault{phv}\n")
        stream.write("\\renewcommand\\ttdefault{pcr}\n")
        stream.write("\\fi\n")
        stream.write("\\font\\greekcapstenrm=cmr10\n")
        stream.write("\\font\\greekcapssevenrm=cmr7\n")
        stream.write("\\font\\greekcapsfiverm=cmr5\n")
        stream.write("\\newfam\\greekcapsfam\n")
        stream.write("\\textfont\\greekcapsfam=\\greekcapstenrm\n")
        stream.write("\\scriptfont\\greekcapsfam=\\greekcapssevenrm\n")
        stream.write("\\scriptscriptfont\\greekcapsfam=\\greekcapsfiverm\n")
        stream.write("\\let\\tmpLambda=\\Lambda \\def\\Lambda{{\\fam\\greekcapsfam\\tmpLambda}}\n")

        

def writefoot(stream):
    stream.write("\\end{document}\n")

def writeATLAS(outfile,label="\\ATLASlabel",inside=True,
               labels=["\\scriptsize{$\\sqrt{s}=$13 TeV, 139 fb$^{\\scriptsize{-1}}$}",
                       "\\scriptsize{$m_{\\scalebox{.9}{$\\scriptstyle H$}}=$ 125.09 GeV, $|y_{\\scalebox{.9}{$\\scriptstyle H$}}|$ $<$ 2.5}"],labelspread=1.8):
    from RooFitUtils.io import texify
    if inside:
        outfile.write("\\node (atlas) at (rel axis cs:0.025,0.975) [anchor=north west]{\\textsf{\\textbf{\\textit{ATLAS}}}};\n")
        outfile.write("\\node (atlaslabel) at (atlas.east) [anchor=west]{\\textsf{"+label+"}};\n")        
    else:
        outfile.write("\\node (atlas) at (rel axis cs:0,1) [scale=2,above right]{\\textsf{\\textbf{\\textit{ATLAS}}}};\n")
        outfile.write("\\node (atlaslabel) at (atlas.east) [scale=2,anchor=west]{\\textsf{"+label+"}};\n")
    if inside:
        for i in range(0,len(labels)):
            outfile.write("\\node at (atlas.west) [anchor=west,scale=0.8,yshift=-"+str(labelspread*(i+1))+"em]{"+texify(labels[i])+"};\n")
    else:
        outfile.write("\\node at (atlaslabel.east) [anchor=west,scale=1.5]{"+", ".join(map(texify,labels))+"};\n")

def concat(strlist):
    string = ""
    for x in range(0,len(strlist)):
        string = string +","+ strlist[x]
    return string[1:len(string)] # remove comma which is the first character

def writepoiset(idx,poinames,allpois,outfile,style,poiopts,spread,printvalues):
    from math import isnan
    from RooFitUtils.util import formatNumberPDG,formatPDG,isdict
    from RooFitUtils.io import texify,readparameter
    color = style.get("color","black")
    count = 0
    shift = style.get("labelshift",str(1+4*idx)+"em")
    if printvalues and "label" in style.keys():
        outfile.write("  \\node[xshift="+shift+",anchor=base] at ({rel axis cs:1,0}|-{axis cs:0,"+str(float(spread*len(poinames)))+"}) {\\textsf{"+texify(style.get("label",""))+"}};\n")
    for poi in poinames:
        opts = poiopts.get(poi,{})
        if isdict(opts):
            scale= opts.get("scale",1.)        
        if not poi in allpois.keys(): continue
        cvs,intervals = readparameter(allpois[poi])
        
        for lo,hi in intervals:
            outfile.write("  \\draw[color="+color+","+style.get("style","solid")+"] (axis cs:{:.5f},{:.2f})--(axis cs:{:.5f},{:.2f});\n".format(scale*lo,spread*count,scale*hi,spread*count))
        for v in cvs:
            outfile.write("  \\node[circle,fill,inner sep=2pt,color="+color+","+style.get("style","draw=none")+"] at (axis cs:{:.5f},{:.2f})".format(float(scale*v),float(spread*count))+ "{};\n")
        if printvalues:
            if len(cvs) == 1 and len(intervals) == 1:
                v = cvs[0]
                lo = intervals[0][0]-v
                hi = intervals[0][1]-v
                outfile.write("  \\node[xshift="+shift+",anchor=base west] at ({rel axis cs:1,0}|-{axis cs:0,"+str(float(spread*count))+"}) {")                    
                if style.get("point",True) and style.get("error",True):
                    outfile.write(formatPDG(v,hi,lo,lap=True,show="v+-"))
                elif style.get("error",True):
                    outfile.write(formatPDG(v,hi,lo,lap=True,show="v"))
                elif style.get("point",True):
                    outfile.write(formatNumberPDG(v))
                outfile.write("};\n")
        count = count+1


        
def writepois(atlas,pois,allsets_input,outfilename,plotlabels=[],range=[-2,2],smval=0,spread=1,printvalues=False):
    """write a POI plot to a pgfplots tex file"""
    from RooFitUtils.io import texify,readparameter
    from RooFitUtils.interpolate import inf
    from RooFitUtils.util import parsepois,parsegroups,islist,isdict,parsedict,flattened
    
    if islist(allsets_input):
        allsets = allsets_input
    elif isdict(allsets_input):
        allsets = [(parsedict(key),value) for key,value in allsets_input.items()]
    
    poinames,poiopts = parsepois(pois,options=True)
    poinames.reverse()
    groups = parsegroups(pois)

    if range[0] != None:
        minval = range[0]
    else:
        minval = inf
        for key,poiset in allsets:
            for poi in poinames:
                if not poi in poiset.keys(): continue
                cvs,intervals = readparameter(poiset[poi])
                minval = min(minval,min(cvs + list(flattened(intervals))))
    if range[1] != None:
        maxval = range[1]
    else:
        maxval = -inf
        for key,poiset in allsets:
            for poi in poinames:
                if not poi in poiset.keys(): continue
                cvs,intervals = readparameter(poiset[poi])
                maxval = max(maxval,max(cvs + list(flattened(intervals))))
    
    with open(outfilename,"w") as outfile:
        writehead(outfile)
        outfile.write("\\begin{tikzpicture}\n")
        outfile.write("\\begin{axis}[\n")
        outfile.write("    width = 0.8\\textwidth,\n")        
        outfile.write("    width = 0.8\\textwidth,\n")
        outfile.write("    y=1cm,\n")
        outfile.write("    xlabel = {\\textsf{Parameter Value}}, \n")
        outfile.write("    clip = false,\n")
        outfile.write("    ymin=-0.5,\n")
        ymax = spread*(len(poinames)-0.5)
        if atlas:
            ymax += 1+len(plotlabels)
        outfile.write("    ymax= "+str(ymax)+ ",\n")
        outfile.write("    xmin="+str(minval)+",\n")
        outfile.write("    xmax="+str(maxval)+",\n")
        outfile.write("    minor tick num=4,\n")
        outfile.write("    ytick style={draw=none},\n")
        outfile.write("    yticklabels=\empty,\n")
        outfile.write("    enlarge x limits=true,\n")
        outfile.write("    xticklabel style={/pgf/number format/fixed,font=\sffamily},\n")
        outfile.write("    scaled ticks=false,\n")
        outfile.write("]\n")
        if atlas: writeATLAS(outfile,atlas,inside=True,labels=plotlabels)            
        count = 0
        if smval > minval and smval < maxval:
            outfile.write("\\draw ("+str(smval)+",-1) -- ("+str(smval)+","+str(spread*(len(poinames)-0.5))+");\n")
        for poi in poinames:
            outfile.write("\\node at ({rel axis cs:0,0}|-{axis cs:0,"+ str(spread*count)+"}) [anchor = east]{\\textsf{"+texify(poi))
            opts = poiopts.get(poi,{})
            if isdict(opts):
                scale= opts.get("scale",1.)
            else:
                scale = 1.
            if scale != 1.:
                outfile.write(" ($\\times {:g}$)".format(scale))
            outfile.write("}};\n")
            count += 1
        count = 0
        
        for key,poiset in allsets:
            writepoiset(count,poinames,poiset,outfile,key,poiopts,spread,printvalues)
            count += 1
                
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)
        print("wrote "+outfilename)

def makelabels(listofstrings):
    x = [ i.replace("_","\_") for i in listofstrings]
    return x

def guessanchor(angle):
    if angle < 45:
        return "north"
    else:
        return "north east"
    
def writematrix(atlas,xcoords_orig,ycoords_orig,allvalues,outfilename,minval=None,maxval=None,rotatelabels=90,plotlabels=[],showall=False,flip=False,axlabel=None,centralval=None,xlabels=None,ylabels=None):
    """write a correlation matrix to a pgfplots tex file"""
    from RooFitUtils.util import flipped,parsegroups,parsepois,flattened
    from RooFitUtils.io import texify
    xgroups = None
    ygroups = None
    xgroups = parsegroups(xcoords_orig)
    xcoords_orig = parsepois(xcoords_orig)
    ygroups = parsegroups(ycoords_orig)
    ycoords_orig = parsepois(ycoords_orig)
    if not xlabels:
        xlabels = [ texify(i) for i in xcoords_orig ]
    if not ylabels:
        ylabels = [ texify(i) for i in flipped(ycoords_orig,flip) ]
    xcoords = [ i.replace("_","") for i in xcoords_orig]
    ycoords = [ i.replace("_","") for i in flipped(ycoords_orig,flip) ]
    if centralval:
        ex = max(map(lambda x:abs(x-centralval),flattened(allvalues)))
        minval = centralval-ex
        maxval = centralval+ex        
    else:
        if not minval: minval = min(flattened(allvalues))
        if not maxval: maxval = max(flattened(allvalues))
    if len(ycoords) != len(allvalues):
        raise RuntimeError("incompatible lengths in y")
    if len(xcoords) != len(allvalues[0]):
        raise RuntimeError("incompatible lengths in x")
    with open(outfilename,"w") as outfile:
        writehead(outfile)
        outfile.write("\\begin{tikzpicture}[\n")
        if atlas:
            outfile.write("  font={\sffamily},\n")
        from RooFitUtils.util import sgnstr,formatNumber
        outfile.write("  map color/.code={\\pgfmathparse{"+str(int(-minval*(1000./(maxval-minval))))+" + "+str(int(1000./(maxval-minval)))+"*#1}\\pgfplotscolormapdefinemappedcolor{\\pgfmathresult}},\n")
        outfile.write("  meta/.style={map color=#1,minimum size=3em,fill=mapped color}\n")        
        outfile.write("]\n")
        outfile.write("\\begin{axis}[\n")
        outfile.write("    colormap={bluewhiteyellow}{color=(myyellow) color=(white) color=(myblue)},\n")
        outfile.write("    colormap={bluewhiteyellow-texthighlighting}{rgb(0pt)=(0.,0.,0.);rgb(900pt)=(0.,0.,0.);rgb(901pt)=(1.,1.,1.);rgb(1000pt)=(1.,1.,1.);},\n")
        outfile.write("    nodes near coords style={anchor=center,font=\\footnotesize,/pgf/number format/fixed,/pgf/number format/precision=2,color of colormap=\pgfplotspointmetatransformed of bluewhiteyellow-texthighlighting")
        if showall: outfile.write(",/pgf/number format/fixed zerofill")
        outfile.write("},\n")
        outfile.write("    clip = false,\n")
        outfile.write("    colorbar,\n")
        outfile.write("    colormap name={bluewhiteyellow},\n")
        outfile.write("    x=3em,\n")
        outfile.write("    y=3em,\n")
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
        outfile.write("    xticklabels={,"+ concat(xlabels) + "},\n") 
        outfile.write("    yticklabels={,"+ concat(ylabels) + "},\n") 
        outfile.write("    axis on top,\n")
        outfile.write("    layers/decorated plot/.define layer set={background,main,foreground}{},\n")
        outfile.write("    clip mode=individual,\n")
        outfile.write("    set layers=decorated plot,\n")
        outfile.write("    x tick label style={font={\sffamily},scale=1.5,anchor="+guessanchor(rotatelabels)+",rotate="+str(rotatelabels)+"},\n")
        outfile.write("    y tick label style={font={\sffamily},scale=1.5},\n")
        outfile.write("    colorbar style={y tick label style={scale=1.5,/pgf/number format/fixed}")
        if axlabel: outfile.write(",ylabel={"+axlabel+"},y label style={anchor=east,yshift=-2cm,scale=2,at={(axis description cs:1,1)}}")
        outfile.write("},\n")        
        outfile.write("    tick style={draw=none}\n")
        outfile.write("]\n")
        outfile.write("\\addplot[only marks,mark=square*,scatter,mark size=1.5em,scatter src=explicit,nodes near coords*] coordinates {\n")
        for x in range(0,len(xcoords)):
            for y in range(0,len(ycoords)):
                if showall or abs(allvalues[y][x]) > 0.005:
                    outfile.write("  ("+str(xcoords[x])+","+str(ycoords[y])+") ["+str(allvalues[len(ycoords) - y - 1 if flip else y][x])+"]\n")
        outfile.write("};\n")
        if ygroups:
            iy = 0
            if flip:
                ycoord = len(ycoords)-0.5
            else:
                ycoord = 0.5
            for g in ygroups:
                if flip:
                    step = -g[-1]
                else:
                    step = g[-1]                
                if iy > 0:
                    outfile.write("\\draw[ultra thick,dashed,on layer=foreground] (axis cs:{[normalized]"+str(len(xcoords)-0.5)+"},{[normalized]"+str(ycoord)+"}) -- (axis cs:{[normalized]-0.5},{[normalized]"+str(ycoord)+"})  -- ++(-30em,0em);\n")
                outfile.write("\\node[on layer=foreground,xshift=-30em,anchor=north,scale=2,rotate=90] at (axis cs:{[normalized]-0.5},{[normalized]"+str(ycoord+0.5*step)+"}) {"+texify(g[0])+"};\n")                    
                iy += 1
                ycoord += step
        if xgroups:
            ix = 0
            xcoord = -0.5
            for g in xgroups:
                step = g[-1]                
                if ix > 0:
                    outfile.write("\\draw[ultra thick,dashed,on layer=foreground] (axis cs:{[normalized]"+str(xcoord)+"},{[normalized]"+str(len(ycoords)-0.5)+"}) -- (axis cs:{[normalized]"+str(xcoord)+"},{[normalized]-0.5});\n")
                ix += 1
                xcoord += step                
        if atlas: writeATLAS(outfile,atlas,inside=False,labels=plotlabels)                        
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)
        print("wrote "+outfilename)

def writecorrmatrix(atlas,parslist,allcorrs,outfilename):
    writematrix(atlas,parslist,parslist,allcorrs,outfilename,minval=-1,maxval=1,rotatelabels=45,axlabel="$\\rho(X,Y)$",flip=True,showall=True)

def writescans1d(atlas,par,allscans,outfilename,percent_thresholds=None,drawpoints=False,ymax=None,rangex=None,plotlabels=[],otherscans1d=[],axis_options=[],append=None):
    from RooFitUtils.util import make1dscan

    for otherscan in otherscans1d:
        from RooFitUtils.util import mergescans1d
        if len(otherscan) > 0: allscans = mergescans1d(allscans,otherscan)

    """write a bunch of 1d scans to a pgfplots tex file"""
    with open(outfilename,"w") as outfile:
        writehead(outfile)
        allvals = [ v[0] for curves in allscans.values() for scan in curves.values() for v in scan.keys()]
        domain = "domain={0:f}:{1:f}".format(min(allvals),max(allvals))
        outfile.write("\\begin{tikzpicture}\n")
        outfile.write("\\begin{axis}[\n")
        outfile.write("    ymin=0,\n")
        for opt in axis_options:
            outfile.write("    "+opt+",\n")            
        if ymax:
            outfile.write("    ymax="+str(ymax)+",\n")
        outfile.write("    "+domain+",\n")
        outfile.write("    legend pos=north east,legend style={draw=none},\n")
        outfile.write("    xlabel=${0:s}$, ylabel=$-2\\ \\ln \\Lambda$,\n".format(par.strip("$")))
        outfile.write("    every axis x label/.style={at={(axis description cs:1.0,-0.1)},anchor=north east},\n")
        outfile.write("    every axis y label/.style={at={(axis description cs:-0.1,1.0)},rotate=90,anchor=south east},\n")
        if rangex:
            outfile.write("    xmin={0:f},xmax={1:f}\n".format(*rangex))
        else:
            outfile.write("    xmin={0:f},xmax={1:f}\n".format(min(allvals),max(allvals)))            
        outfile.write("]\n")
        if atlas: writeATLAS(outfile,atlas,inside=True,labels=plotlabels)
        for pnamelist,curve in allscans.items():
            for options,scan in curve.items():
                print("writing scan for "+pnamelist[0])
                writescan1d(pnamelist[0],par,make1dscan(scan),options,outfile,percent_thresholds,drawpoints,ymax)
        outfile.write("\\addplot[gray,densely dashed,thick] {1};\n")
        outfile.write("\\addplot[gray,densely dashed,thick] {4};\n")
        if append:
            for line in append:            
                outfile.write(append)
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)
        print("wrote "+outfilename)

def writescan1d(parname,parlabel,allpoints,options,outfile,percent_thresholds,drawpoints=False,ymax=None):
    """write a single 1d sncan to a pgfplots tex file"""
    from math import isnan
    from RooFitUtils.interpolate import findminimum,findcrossings
    from RooFitUtils.util import graphrescale,formatPDG,parsedict
    # obtaining the - 2 log Lambda(x), where Lambda = LL(x)/LL(x0)
    nllmin = findminimum(allpoints)
    points = graphrescale(allpoints,nllmin,2)
    style = parsedict(options)
    outfile.write("\\addplot["+options+",very thick,mark=none,smooth] coordinates {\n")
    for x,y in points:  outfile.write("    ({0:f},{1:f})\n".format(x,y))
    outfile.write("};\n")
    if drawpoints:
        outfile.write("\\addplot[draw=none,mark=x] coordinates {\n")
        for x,y in points:  outfile.write("    ({0:f},{1:f})\n".format(x,y))
        outfile.write("};\n")
 
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
                    if "title" in style.keys():
                        title = style["title"]
                    else:
                        title = parlabel.strip("$")
                    outfile.write("\\addlegendentry{{${:s} = {:s}$}}".format(title,formatPDG(cv,up,down)))
                s = s + ", {:.3f}% CL = +{:f} -{:f}".format(100*percent_thresholds[i],abs(up),abs(down))
       # print(s)

def ishexcolor(s):
    if len(s) == 6 and all(c in "1234567890ABCDEFabcdef" for c in s):
        return True
    return False
       
def writescans2d(atlas,labels,scans2d,outfilename,extrapoints,npoints,percent_thresholds,plotlabels=[],otherscans2d=[],flipAxes=False,contourAlg="skimage",smooth=False,axis_options=[],append=None):
    """write a bunch of 2d scans to a pgfplots tex file"""
    from RooFitUtils.util import parsedict
    from RooFitUtils.io import texify    
    with open(outfilename,"w") as outfile:
        writehead(outfile)
        for pnamelist,scan in scans2d.items():
            for drawopts in scan.keys():
                opts = parsedict(drawopts)
                if "color" in opts.keys():
                    color = opts["color"]
                    if ishexcolor(color):
                        outfile.write("\\definecolor{"+color+"}{HTML}{"+color.upper()+"}\n")
        outfile.write("\\begin{tikzpicture}[\n")
        if atlas:
            outfile.write("  font={\sffamily}\n")
        outfile.write("]\n")
        outfile.write("\\begin{axis}[area legend,\n")
        for opt in axis_options:
            outfile.write("  "+opt+",\n")
        if atlas:
            outfile.write("legend pos=north east,legend style={draw=none,fill=none},legend cell align={left},\n")
            outfile.write("xticklabel={\\pgfmathprintnumber[fixed,assume math mode=true]{\\tick}},\n")
            outfile.write("yticklabel={\\pgfmathprintnumber[fixed,assume math mode=true]{\\tick}},\n")
        if len(labels) == 2:
            if flipAxes:
                outfile.write("    ylabel={"+texify(labels[0],math=True)+"},\n")
                outfile.write("    xlabel={"+texify(labels[1],math=True)+"},\n")
            else:                                           
                outfile.write("    xlabel={"+texify(labels[0],math=True)+"},\n")
                outfile.write("    ylabel={"+texify(labels[1],math=True)+"},\n")
        outfile.write("]\n")
        if atlas: writeATLAS(outfile,atlas,inside=True,labels=plotlabels)
        for pnamelist,scan in scans2d.items():
            for drawopts,points in scan.items():
                morepoints = [ s[pnamelist][drawopts] for s in otherscans2d if pnamelist in s.keys() ]
                writescan2d(points,outfile,percent_thresholds,parsedict(drawopts),npoints,morepoints=morepoints,flipAxes=flipAxes,contourAlg=contourAlg,smooth=smooth)
        for drawopts,points in extrapoints.items():
            writepoints2d(points,outfile,parsedict(drawopts),flipAxes=flipAxes)
        if append:
            for line in append:
                outfile.write(line)            
        outfile.write("\\end{axis}\n")
        outfile.write("\\end{tikzpicture}\n")
        writefoot(outfile)
        print("wrote "+outfilename)

def writepoints2d(points,outfile,style,flipAxes=False,keys=None):
    outfile.write("\\addplot[mark=x,mark options={scale=.5},only marks,draw="+style.get("color","black")+"] coordinates {\n")
    i = 0
    red = int(style.get("reduce",0))
    from random import randint
    for point in points:
        i = i + 1
        if red > 0 and randint(0, red) != 0: continue
        if not keys:
            keys = list(point.keys())
        if flipAxes:
            outfile.write("    ({:f},{:f})\n".format(point[keys[1]],point[keys[0]]))
        else:
            outfile.write("    ({:f},{:f})\n".format(point[keys[0]],point[keys[1]]))
    outfile.write("};\n")


def writescan2d(allpoints,outfile,percent_thresholds,style,npoints,morepoints=[],flipAxes=False,contourAlg="skimage",smooth=False):
    """write a single 2d scan to a pgfplots tex file"""
    from RooFitUtils.interpolate import findcontours
    thresholds = [ getThreshold(p,2)*0.5 for p in percent_thresholds ]
    contours,minimum = findcontours(allpoints,thresholds,smooth,npoints,algorithm=contourAlg,morepoints=morepoints)
    if not contours:
        print("no contours found!")
    outfile.write("\\draw (axis cs:")
    if flipAxes:
        outfile.write("{:f},{:f}".format(minimum[1],minimum[0]))
    else:
        outfile.write("{:f},{:f}".format(minimum[0],minimum[1]))
    outfile.write(") node[cross,color="+style.get("color","black")+"] {};\n")

    first = True
    i = 0
    ignore = [ c.split("/") for c in style.get("ignore_contours","").split(";") if c ]
    for v,conts in zip(thresholds,contours):
        icont = 0
        skip = False
        for c in conts:
            skip = False
            for elem in ignore:
                if (elem[0] == '*' or int(elem[0]) == i) and (elem[1] == '*' or int(elem[1]) == icont):
                    skip = True
            if skip:
                outfile.write("% contour {:d} of {:f} skipped\n".format(icont,v))                
                continue
            icont = icont+1
            outfile.write("% contour {:d} of {:f}\n".format(icont-1,v))
            outfile.write("\\addplot[color="+style.get("color","black")+","+thresholdStyles[i]+",mark=none,smooth")
            if not first: outfile.write(",forget plot")
            outfile.write("] coordinates {\n")
            for x,y in c:
                if flipAxes:
                    outfile.write("    ({:f},{:f})\n".format(y,x))
                else:
                    outfile.write("    ({:f},{:f})\n".format(x,y))
            outfile.write("} ;\n")
            if first and "title" in style.keys():
                outfile.write("\\addlegendentry{"+style["title"]+"};\n")
            first=False
        i = i+1

def getColorDefStringLaTeX(name,color):
    c = ROOT.gROOT.GetColor(color)
    return getColorDefStringLaTeX(name,c)

def getColorDefStringLaTeX(name,color):
  if not color: return ""
  r,g,b = 0.,0.,0.
  color.GetRGB(r,g,b);
  return "\\definecolor{"+name+"}{rgb}{"+str(r)+","+str(g)+","+str(b)+"}"

def writepull(args,results,outfile,parname,ypos,offset=True,style="red"):
    from math import isnan
    res = results[parname]
    ires = 0
    cv = res['val']
    if "eUp" in res.keys():
        eup = res["eUp"]
    else:
        eup = res["err"]
    if "eDn" in res.keys():
        edn = res["eDn"]
    else:
        edn = res["err"]
    if isnan(edn): edn = 0
    if isnan(eup): eup = 0
    ires = ires+1
    if offset:
        voffset = 0
    else:
        voffset = 0
    if cv+abs(eup) > args.range[1] or cv-abs(edn) < args.range[0]:
        print("unable to print parameter "+parname+", dimension too large: "+str(cv)+" +"+str(eup)+" "+str(edn))
    else:
       outfile.write("  \\draw[pull,"+style+"] ({:.5f},{:.2f})--({:.5f},{:.2f});".format(cv-abs(edn),ypos+voffset,cv+abs(eup),ypos+voffset))
       outfile.write("  \\node[dot,"+style+"] at ({:.5f},{:.2f})".format(cv,ypos+voffset)+ "{};\n")

def writeranking(args,ranking,outfile,ypos,impscale=1,ysize="5pt"):
    from math import isnan
    up,dn = ranking
    if isnan(up): up=0
    if isnan(dn): dn=0
    outfile.write("  \\fill[myblue]   ($({:.5f},{:.2f})+(0,-{:s})$) rectangle ($({:.5f},{:.2f})+(0,{:s})$);".format(float(impscale*dn),ypos,ysize,0, ypos,ysize))
    outfile.write("  \\fill[myyellow] ($({:.5f},{:.2f})+(0,-{:s})$) rectangle ($({:.5f},{:.2f})+(0,{:s})$);".format(0, ypos,ysize,float(impscale*up),ypos,ysize))
    outfile.write("\n")
        
def writepullnumbers(args,results,outfile,parname,ypos,hoffset,numbers=False):
    from math import isnan
    res = results[parname]
    ires = 0
    cv = res['val']
    if "eUp" in res.keys():
        eup = res["eUp"]
    else:
        eup = res["err"]
    if "eDn" in res.keys():
        edn = res["eDn"]
    else:
        edn = res["err"]
    if not isnan(eup): seup = "+{:.4f}".format(abs(eup))
    else:          seup = "?"
    if not isnan(edn): sedn = "-{:.4f}".format(abs(edn))
    else:          sedn = "?"
    snom = "{:.4f}".format(cv)
    ires = ires+1
    hpos = args.range[1]
    hoffset += 5*float(ires)
    outfile.write("  \\node[lbl,xshift="+str(hoffset)+"cm,anchor=east] at ("+str(hpos)+","+str(ypos)+") ")
    if numbers:
        outfile.write(" {$"+snom+"^{"+seup+"}_{"+sedn+"}$};")
    else:
        outfile.write(" {};")

def writerankinghead(args,outfile,allpars,zoom=10,poiname="\mu"):
    from RooFitUtils.util import roundAutoUp,sign
    npar = len(allpars)
    from numpy import arange
    outfile.write("% top axis\n")
    outfile.write("  \\draw[blue] (" +str(int(args.range[0])-0.1)+ "," +str(-0.5-npar)+ ") -- (" +str(int(args.range[1])+0.1)+ "," +str(-0.5-npar)+ ") node [pos=1,anchor=north east,yshift=1.1cm]{$\Delta "+poiname+"$};\n")
    ii = 0
    x0,x1 = round(args.range[0]/zoom,1),round(args.range[1]/zoom,1)
    step = round(0.2*abs(x1-x0),1)
    for x in arange(x0,x1+step,step=step):
        outfile.write("  \\draw[blue] (" +str(zoom*x)+ "," + str(-0.5-npar) + ") -- ++ (0,0.2) node [axlbl,above=3pt]{" + "{:.1f}".format(x)+ "};\n")
        ii = ii + 1

def writeparametershead(args,outfile):
    xunit = 3
    yunit = .5
    outfile.write("\\begin{tikzpicture}[x="+str(xunit)+"cm,y=-"+str(yunit)+"cm,%\n")
    outfile.write("  lbl/.style={scale=1,anchor=west},%\n")
    outfile.write("  axlbl/.style={scale=0.7,anchor=center},%\n")
    outfile.write("  pull/.style={{|[scale=1]}-{|[scale=1]}},%\n")
    outfile.write("  dot/.style={circle,fill,inner sep=1pt},\n")
    outfile.write("  every node/.append style={font=\\sffamily}\n")
    outfile.write("]\n")
    outfile.write("\\pgfdeclarelayer{background}\\pgfsetlayers{background,main}\n")
    
def writeparameters(args,outfile,allpars,labels):    
    npar = len(allpars)
    for np in range(0,npar):
        text = allpars[np]
        if labels == "r":
            outfile.write("  \\node[lbl,xshift=1cm,anchor=west] at ("+str(args.range[1])+","+str(np-npar)+") ")
        else:
            outfile.write("  \\node[lbl,xshift=-1cm,anchor=east] at ("+str(args.range[0])+","+str(np-npar)+") ")            
        outfile.write("{")            
        if "{" in text and not "$" in text:
            outfile.write("\\ensuremath{"+text+"}")
        else:
            outfile.write(text.replace("_","\\_"))
        outfile.write("};\n")

def writepullsfoot(args,outfile,allpars):
    npar = len(allpars)
    from math import floor,ceil
    from numpy import arange
    outfile.write("% bottom axis\n")
    for x in arange(floor(args.range[0]),ceil(args.range[1])+.1,step=0.25):
        outfile.write("  \\draw[black] (" +str(x)+ "," + str(-0.5+0.2) + ") -- (" +str(x)+ "," +str(-0.5)+ ") node [axlbl,below=3pt]{" +str(x)+ "};\n")
    outfile.write("  \\draw[black] (" +str(int(args.range[0])-0.1)+ "," +str(-0.5)+ ") -- (" +str(int(args.range[1])+0.1)+ "," +str(-0.5)+ ") node [pos=1,anchor=north east,yshift=-.5cm]{$(\hat{\\theta}-\\theta_{0})/\Delta\\theta$};\n")
    for x in range(int(floor(args.range[0])),int(ceil(args.range[1]+1))):
        outfile.write("  \\draw[dashed,black] (" +str(x)+ "," +str(-0.5)+ ") -- (" +str(x)+ "," +str(-0.5-npar)+ ");\n")
    outfile.write("% highlighting\n")        
    outfile.write("\\begin{pgfonlayer}{background}\n")
    outfile.write("  \\foreach \\i in ")
    if npar>3: outfile.write("{-1,-3,...,"+str(-2*((npar+1)/2))+"}")
    else: outfile.write("{1}")
    outfile.write("{\\fill[fill=white!90!black] let \\p1 = (current bounding box.east), \\p2 = (current bounding box.west) in (\\x1,\\i-0.5) rectangle (\\x2,\\i+0.5); }\n")
    outfile.write("\\end{pgfonlayer}\n")
    outfile.write("\\end{tikzpicture}\n")
        
def writepulls(args,results,outfile,allpars=None,offset=True,labels="r",numbers=False):
    if not allpars:
        allpars = sorted(list(results.keys()))
    writehead(outfile)
    writeparametershead(args,outfile)    
    writeparameters(args,outfile,allpars,labels)
    npar = len(allpars)
    parwidth = max(map(len,allpars))
    for np in range(0,npar):
        writepull(args,results,outfile,allpars[np],np-npar,offset)
        writepullnumers(args,results,outfile,allpars[np],np-npar,0.1*parwidth,numbers)
    writepullsfoot(args,outfile,allpars)        
    writefoot(outfile)   

def plotBarPanel(outfile,data,layout,showlabels,xwidth):
    """plot a single panel of a multipanel bar chart"""
    from RooFitUtils.util import flattened,keys
    from RooFitUtils.io import texify
    allbins = list(keys(layout["bins"]))
    coords = list(map(lambda x:x.replace("_",""),allbins))
    ndata = len(keys(layout["data"]))
    # actually make the plots
    ymin = 0
    ymax = 0    
    plots = []
    i = 0
    for p in keys(layout["data"]):
        style = "[draw=none,fill="+histocolors[i]+"]"
        i += 1
        if type(layout["data"]) == dict:
            if "style" in layout["data"][p].keys():
                style = "[" + layout["data"][p]["style"] + "]"
        plot = ["\\addplot" + style + " coordinates {"]
        scale = 1
        if type(layout["data"]) == dict:
            scale = layout["data"][p].get("scale",1)
        for b in range(0,len(allbins)):
            if allbins[b] in data[p]:
                yval = data[p][allbins[b]] * scale
                ymin = min(yval,ymin)
                ymax = max(yval,ymax)
                plot.append("    ("+coords[b]+","+ str(yval) +")")
        plot.append("};")
        if type(layout["data"]) == dict and "label" in layout["data"][p].keys():
            plot.append("\\addlegendentry{$"+layout["data"][p]["label"]+"$};\n")
        plots.append(plot)
    # generate axes and everything else
    x = layout.get("x",xwidth)
    labels = [ layout["bins"][b].get("label",texify(b)) for b in allbins ]
    outfile.write("\\begin{axis}[\n")
    if ndata == 1:
        outfile.write("    ybar=0pt,\n")
        outfile.write("    bar width="+x+",\n")
    else:
        outfile.write("    ybar=-2pt,\n")
        outfile.write("    bar width=3pt,\n")
    outfile.write("    height="+layout.get("height","3cm")+",\n")    
    outfile.write("    ymin="+str(1.1*ymin)+",ymax="+str(1.1*ymax)+",\n")
    outfile.write("    xmin="+coords[0]+",xmax="+coords[-1]+",\n")    
    outfile.write("    scale only axis,axis on top,\n")        
    outfile.write("    enlarge x limits={abs=.5em},\n")    
    outfile.write("    xminorgrids=true,minor tick num=1,minor grid style={line width=.2pt,draw=gray!50,dashed},\n")
    outfile.write("    legend style={at={(1.0,0.5)},anchor=west,legend columns=1,draw=none,fill=none,nodes={scale=0.5}},legend cell align={left},\n")
    if "y" in layout.keys():
        outfile.write("    ylabel={$"+layout["y"]+"$},\n")
    outfile.write("    symbolic x coords={" + ",".join(coords)+"},\n")
    outfile.write("    xtick={" + ",".join(coords)+"},\n")    
    outfile.write("    scaled ticks=false,\n")
    outfile.write("    y label style={yshift=1.5em,at={(0,0.5)},minimum size=1em,anchor=base,font=\\footnotesize},\n")        
    outfile.write("    y tick label style={/pgf/number format/fixed,scale=0.5,font=\sffamily},\n")    
    outfile.write("    x tick label style={anchor=north east,rotate=60,major tick length=0pt,minor tick length=0pt,scale=0.7,font=\sffamily},\n")
    if showlabels:
        outfile.write("    xticklabels={" + ",".join([ "{$"+b+"$}" if b else "" for b in labels]) + "},\n")
    else:
        outfile.write("    xticklabels=\empty,\n")
    outfile.write("    x="+x+",\n")
    outfile.write("    ]\n")
    for plot in plots:
        outfile.write("\n".join(plot)+"\n")
    for b,style in layout["bins"].items():
        if "box" in style.keys():
            outfile.write("\\draw[xshift=-0.5*"+x+","+style["box"]+"] ({axis cs:"+b+",0}|-{rel axis cs:0,0}) rectangle ++($({"+x+",0pt}|-{rel axis cs:0,1})+(0pt,"+x+")$);\n")
    outfile.write("\\draw (axis cs:"+coords[0]+",0) ++ (-0.5*"+x+",0pt) -- (axis cs:"+coords[-1]+",0) -- ++("+x+",0pt);\n")
    outfile.write("\\end{axis}\n")
    
def plotBars(outfilename,data,layout,xwidth="1em"):
    """plot a multipanel bar chart"""
    with open(outfilename,"wt") as outfile:
        writehead(outfile,varwidth="20cm")
        i = 0
        for panel in layout.values():
            islast = i==len(layout)-1
            outfile.write("\\pgfplotsset{/pgfplots/ybar legend/.style={/pgfplots/legend image code/.code={\\fill[##1] (0cm,.25em) rectangle (\pgfplotbarwidth,-0.25em);}}}\n")
            outfile.write("\\begin{tikzpicture}[]\n")
            plotBarPanel(outfile,data,panel,islast,xwidth)
            outfile.write("\\pgfresetboundingbox\n")
            outfile.write("\\path[use as bounding box] (rel axis cs:-0.05," + str(-2 if islast else 0.)+") rectangle (rel axis cs:1.2," + str(1.05 if i==0 else 1.0)+");\n")
            outfile.write("\\end{tikzpicture}\n\n")
            i += 1
        writefoot(outfile)
        
    
