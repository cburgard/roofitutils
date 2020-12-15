
import os

latex_aliases = {}

def texdict(d):
    global latex_aliases
    latex_aliases.update(d)

def texdef(a,b):
    global latex_aliases
    latex_aliases[a] = b

def texify(string):
    if string in latex_aliases.keys():
        return latex_aliases[string]
    return string.replace("_","\\_")

def writeResult(out,result,writecorrmat):
    if result.min.nll:
        out.write("Minimization: minNll = ")
        out.write(str(result.min.nll))
        out.write("\n")
        for p in result.parameters:
            out.write("{0:s} = {1:g} - {2:g} + {3:g}\n".format(p.name,p.value,abs(p.errLo),abs(p.errHi)))
    if result.fit and writecorrmat:
        matrix = result.fit.correlationHist()
        out.write("Correlations {:d}\n".format(matrix.GetNbinsX()))
        for i in range (0,matrix.GetXaxis().GetNbins()):
            out.write(matrix.GetXaxis().GetBinLabel(i+1)+" ")
        out.write("\n")
        for i in range (0,matrix.GetNbinsX()):
            for j in range(0,matrix.GetNbinsY()):
                out.write(matrix.GetXaxis().GetBinLabel(i+1)+" ")
                out.write(matrix.GetYaxis().GetBinLabel(j+1)+" ")
                out.write("{:.6f} ".format(matrix.GetBinContent(i+1,j+1)))
                out.write("\n")
    for scan in result.scans:
        out.write((" ".join(map(str,scan.parNames))) + " nll status\n")
        for i in range(0,len(scan.nllValues)):
            out.write((" ".join([ str(scan.parValues[i][j]) for j in range(0,len(scan.parNames)) ]))+" "+str(scan.nllValues[i])+" "+str(scan.fitStatus[i])+"\n")

def collectpoints(points,files,label):
    import glob
    filenames = []
    for expression in files:
        filenames.extend(glob.glob(expression))
    if len(filenames) == 0:
        print("no points found in "+",".join(files))
        return
    from RooFitUtils.util import parsedict
    allpoints = []
    for filename in filenames:
        if os.path.isfile(filename):
            with open(filename,'r') as infile:
                for line in infile:
                    point = parsedict(line,float)
                    allpoints.append(point)
    points[label] = allpoints

def reduceparams(allpars,parfilter):
    import re
    skimpars = ""
    parpat = re.compile("^"+parfilter)
    for x in allpars:
        if x.startswith(parfilter):
            skimpars += " "+ texify(x)
    return skimpars[1:len(skimpars)]

def getCovariance_ROOT(infilename,frname,parlist):
    import ROOT
    from os.path import isfile
    if not isfile(infilename):
        raise RuntimeError("unable to open file '{:s}'".format(infilename))    
    infile = ROOT.TFile.Open(infilename,"READ")
    if not infile or not infile.IsOpen():
        raise RuntimeError("unable to open file '{:s}'".format(infilename))
    if frname:
        fr = infile.Get(frname)
    else:
        for obj in infile.GetListOfKeys():
            if obj.GetClassName() == "RooFitResult":
                fr = obj.ReadObj()
                break
    if parlist:
        rooargs = ROOT.RooArgList()
        for p in parlist:
            rooargs.add(fr.floatParsFinal().find(p))
        return fr.reducedCovarianceMatrix(rooargs)
    else:
        return fr.covarianceMatrix()        

def getCovariance(infilename,frname=None,parameters=None):
    mat = getCovariance_ROOT(infilename,frname,parlist=parameters)
    from numpy import array
    return array([ [ mat(i,j) for i in range(0,mat.GetNcols()) ] for j in range(0,mat.GetNrows()) ])

def getCorrelation(infilename,frname=None,parameters=None):
    mat = getCovariance_ROOT(infilename,frname,parlist=parameters)
    from numpy import array
    from math import sqrt
    return array([ [ mat(i,j) / sqrt(mat(i,i) * mat(j,j)) for i in range(0,mat.GetNcols()) ] for j in range(0,mat.GetNrows()) ])

def collectcorrelations(results,filename,parfilter):
    import re
    ncorrpat = re.compile("Correlations[ ](\w+)")
    import glob
    ncorr = 0
    parnames = []
    if os.path.isfile(filename):
        with open(filename,"r") as infile:
             lines = [line for line in infile ]
             for lineno in range(0,len(lines)):
                line = lines[lineno]
                ncorrmatch = ncorrpat.match(line)
                if ncorrmatch:
                    ncorr = int(ncorrmatch.group(1))
                    lineno = lineno + 1
                    parline = lines[lineno]
                    allpars = parline.split(" ")
                    allpars = allpars[0:len(allpars)-1]
                    if len(parfilter) != 0:
                        redpars = reduceparams(allpars,parfilter)
                        parnames.append(redpars)
                    else:
                        for x in allpars: parnames.append(texify(x))
                    pattern = ""
                    for x in range(0, ncorr):
                        pattern += "([-]*\d+.\d+)[ ]"
                    rowpat = re.compile("([a-zA-Z0-9_.]+)[ ]([a-zA-Z0-9_.]+)[ ]"+"([-]*\d+.\d+)")
                    for lineno in range(lineno,len(lines)):
                        rowmatch = rowpat.match(lines[lineno])
                        if rowmatch:
                            par1 = (rowmatch.group(1),rowmatch.group(2))
                            if (len(parfilter) != 0 and par1[0].startswith(parfilter) and par1[1].startswith(parfilter)) or len(parfilter) == 0:
                                results.append( texify((par1[0])) +" "+ texify((par1[1]))+" "+rowmatch.group(3))
    return parnames[0]

def readtables(infilenames):
    from RooFitUtils.util import isstr,striplatex
    if isstr(infilenames):
        infilenames = [infilenames]
    import re
    begintable = re.compile("\\\\begin[ ]*{tabular}[ ]*{(.*)}")
    endtable = re.compile("\\\\end[ ]*{tabular}")
    tables = []
    for infilename in infilenames:
        with open(infilename,"rt") as infile:
            intable = False
            columns = None
            currenttext = ""
            for line in infile:
                if not intable:
                    begin =  begintable.search(line)
                    if begin:
                        intable = True
                        columns = begin.group(1)
                else:
                    end = endtable.search(line)
                    if end:
                        tables.append([list(map(striplatex,row.split("&"))) for row in currenttext.split("\\\\")])
                        intable=False
                    else:
                        currenttext += line
    return tables

def readsummary(infilename,form={"cv":("Central","Tot lo","Tot hi"),"stat":("Central","Stat lo","Stat hi"), "sys":("Central","Syst lo","Syst hi")}):
    import csv
    translations = {}
    results = {k:{} for k in form.keys()}
    with open(infilename,"rt") as infile:
        for line in infile:
            parts = line.split(",")
            if len(translations) == 0:
                for i in range(0,len(parts)):
                    translations[parts[i].strip()] = i
                continue
            parname = parts[translations["POI"]]
            for k in form.keys():
                results[k][parname] = tuple([float(parts[translations[key]]) for key in form[k]])
    return results
            
def collectresults(scans,results,files,label):
    """collect a set of results files and return the contents as a dictionary"""
    import re
    parpat = re.compile("([a-zA-Z0-9_.-]+)[ ]*=[ ]*([0-9.naife+-]+)[ ]*-[ ]*([0-9.naife-]+)[ ]*\+[ ]*([0-9.naife+-]+)[ ]*")
    parpat_legacy = re.compile("[ \t]*([a-zA-Z0-9_.-]+)[ \t=]+([0-9.naife+-]+)[ \t]*\+/-[ \t]*([0-9.naife+-]+).*")
    nllpat = re.compile("Minimization:[ ]*minNll[ ]*=[ ]*([0-9.naife+-]+)")
    if isinstance(files, str): files = [files]
    import glob
    filenames = []
    for expression in files:
        filenames.extend(glob.glob(expression))
    if len(filenames) == 0:
        print("no points found in "+",".join(files))
        return
    for filename in filenames:
        if os.path.isfile(filename):
            if filename.endswith(".root"):
                import ROOT
                infile = ROOT.TFile.Open(filename,"READ")
                for key in infile.GetListOfKeys():
                    obj = key.ReadObj()
                    if obj.InheritsFrom(ROOT.TGraph.Class()):
                        key = (obj.GetXaxis().GetTitle(),)
                        # try to guess if the graph is 2*NLL or just NLL
                        scale = 2
                        if "2" in obj.GetYaxis().GetTitle():
                            scale = 1
                        if key not in scans.keys():
                            scans[key] = {}
                        if label not in scans[key].keys():
                            scans[key][label] = {}                        
                        for i in range(0,obj.GetN()):
                            scans[key][label][(obj.GetX()[i],)] = obj.GetY()[i]*scale
            else:
                with open(filename,'r') as infile:
                    lines = [ line for line in infile ]
                    for lineno in range(0,len(lines)):
                        line = lines[lineno]
                        try:
                            parts = line.split()
                            nllmatch = nllpat.match(line)
                            match = parpat.match(line)
                            match_legacy = parpat_legacy.match(line)
                            if nllmatch:
                                minnll = float(nllmatch.group(1))
                            elif match or "BkgTheory" in line: # todo: cleanup!
                                pname,cv,ed,eu = match.group(1).strip(),match.group(2),match.group(3),match.group(4)
                                result = (float(cv),-float(ed),float(eu))
                                if not pname in results.keys():
                                    results[pname] = {}
                                results[pname][label] = result
                            elif match_legacy:
                                pname,cv,e = match_legacy.group(1).strip(),match_legacy.group(2),match_legacy.group(3)
                                result = (float(cv),-float(e),float(e))
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
                                if int(parts[-1]) == 0: # check if fit ended with status 0
                                    scans[key][label][pvals] = nllval
                        except (KeyError,ValueError) as err:
                            if nllmatch:
                                print("unable to parse line '"+line.strip()+"' in file '"+filename+"', attempted to parse as nll value. "+str(err))
                            elif match:
                                print("unable to parse line '"+line.strip()+"' in file '"+filename+"', attempted to parse as parameter value. "+str(+err))
                            else:
                                print("unable to parse line '"+line.strip()+"' in file '"+filename+"', "+str(err))
                            continue

                for p in scans.keys():
                    if p in results.keys():
                        scans[p][label][results[p][label][0]]=minnll;

def readcsv2dict(filename):
    import csv
    outdict = {}
    with open(filename,"rt") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            outdict[row["name"].strip()] = {k.strip():v.strip() for k,v in row.items()}
    return outdict

def collecthistograms(histograms,cfg,parameters=None):
    from RooFitUtils.util import names
    import re
    infilename = cfg["input"]
    if infilename.endswith(".root"):
        import ROOT
        tfile = ROOT.TFile.Open(infilename,"READ")
        if parameters:
            pars = parameters
        else:
            pars = names(tfile.GetListOfKeys())
        pattern = re.compile(cfg.get("pattern",".*_(?P<p>[A-z0-9]*)_.*"))
        for obj in tfile.GetListOfKeys():
            match = pattern.match(obj.GetName())
            if match:
                histo = obj.ReadObj()
                n = histo.GetXaxis().GetNbins()
                if "p" in match.groupdict().keys():
                    p = match.group("p")
                else:
                    p = cfg["parameter"]
                if not p in histograms.keys():
                    histograms[p] = {}
                if "label" in match.groupdict().keys():
                    label = match.group("label")
                    idx = cfg.get("bin",1)
                    histograms[p][label] = histo.GetBinContent(idx)
                elif "labels" in cfg.keys():
                    for b in cfg["labels"]:
                        idx = histo.GetXaxis().FindBin(b)
                        histograms[p][b] = histo.GetBinContent(idx)
                elif histo.GetXaxis().IsAlphanumeric():
                    for b in range(0,n):
                        histograms[p][histo.GetXaxis().GetBinLabel(b+1)] = histo.GetBinContent(b+1)
    elif infilename.endswith(".xml"):
        import xml.etree.ElementTree as ET
        tree = ET.parse(infilename)
        root = tree.getroot()
        pattern = re.compile(cfg["pattern"])
        subpattern = re.compile(cfg.get("subpattern","(?P<c>[+-]?[ ]*[0-9.]+)[ *]*(?P<p>[A-z][A-z0-9]+)"))
        for node in root:
            if node.tag == "Item":
                text = node.attrib["Name"]
                match = pattern.match(text)
                if not match: continue
                if "label" in match.groupdict().keys():
                    b = match.group("label")
                else:
                    continue
                if "expr" in match.groupdict().keys() and "p" in match.groupdict().keys():
                    fmt = re.sub("@([0-9]+)","{:s}",match["expr"])
                    parts = [ e.strip() for e in match["p"].split(",") ]
                    expr = fmt.format(*parts)
                    matches = subpattern.finditer(expr)
                    for match in matches:
                        p = match.group("p")
                        if not p in histograms.keys():
                            histograms[p] = {}                        
                        histograms[p][b] = float(match.group("c").replace(" ",""))

                        
