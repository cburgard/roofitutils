
import os
import math

latex_aliases = {}

def texdict(d):
    global latex_aliases
    latex_aliases.update(d)

def texdef(a,b):
    global latex_aliases
    latex_aliases[a] = b

def texify(string,math=False):
    if string in latex_aliases.keys():
        return latex_aliases[string]
    if math:
        return "\\ensuremath{"+string.strip("$")+"}"
    else:
        return string.replace("_","\\_")

def dict2mat(mat,parnames=None):
    if not parnames:
        parnames = sorted(list(set(list(mat.keys()))))
    m = [ [ mat[pi][pj] if pj in mat[pi].keys() else mat[pj][pi] for pi in parnames] for pj in parnames ]
    from numpy import array
    return array(m)

def mat2dict(parameters,mat):
    return {parameters[j].name:{parameters[i].name:mat[i][j] for i in range(0,j+1)} for j in range(0,len(parameters))}

def list2dict(l,key="name"):
    return { v[key]:v for v in l }

def result2dict(result,addcorrmat=True):
    return result

def result2dict_extendedResult(result,addcorrmat=True):
    from RooFitUtils.util import isclose
    mini_info = {"type":result.min.config.MinimizerType(),"strategy":result.min.strategy,"status":result.min.status}
    d = {"MLE":{"minimizer":mini_info,"nll":result.min.nll, "parameters":[]}}
    for p in result.min.parameters:
        par={"name":str(p.name),"val":p.value}
        if isclose(abs(p.errHi),abs(p.errLo)):
            par["err"]=p.errHi
        else:
            par["eUp"]=p.errHi
            par["eDn"]=p.errLo
        d["MLE"]["parameters"].append(par)
    if result.min.fit and addcorrmat:
        npar = len(result.min.parameters)
        d["MLE"]["cov"] = { "parameter_names":[ p.name for p in result.min.parameters ],"matrix":[ [ result.min.fit.covarianceMatrix()[i][j] for j in range(0,npar) ] for i in range(0,npar) ] }
    d["scans"] = []
    for scan in result.scans:
        d["scans"].append({"label":scan.name,"points":[{"nll":scan.nllValues[i],"minimizer":{"status":scan.fitStatus[i]},"parameters":[{"name":scan.parNames[j],"val":scan.parValues[i][j]} for j in range(len(scan.parNames))]} for i in range(len(scan.nllValues))]})
    return d
        
def writeResultJSON(out,result,writecorrmat):
    if type(result) is dict:
        js = result2dict(result,writecorrmat)
    else:
        js = result2dict_extendedResult(result,writecorrmat)
    from datetime import datetime
    now = datetime.now() # current date and time
    timestamp = now.strftime("%Y-%m-%d %H:%M:%S")
    import getpass
    import platform
    import socket
    userstamp = getpass.getuser() + "@" + socket.gethostname() + " " + platform.system() + " " + platform.release()
    import json
    js["stamp"] = timestamp+" "+userstamp
    json.dump(js,out,sort_keys=True,indent=4)

    
def writeResultTxt(out,result,writecorrmat):
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
        out.write(" ".join(map(str,scan.parNames)) + " nll status "+ " ".join(map(str,scan.extraParNames)) + "\n")
        for i in range(0,len(scan.nllValues)):
            out.write(" ".join([ str(scan.parValues[i][j]) for j in range(0,len(scan.parNames)) ]) + " "+str(scan.nllValues[i])+" "+str(scan.fitStatus[i]) + " " + " ".join([ str(scan.extraParValues[i][j]) for j in range(0,len(scan.extraParNames)) ]) + "\n")

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
        from os.path import isfile
        if isfile(filename):
            with open(filename,'r') as infile:
                for line in infile:
                    point = parsedict(line,float)
                    allpoints.append(point)
    points[label] = allpoints

def reduceparams(allpars,parfilter):
    import re
    skimpars = ""
    parpat = re.compile(r"^"+parfilter)
    for x in allpars:
        if x.startswith(parfilter):
            skimpars += " "+ texify(x)
    return skimpars[1:len(skimpars)]

def getFitResult(infilename,key=None):
    import ROOT
    file0 = ROOT.TFile.Open(infilename,"READ")
    if not file0 or not file0.IsOpen():
        raise RuntimeError("unable to open file "+infilename)
    if key:
        fitresult = file0.Get(key)
    else:
        for k in file0.GetListOfKeys():
            if k.GetClassName() == "RooFitResult":
                fitresult = k.ReadObj()
    file0.Close()
    return fitresult

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
    return cov2corr(mat)

def cov2corr(mat):
    from numpy import array
    from math import sqrt
    try:
        return array([ [ mat(i,j) / sqrt(mat(i,i) * mat(j,j)) for i in range(0,mat.GetNcols()) ] for j in range(0,mat.GetNrows()) ])
    except:
        return array([ [ mat[i][j] / sqrt(mat[i][i] * mat[j][j]) for i in range(0,len(mat)) ] for j in range(0,len(mat)) ])        
    
def readtables(infilenames):
    from RooFitUtils.util import isstr,striplatex
    if isstr(infilenames):
        infilenames = [infilenames]
    import re
    begintable = re.compile(r"\\\\begin[ ]*{tabular}[ ]*{(.*)}")
    endtable = re.compile(r"\\\\end[ ]*{tabular}")
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

def collectresult_json(results,filename,label,**kwargs):
    from RooFitUtils.util import isstr
    if not "scans" in results.keys():
        results["scans"] = {}
    scans = results["scans"]
    if not "MLE" in results.keys():
        results["MLE"] = {}
    if not "cov" in results.keys():
        results["cov"] = {}    
    with open(filename,"rt") as infile:
        import json
        js = json.load(infile)
        if "scans" in js.keys():
            for scan in js["scans"]:
                if isstr(scan["label"]):
                    key = tuple(map(str,scan["label"].split(",")))
                else:
                    key = tuple(scan["label"])
                if not key in scans.keys():
                    scans[key] = {}
                scans[key][label] = {}
                for point in scan["points"]:
                    if point["minimizer"]["status"] != 0 and point["minimizer"]["status"] != 1:
                        continue
                    if "filterScans" in kwargs.keys() and kwargs["filterScans"]:
                        pvals = tuple([ p["val"] for p in point["parameters"] if p["name"] in kwargs["filterScans"] ])
                    else:
                        pvals = tuple([ p["val"] for p in point["parameters"]])
                    nll = point["nll"]
                    scans[key][label][pvals] = nll
        if "MLE" in js.keys():
            MLE = results["MLE"]
            for par in js["MLE"]["parameters"]:
                pname = par["name"]
                if not label in MLE.keys():
                    MLE[label] = {}
                MLE[label][pname] = par
            if "cov" in js["MLE"].keys():
                matrix = js["MLE"]["cov"]["matrix"]
                params = js["MLE"]["cov"]["parameter_names"]        
                results["cov"][label] = { params[i]:{params[j]:matrix[i][j] for j in range(len(params))} for i in range(len(params)) }


def collectresult_root(results,filename,label,**kwargs):
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
        if obj.InheritsFrom(ROOT.TTree.Class()):
            from os.path import basename
            from RooFitUtils.util import allkeys
            rdf = ROOT.RDataFrame(obj)
            if "filterScans" in kwargs.keys() and kwargs["filterScans"]:
                pois = kwargs["filterScans"]
                key = tuple(pois)
                nparr = rdf.AsNumpy(columns=pois + ["nll","status"])                
            else:
                key = basename(filename)
                nparr = rdf.AsNumpy()
            if obj.GetName() == "nllscan":
                if not "scans" in results.keys():
                    results["scans"] = []
                scans = results["scans"]            
                scan = { "label":key, "points":[] }
                for pname,values in nparr.items():
                    for i in range(0,len(values)):
                        if len(scan["points"]) <= i:
                            scan["points"].append({})
                        point = scan["points"][i]
                        if pname == "nll":
                            point["nll"] = float(values[i])
                        elif pname == "status":
                            point["minimizer"] = {"status":int(values[i])}
                        else:
                            if not "parameters" in point.keys():
                                point["parameters"] = []
                            point["parameters"].append({"name":pname,"val":float(values[i])})
                scans.append(scan)
            else:
                if not "limits" in results.keys():
                    results["limits"] = {}
                limits = results["limits"]                        
                limit = { k:float(arr[0]) for k,arr in nparr.items() if "upperlimit" in k }
                if "fit_status" in nparr.keys():
                    limit["info"] = {"status":float(nparr["fit_status"])}
                limit["obs"] = { k:float(arr[0]) for k,arr in nparr.items() if "_obs" in k and not "param_" in k and not "hat" in k}
                limit["med"] = { k:float(arr[0]) for k,arr in nparr.items() if "_med" in k and not "param_" in k and not "hat" in k}
                limit["parameters_med"] = { k:float(arr[0]) for k,arr in nparr.items() if "_med" in k and "param_" in k}
                limit["parameters_hat"] = { k:float(arr[0]) for k,arr in nparr.items() if "_hat" in k and "param_" in k}            
                limit["parameters_const"] = { k:float(arr[0]) for k,arr in nparr.items() if not k in allkeys(limit) }
                limits[key] = {label:limit}


def collectresult_txt(results,filename,label,**kwargs):
    if not "scans" in results.keys():
        results["scans"] = {}
    scans = results["scans"]
    if not "MLE" in results.keys():
        results["MLE"] = {}
    MLEs = results["MLE"]
    import re
    parpat = re.compile(r"([a-zA-Z0-9_.-]+)[ ]*=[ ]*([0-9.naife+-]+)[ ]*-[ ]*([0-9.naife-]+)[ ]*\+[ ]*([0-9.naife+-]+)[ ]*")
    parpat_legacy = re.compile(r"[ \t]*([a-zA-Z0-9_.-]+)[ \t=]+([0-9.naife+-]+)[ \t]*\+/-[ \t]*([0-9.naife+-]+).*")
    nllpat = re.compile(r"Minimization:[ ]*minNll[ ]*=[ ]*([0-9.naife+-]+)")
    with open(filename,'r') as infile:
        lines = [ line for line in infile ]
        for lineno in range(0,len(lines)):
            line = lines[lineno]
            if ".py" in line:
                print("skipping file "+filename+", probably a job definition file")
                break
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
                    if not pname in MLEs.keys():
                        MLEs[pname] = {}
                    MLEs[pname][label] = result
                elif match_legacy:
                    pname,cv,e = match_legacy.group(1).strip(),match_legacy.group(2),match_legacy.group(3)
                    result = (float(cv),-float(e),float(e))
                    if not pname in MLEs.keys():
                        MLEs[pname] = {}
                    MLEs[pname][label] = result                            
                elif "nll" in parts:
                    npars = parts.index("nll")
                    scanps = parts[0:npars]
                    key = tuple(scanps)
                    if key not in scans.keys():
                        scans[key] = {}
                    if label not in scans[key].keys():
                        scans[key][label] = {}
                else:
                    pvals   = tuple([float(parts[i]) for i in range(0,npars)])
                    nllval = float(parts[npars])
                    if int(parts[npars+1]) == 0: # check if fit ended with status 0
                        scans[key][label][pvals] = nllval
            except UnboundLocalError as err:
                print("unable to parse line '"+line.strip()+"' in file '"+filename+"', trying to read as nll scan, but no parameter list found. ")
            except (KeyError,ValueError) as err:
                if nllmatch:
                    print("unable to parse line '"+line.strip()+"' in file '"+filename+"', attempted to parse as nll value. "+str(err))
                elif match:
                    print("unable to parse line '"+line.strip()+"' in file '"+filename+"', attempted to parse as parameter value. "+str(+err))
                else:
                    print("unable to parse line '"+line.strip()+"' in file '"+filename+"', "+str(err))
                continue
    for p in scans.keys():
        if p in MLEs.keys():
            scans[p][label][MLEs[p][label][0]]=minnll;
    
def file_starts_with(fname,sequence):
    with open(fname,"rt") as thefile:
        return next(thefile).startswith(sequence)
        
    
def collectresult(results,filename,label,**kwargs):
    import re
    if os.path.isfile(filename):
        if filename.endswith(".json") or file_starts_with(filename,"{"):
            collectresult_json(results,filename,label,**kwargs)
        elif filename.endswith(".root"):
            collectresult_root(results,filename,label,**kwargs)
        else:
            collectresult_txt(results,filename,label,**kwargs)

def collectfilenames(files):
    if isinstance(files, str): files = [files]
    import glob
    filenames = []
    for expression in files:
        filenames.extend(glob.glob(expression))
    if len(filenames) == 0:
        raise RuntimeError("no files found matching "+",".join(files))
    return filenames

def collectimpacts(rankings,files,poiname):
    nan = float("nan")
    filenames = collectfilenames(files)
    import re
    fpat = re.compile(r"[^.]*\.([^/^.]*)\.([^.]*)\.json")
    from os.path import basename
    results = {}
    allvars = []
    for filename in filenames:
        fname = basename(filename)
        if "impact.nominal" in fname:
            collectresult(results,filename,"nominal")
        else:
            match = fpat.match(fname)
            if not match: continue
            npname = match.group(1)
            allvars.append(npname)
            direction = match.group(2)
            collectresult(results,filename,(npname,direction))
    if not 'MLE' in results.keys():
        raise RuntimeError("did not find any impact results!")
    fits = results["MLE"]
    if not "nominal" in fits.keys():
        raise RuntimeError("cannot create ranking without nominal result!")
    try:
        nomval = fits["nominal"][poiname]["val"]
    except KeyError:
        raise RuntimeError("unable to retrieve result for POI '{:s}', available POIs are: {:s}".format(poiname,",".join(fits["nominal"].keys())))
    if not poiname in rankings.keys():
        rankings[poiname] = {}
    for par in allvars:
        if par == poiname: continue
        try:
            upval = fits[(par,"up")][poiname]["val"]
        except KeyError:
            upval = nan
        try:
            dnval = fits[(par,"dn")][poiname]["val"]
        except KeyError:
            dnval = nan
        rankings[poiname][par] = (upval-nomval,dnval-nomval)
    return allvars

def collectbreakdowns(rankings,files,poiname):
    filenames = collectfilenames(files)
    import re
    fpat = re.compile(r"breakdown\.([^/^.]*)\.json")
    from os.path import basename
    results = {}
    allvars = []
    for filename in filenames:
        fname = basename(filename)
        match = fpat.match(fname)
        if not match: continue
        name = match.group(1)
        collectresult(results,filename,name)
    if not "MLE" in results.keys():
        raise RuntimeError("unable to retrieve result")
        return allvars
    fits = results["MLE"]
    if not "nominal" in fits.keys():
        raise RuntimeError("cannot create ranking without nominal result!")
    if not poiname in rankings.keys():
        rankings[poiname] = {}
    for key in fits.keys():
        if key == "nominal": continue
        upvar = + math.sqrt(pow(fits["nominal"][poiname]["eUp"],2) - pow(fits[key][poiname]["eUp"],2))
        dnvar = - math.sqrt(pow(fits["nominal"][poiname]["eDn"],2) - pow(fits[key][poiname]["eDn"],2))
        rankings[poiname][key] = (upvar,dnvar)
        allvars.append(key)
    return allvars
        
    
def collectresults(results,files,label="default",**kwargs):
    """collect a set of results files and return the contents as a dictionary"""
    filenames = collectfilenames(files)
    for filename in filenames:
        collectresult(results,filename,label,**kwargs)

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
        pattern = re.compile(cfg.get("pattern",r".*_(?P<p>[A-z0-9]*)_.*"))
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
        subpattern = re.compile(cfg.get("subpattern",r"(?P<c>[+-]?[ ]*[0-9.]+)[ *]*(?P<p>[A-z][A-z0-9]+)"))
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

                        
def readparameter(elem,style={}):
    from RooFitUtils.util import isdict
    cvs = []
    intervals = []
    if isdict(elem):
        try:
            cvs = [ v for v in elem["val"] ]
        except TypeError:
            cvs = [elem["val"]]
        if "intervals" in elem.keys():
            intervals = elem["intervals"]
        elif "eDn" in elem.keys() and "eUp" in elem.keys():
            intervals = [[elem["val"]-abs(elem["eDn"]),elem["val"]+abs(elem["eUp"])]]
        elif "err" in elem.keys():
            intervals = [[elem["val"]-abs(elem["err"]),elem["val"]+abs(elem["err"])]]
    elif style.get("interval",False):
        intervals = [(lo,hi) for lo,hi in elem]
    elif style.get("error",True):
        cv,lo,hi = elem
        cvs = [cv]
        intervals = [[cv-abs(lo),cv+abs(hi)]]
    elif style.get("point",True):
        cvs = [elem[0]]
    return cvs,intervals

