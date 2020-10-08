
import os

def texprep(string):
    import re
    x = ""
    wilcoeffpat = re.compile("(^c\w+)")
    mwilcoeff = wilcoeffpat.search(string)
    if mwilcoeff:
        x = string.replace("c","$c_{") + "}$"
        if "box" in x: x = x.replace("box","\\Box")
    else : x = string
    return x

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
        exit(0)
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
            skimpars += " "+ texprep(x)
    return skimpars[1:len(skimpars)]

def collectcorrelationmatrix(parnames,results,filename,parfilter,label):
    import re
    ncorrpat = re.compile("Correlations[ ](\w+)")
    import glob
    ncorr = 0
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
                        for x in allpars: parnames.append(texprep(x))
                    pattern = ""
                    for x in range(0, ncorr):
                        pattern += "([-]*\d+.\d+)[ ]"
                    rowpat = re.compile("([a-zA-Z0-9_.]+)[ ]([a-zA-Z0-9_.]+)[ ]"+"([-]*\d+.\d+)")
                    for lineno in range(lineno,len(lines)):
                        rowmatch = rowpat.match(lines[lineno])
                        if rowmatch:
                            par1 = (rowmatch.group(1),rowmatch.group(2))
                            if (len(parfilter) != 0 and par1[0].startswith(parfilter) and par1[1].startswith(parfilter)) or len(parfilter) == 0:
                                results.append( texprep((par1[0])) +" "+ texprep((par1[1]))+" "+rowmatch.group(3))
    parnames = (parnames[0])

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
                        match_legacy = parpat_legacy.match(line)
                        if nllmatch:
                            minnll = float(nllmatch.group(1))
                        elif match or "BkgTheory" in line: # todo: cleanup!
                            pname,cv,ed,eu = match.group(1).strip(),match.group(2),match.group(3),match.group(4)
                            result = (float(cv),float(ed),float(eu))
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
