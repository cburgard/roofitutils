#!/bin/evn python

import os

def texprep(string):
    import re
    x = ""
    wilcoeffpat = re.compile("(^c\w+)")
    mwilcoeff = wilcoeffpat.search(string)
    if mwilcoeff:
        x = string.replace("c","$c_{") + "}$"
        if "box" in x: x = x.replace("box","\\Box")
    else : x = string.replace("_","\_")
    return x

def writeResult(out,result,writehesse):
    if result.min.nll:
        out.write("Minimization: minNll = ")
        out.write(str(result.min.nll))
        out.write("\n")
        for p in result.parameters:
            out.write("{0:s} = {1:g} - {2:g} + {3:g}\n".format(p.name,p.value,abs(p.errLo),abs(p.errHi)))
    if result.fit and writehesse:
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
        out.write((" ".join(scan.parNames)) + " nll status\n")
        for i in range(0,len(scan.nllValues)):
            out.write((" ".join([ str(scan.parValues[i][j]) for j in range(0,len(scan.parNames)) ]))+" "+str(scan.nllValues[i])+" "+str(scan.fitStatus[i])+"\n")

def collectpoints(points,files,label):
    import glob
    filenames = []
    for expression in files:
        filenames.extend(glob.glob(expression))
    if len(filenames) == 0:
        print("no points found in "+expression)
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
        if x.startswith(parfilter) : 
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
                           parnames.append(texprep(parline))
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
           
def collectresults(scans,results,files,label):
    """collect a set of results files and return the contents as a dictionary"""
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
                    except KeyError:
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
