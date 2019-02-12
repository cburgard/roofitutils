#!/bin/evn python

import os

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


def getvals(d,nllmin):
    """transform the values of a curvey by sorting them in x and subtracting a fixed value in y"""
    xvals = sorted(d.keys())
    yvals = [ max(d[k] - nllmin,0) for k in xvals ]
    return list(zip(xvals,yvals))
        

def parsedict(s):
    """parse a string of the format "a=b,c=d" into a dictionary"""
    d = {}
    for kv in s.split(","):
        k,v = kv.split("=")
        d[k] = v
    return d
