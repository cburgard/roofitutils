#!/bin/evn python

import os

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
            out.write(matrix.GetYaxis().GetBinLabel(i+1)+" ")
        out.write("\n")
        for i in range (0,matrix.GetNbinsX()):
            out.write(matrix.GetYaxis().GetBinLabel(i+1)+" ")
            for j in range(0,matrix.GetNbinsY()):
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
                    
def collectresults(scans,results,files,label):
    """collect a set of results files and return the contents as a dictionary"""
    import re
    parpat = re.compile("([a-zA-Z0-9_.a-]+)[ ]*=[ ]*([0-9.naife-]+)[ ]*-[ ]*([0-9.naife-]+)[ ]*\+[ ]*([0-9.naife-]+)[ ]*")
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
#		print lines
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
			    status = int(parts[-1])
		            #if -0.001< pvals[1] < 0.001 and -0.001< pvals[0] < 0.001: 
			    #    print str(pvals)+str(nllval)
#			    if pvals[1] < 0.01 and pvals[0] < 0: print nllval
#			    if pvals[1] < 0.01 and pvals[0] < 0: print status
		            if status == 0: 
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
