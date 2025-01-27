def writecsv(outfilename,thetable):
    import csv
    with open(outfilename, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',',quotechar='"',quoting=csv.QUOTE_MINIMAL)
        for row in thetable:
            writer.writerow(row)

def getmatrix(xcoords,ycoords,allvalues):       
    """convert a correlation matrix to a csv table"""
    tab = [["parameter"] + xcoords]
    for j in reversed(range(0,len(ycoords))):
        row = []
        row.append(ycoords[j])
        for i in range(0,len(xcoords)):
            row.append(float(allvalues[j][i]))
        tab.append(row)
    return tab

def dict2matrix(allvalues):       
    """convert a nested dict to a csv table"""
    from RooFitUtils.util import flattened
    allkeys = sorted(list(allvalues.keys()))
    allsubkeys = sorted(list(set(flattened([d.keys() for d in allvalues.values()]))))
    tab = [["parameter"] + allsubkeys]
    for k1 in allkeys:
        row = [k1]
        for k2 in allsubkeys:
            row.append(float(allvalues[k1].get(k2,0)))
        tab.append(row)
    return tab

def getvalues(pois,allsets):
    """convert postfit values and confidence intervals to a hepdata dictionary"""
    from math import isnan    
    from RooFitUtils.util import parsepois
    poinames = parsepois(pois,options=False)
    row = ["POI"]
    for options,poiset in allsets:
        if options.get("point",True) and options.get("error",True):
            row.append(options.get("label","default"))            
            row.append(options.get("label","default")+" up")            
            row.append(options.get("label","default")+" down")            
        elif options.get("point",True):
            row.append(options.get("label","default"))            
        elif options.get("error",True):
            row.append(options.get("label","default")+" up")            
            row.append(options.get("label","default")+" down")                            
    tab = [row]
    for x in poinames:
        row = [x]
        for options,poiset in allsets:
            if not x in poiset.keys(): continue
            tup = poiset[x]
            if options.get("point",True):
                if options.get("error",True):
                    cvs = [tup[0]]
                else:
                    try:
                        cvs = [v for v in tup]
                    except TypeError:
                        cvs = [tup]                
                for cv in cvs:
                    row.append(cv)
            if options.get("interval",False):
                row.append(",".join([str(lo)+":"+str(hi) for lo,hi in tup]))
            if options.get("error",True):
                cv,lo,hi = tup
                if isnan(lo) or isnan(hi):
                    if isnan(lo): lo=hi
                    if isnan(hi): hi=lo
                row.append("+"+str(abs(hi)))
                row.append("-"+str(abs(lo)))
        tab.append(row)
    return tab

def writematrix(xcoords,ycoords,allvalues,outfilename):
    """write a correlation matrix to a csv file"""    
    tab = getmatrix(xcoords,ycoords,allvalues)
    writecsv(outfilename,tab)

def writepois(pois,allsets,outfilename):
    """write postfit values and confidence intervals of the POIs to a csv file"""    
    tab = getvalues(pois,allsets)
    writecsv(outfilename,tab)
