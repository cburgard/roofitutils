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

def getvalues(pois,allsets):
    """convert postfit values and confidence intervals to a hepdata dictionary"""
    from RooFitUtils.util import parsepois
    poinames,poiopts = parsepois(pois)
    row = ["POI"]
    for options,poiset in allsets:
        row.append(options.get("label","default"))
    tab = [row]
    for x in poinames:
        row = [x]
        for options,poiset in allsets:
            if not x in poiset.keys(): continue
            tup = poiset[x]
            if options.get("interval",False):
                row.append(",".join([str(lo)+":"+str(hi) for lo,hi in tup]))
            if options.get("error",True):
                cv,lo,hi = tup
                if isnan(lo) or isnan(hi):
                    if isnan(lo): lo=cv
                    if isnan(hi): hi=cv
                row.append(str(cv-abs(hi))+"/"+str(cv+abs(lo)))                    
            if options.get("point",True):
                try:
                    cvs = [v for v in tup]
                except TypeError:
                    cvs = [tup]
                for cv in cvs:
                    row.append(cv)
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
