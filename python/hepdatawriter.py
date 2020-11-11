def writeyaml(outfilename,thedict):
    import yaml
    yml = yaml.dump(thedict, sort_keys=True)
    with open(outfilename,"wt") as outfile:
        outfile.write(yml)

def getmatrix(xcoords,ycoords,allvalues,label):       
    """convert a correlation matrix to a hepdata dictionary"""
    mat = {"independent_variables":[
        {"header": {"name": "parameter_1"}, "values":[]},
        {"header": {"name": "parameter_2"}, "values":[]},
    ],
              "dependent_variables":[
                  {"header": {"name":label}, "values":[]}
              ]
              }
    for i in range(0,len(xcoords)):
        for j in range(0,len(ycoords)):
            mat["independent_variables"][0]["values"].append(xcoords[i])
            mat["independent_variables"][1]["values"].append(ycoords[j])
            mat["dependent_variables"][0]["values"].append(float(allvalues[j][i]))
    return mat

def getvalues(pois,allsets):
    """convert postfit values and confidence intervals to a hepdata dictionary"""
    from RooFitUtils.util import parsepois
    from math import isnan    
    poinames = parsepois(pois,options=False)    
    values = {"independent_variables":[{"header":{"name":"POI"},"values":poinames}],"dependent_variables":[]}
    for x in poinames:
        for options,poiset in allsets:
            label = options.get("label","default")
            variable = {"header":{"name":x,"label":label},"values":[],"qualifiers":[{"POI":x,"label":label}]}
            if not x in poiset.keys(): continue
            tup = poiset[x]
            if options.get("interval",False):
                for lo,hi in tup:
                    variable["values"].append({"low":lo,"high":hi})
            if options.get("error",True):
                cv,lo,hi = tup
                if isnan(lo) or isnan(hi):
                    if isnan(lo): lo=cv
                    if isnan(hi): hi=cv
                variable["values"].append({"high":cv+abs(hi),"low":cv-abs(lo)})
            if options.get("point",True):
                if options.get("error",True):
                    cvs = [tup[0]]
                else:
                    try:
                        cvs = [v for v in tup]
                    except TypeError:
                        cvs = [tup]                
                for cv in cvs:
                    variable["values"].append({"value":cv})
            values["dependent_variables"].append(variable)
    return values

def writematrix(xcoords,ycoords,allvalues,outfilename,label="correlation"):
    """write a correlation matrix to a hepdata yml file"""    
    yml = getmatrix(xcoords,ycoords,allvalues,label)
    writeyaml(outfilename,yml)

def writepois(pois,allsets,outfilename):
    """write postfit values and confidence intervals of the POIs to a hepdata yml file"""    
    yml = getvalues(pois,allsets)
    writeyaml(outfilename,yml)
