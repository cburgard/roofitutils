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


def writematrix(xcoords,ycoords,allvalues,outfilename,label="correlation"):
    """write a correlation matrix to a hepdata yml file"""    
    yml = getmatrix(xcoords,ycoords,allvalues,label)
    writeyaml(outfilename,yml)
