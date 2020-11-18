import matplotlib
import numpy as np

def parsecsv(infile,delimiter=',',quotechar='|'):
    # requires a csv form where the first column serves as a label 
    import csv
    with open(infile) as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter, quotechar=quotechar)
        header = next(reader, None)
        tabledict = {}
        for row in reader:
            tabledict[row[0]] = {}
            for icol in range(1,len(header)):
                tabledict[row[0]][header[icol]] = float(row[icol])
                
    return tabledict 

def dicttoarr(indict,names=False):
    # convert a dict of dict of a matrix to an array
    # pass names as True if the names of the parameters
    # is to be returned
    matarr = []
    for ipar1 in sorted(indict.keys()):
        tmparr = []
        for ipar2 in sorted(indict.keys()):
            tmparr.append(indict[ipar1][ipar2])
        matarr.append(tmparr)
    if names:
        return sorted(indict.keys()), matarr
    else: return matarr
    
def makecovmat(corarr,uncer):
    # make array of a covariance matrix given a
    # array of correlation matrix and uncertainities
    import math
    covarr, npar = [], len(uncer)
    for ipar1 in range(0,npar):
        val1 = uncer[ipar1]
        tmparr = []
        for ipar2 in range(0,npar):
            val2 = uncer[ipar2]
            tmparr.append(corarr[ipar1][ipar2]*val1*val2)
        covarr.append(tmparr)
        
    return covarr

def arrtotmat(arr):
    # convert array to a TMatrix object
    import ROOT
    
    npar = len(arr)
    tmat = ROOT.TMatrixDSym(npar)
    for ipar1 in range(0,npar):
        for ipar2 in range(0,npar):
            tmat[ipar1][ipar2] = arr[ipar1][ipar2]
    
    return tmat

def getunc(indict,name="asm up"):
    import math
    uncer = []
    for par in sorted(indict.keys()):
        uncer.append(math.sqrt(indict[par][name+"up"]*indict[par][name+"up"]+indict[par][name+"down"]*indict[par][name+"down"]))
    return uncer

summdict = parsecsv('/project/atlas/users/rahulb/project/results/stxssummary_asm.data')
corrdict = parsecsv('/project/atlas/users/rahulb/project/results/correlations-stxs-asimov.data')
cormat = dicttoarr(corrdict)
covmat = makecovmat(cormat,getunc(summdict,name="asm "))

fitRes = mvgpdf.fitTo(dat, ROOT.RooFit.Minimizer('Minuit2', 'Migrad'), ROOT.RooFit.Hesse(True), ROOT.RooFit.Save(True))
import ROOT
cenvals = {}
pois = {}
for i in summdict.keys():
    x = i.replace("mu","x")
    pois[i] = ROOT.RooRealVar(x,x,1,-10,10)
    cenvals[i] = ROOT.RooRealVar(i,i,summdict[i]['asm'])

poilist = ROOT.RooArgList()
poiset = ROOT.RooArgSet()
cvlist = ROOT.RooArgList()
cvset = ROOT.RooArgSet()

for i in cenvals:
    poilist.add(pois[i])
    poiset.add(pois[i]) 
    cvlist.add(cenvals[i])
    cvset.add(cenvals[i])
    
covtmat = ROOT.TMatrixDSym(len(cvlist))
covtmat = arrtotmat(covmat)

mvgpdf = ROOT.RooMultiVarGaussian("mvg","mvg",poilist,cvlist,covtmat)
dat = ROOT.RooDataSet('asimovData_1', '',ROOT.RooArgSet(cvlist))
dat.add(ROOT.RooArgSet(cvlist))

w = ROOT.RooWorkspace("combWS")
mc = ROOT.RooStats.ModelConfig("ModelConfig",w)
mc.SetPdf(mvgpdf)
mc.SetParametersOfInterest(poiset)
mc.SetObservables(cvset)
getattr(w,'import')(mc)
getattr(w,'import')(dat)
w.writeToFile("test.root",True)
