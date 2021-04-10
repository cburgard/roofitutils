from RooFitUtils.util import nodel
from ROOT import RooArgList,RooArgSet,RooRealVar,TMatrixDSym,RooMultiVarGaussian,RooFit

def makeasimov(pdf,muvals,xvals):
    import ROOT
    ws = ROOT.RooWorkspace("ws")
    ws.Import(mvg)
    mc = ROOT.RooStats.ModelConfig(ws)
    mc.SetObservables(xvals)
    mc.SetParametersOfInterest(muvals)
    mc.SetPdf(ws.factory("ProdPdf::mainPdf("+pdf.GetName()+")"))
    data = ROOT.RooStats.AsymptoticCalculator.MakeAsimovData(mc,muvals,ROOT.RooArgSet())
    return data

def comparematrix(m1,m2,scale1=1,scale2=1,tolerance=1e-3):
    for i in range(m1.GetNrows()):
        for j in range(m1.GetNcols()):
            if not abs(scale1*m1[i][j] - scale2*m2[i][j]) < tolerance:
                raise RuntimeError("assertion error: {:g} != {:g}".format(scale1*m1[i][j],scale2*m2[i][j]))

xvals = RooArgList()
muvals = RooArgList()

for i in range(0,3):
    mu = RooRealVar("mu_{:d}".format(i),"mu_{:d}".format(i),1.,-10,10)
    nodel(mu)
    muvals.add(mu)
    x  = RooRealVar("x_{:d}".format(i),"x_{:d}".format(i),1.)
    nodel(x)
    xvals.add(x)

variables = RooArgSet()
variables.add(xvals)

covmat = TMatrixDSym(3)
covmat[0][0] = 1.
covmat[1][1] = 2.
covmat[2][2] = 3.
covmat[0][1] = 0.2
covmat[1][0] = 0.2
covmat[1][2] = 0.6
covmat[2][1] = 0.6
covmat[2][0] = 0.3
covmat[0][2] = 0.3

nevents = 1000.
mvg = RooMultiVarGaussian("mvg","mvg",xvals,muvals,covmat)
data = mvg.generate(variables,nevents,False)

theresult = mvg.fitTo(data,RooFit.Save())
theresultcov = theresult.covarianceMatrix()

comparematrix(theresultcov,covmat,1000)

import ROOT
mini = ROOT.RooFitUtils.ExtendedMinimizer("mini",mvg,data)
mini.minimize(ROOT.RooFit.Hesse(True))
mini.minimize()
minimum = mini.getResult()

cov = minimum.min.hesse.Invert()
comparematrix(cov,covmat,1000)



