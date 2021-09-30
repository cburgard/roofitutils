#!/bin/env python

def getResult_fromROOT(inResult,poinames,tag):
    from RooFitUtils.io import getFitResult
    from RooFitUtils.util import nodel
    import ROOT
    parts = inResult.split(":")
    result = getFitResult(*parts)
    
    poilist, poiset = ROOT.RooArgList(), ROOT.RooArgSet()
    cvlist, cvset = ROOT.RooArgList(), ROOT.RooArgSet()

    for p in poinames:
        poi = result.floatParsFinal().find(p)
        if not poi:
            raise RuntimeError("unable to find parameter named "+p+" in result!")
        poilist.add(poi)
        poiset.add(poi)

        cv = ROOT.RooRealVar(p+tag.strip(),p+tag.strip(),poi.getVal())
        cv.setConstant(True)
        nodel(cv)
        cvlist.add(cv)
        cvset.add(cv)

        
    return {
        "poilist": poilist,
        "poiset": poiset,
        "cvlist": cvlist,
        "cvset": cvset,
        "covmat": result.reducedCovarianceMatrix(poilist)
    }


    
def makeResult_fromCSV(input,poilist,tag):
    from RooFitUtils.util import getunc, arrtomat, makecovmat, dicttoarr

    inSumm,inCorr = input.split(",")
    
    summdict = parsecsv(inSumm)
    corrdict = parsecsv(inCorr)
    cormat = dicttoarr(corrdict)
    covmat = makecovmat(cormat,getunc(summdict,name=tag))
    
    import ROOT
    cenvals = {}
    pois = {}
    for i in sorted(summdict.keys()):
        x = "poi_"+i
        pois[i] = ROOT.RooRealVar(x,x,1,-inf,inf)
        cenvals[i] = ROOT.RooRealVar(i,i,summdict[i][tag.replace(" ","")])
    
    poilist, poiset = ROOT.RooArgList(), ROOT.RooArgSet()
    cvlist, cvset = ROOT.RooArgList(), ROOT.RooArgSet()
    
    for i in sorted(cenvals):
        poilist.add(pois[i])
        poiset.add(pois[i]) 
        cvlist.add(cenvals[i])
        cvset.add(cenvals[i])
        
    return {
        "poilist": poilist,
        "poiset": poiset,
        "cvlist": cvlist,
        "cvset": cvset,
        "covmat": covmat
    }
    
    
def make_mvg(args):
    if ".csv" in args.inResult:
        result = makeResult(args.inResult,args.pois,args.tag)
    elif ".root" in args.inResult:
        result = getResult_fromROOT(args.inResult,args.pois,args.tag)

    import ROOT
    result["poilist"].Print("v")
    mvgpdf = ROOT.RooMultiVarGaussian("mvg","mvg",result["poilist"],result["cvlist"],result["covmat"])
    dat = ROOT.RooDataSet(args.data, '',result["cvset"])
    dat.add(result["cvset"])
    w = ROOT.RooWorkspace(args.workspace)
    mc = ROOT.RooStats.ModelConfig("ModelConfig",w)
    mc.SetPdf(mvgpdf)
    mc.SetParametersOfInterest(result["poiset"])
    mc.SetSnapshot(result["poiset"])
    mc.SetObservables(result["cvset"])    
    mc.SetGlobalObservables(result["cvset"])
    getattr(w,'import')(mc)
    getattr(w,'import')(dat)
    print("writing workspace to "+args.outWS)
    w.writeToFile(args.outWS,True)
    
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("create multi-variate gaussian RooWorkspace")
    parser.add_argument( "-i","--input"           , type=str,     dest="inResult"                 , help="input file with result", required=True, metavar="path/to/file")
    parser.add_argument( "--pois",nargs="+"     , type=str,     dest="pois"                 , help="list of POI names", default=[])
    parser.add_argument( "--datasetName"         , type=str,     dest="data"                   , help="name of output dataset", default="combData")
    parser.add_argument( "--tag"         , type=str,     dest="tag"                   , help="name of the tag to be used", default="asm ")   
    parser.add_argument( "-o","--output"               , type=str,     dest="outWS"                  , help="path to output workspace"     , required=False,default="out.root")
    parser.add_argument( "--workspace",             type=str,  help="name of output workspace"     , required=False,default="combWS")    

    args = parser.parse_args()

    make_mvg(args)


