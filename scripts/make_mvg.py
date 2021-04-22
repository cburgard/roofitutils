from RooFitUtils.util import getunc, arrtomat, makecovmat, dicttoarr

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("create multi-variate gaussian RooWorkspace")
    arglist = []
    arglist.append(parser.add_argument( "--inputCorr"           , type=str,     dest="inCorr"                 , help="csv file of correlation matrix", required=True, metavar="path/to/file"))
    arglist.append(parser.add_argument( "--inputSummary"        , type=str,     dest="inSumm"                 , help="csv file of parameter results", required=True, metavar="path/to/file"))
    arglist.append(parser.add_argument( "--datasetName"         , type=str,     dest="data"                   , help="name of output dataset", required=False, metavar="asimovData_1"))
    arglist.append(parser.add_argument( "--tag"                 , type=str,     dest="tag"                    , help="tag", required=False,default="asm "))
    arglist.append(parser.add_argument( "--outWS"               , type=str,     dest="outWS"                  , help="path to output workspace"     , required=False,default="out.root"))

    args = parser.parse_args()

    summdict = parsecsv(args.inSumm)
    corrdict = parsecsv(args.inCorr)
    cormat = dicttoarr(corrdict)
    covmat = makecovmat(cormat,getunc(summdict,name=args.tag))
    
    import ROOT
    cenvals = {}
    pois = {}
    for i in sorted(summdict.keys()):
        x = "poi_"+i
        pois[i] = ROOT.RooRealVar(x,x,1,-inf,inf)
        cenvals[i] = ROOT.RooRealVar(i,i,summdict[i][args.tag.replace(" ","")])
    
    poilist, poiset = ROOT.RooArgList(), ROOT.RooArgSet()
    cvlist, cvset = ROOT.RooArgList(), ROOT.RooArgSet()
    
    for i in sorted(cenvals):
        poilist.add(pois[i])
        poiset.add(pois[i]) 
        cvlist.add(cenvals[i])
        cvset.add(cenvals[i])
        
    covtmat = ROOT.TMatrixDSym(len(cvlist))
    covtmat = arrtotmat(covmat)
    
    mvgpdf = ROOT.RooMultiVarGaussian("mvg","mvg",poilist,cvlist,covtmat)
    dat = ROOT.RooDataSet(args.data, '',ROOT.RooArgSet(cvlist))
    dat.add(ROOT.RooArgSet(cvlist))
    w = ROOT.RooWorkspace("combWS")
    mc = ROOT.RooStats.ModelConfig("ModelConfig",w)
    mc.SetPdf(mvgpdf)
    mc.SetParametersOfInterest(poiset)
    mc.SetSnapshot(poiset)
    mc.SetObservables(cvset)
    getattr(w,'import')(mc)
    getattr(w,'import')(dat)
    w.writeToFile(args.outWS,True)
