def runfit(filename,workspacename,moreargs={}):
    from RooFitUtils.fit import setup,buildModel,buildMinimizer,fit
    from RooFitUtils.io import result2dict
    args = {"inFileName":filename,"wsName":workspacename}
    args.update(moreargs)
    setup(args)
    model = buildModel(args)
    minimizer = buildMinimizer(args,model)
    result = fit(args,model,minimizer)
    return result2dict(result)

def test_poivals(filename,workspacename):
    return runfit(filename,workspacename)["min"]["parameters"]

def test_poiscan(filename,workspacename,scanpars):
    result = runfit(filename,workspacename,{"scan":scanpars})["scans"]["mu"]
    return [ point["nll"] for point in result ]

def test_covmat(filename,workspacename):
    from RooFitUtils.io import dict2mat
    return dict2mat(runfit(filename,workspacename,{"correlationMatrix":True})["min"]["cov"])

