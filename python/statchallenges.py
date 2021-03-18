# this is the implementation of the ATLAS statistics challenges
#
# First, check out and compile RooFitUtils
#
#     git clone ssh://git@gitlab.cern.ch:7999/cburgard/RooFitUtils.git
#     mkdir build
#     cd build
#     cmake ..
#     make -j4
#     source setup.sh
#     cde ..
#
# Then, run the challenges
#
#     git clone ssh://git@gitlab.cern.ch:7999/cburgard/statchallenges.git
#     cd statchallenges
#     pytest --subjects=../RooFitUtils/python/statchallenges.py -k simplehf

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

def test_simplehf_poivals(filename,workspacename):
    return runfit(filename,workspacename)["min"]["parameters"]

def test_simplehf_poiscan(filename,workspacename,scanpars):
    result = runfit(filename,workspacename,{"scan":scanpars})["scans"]["mu"]
    return [ point["nll"] for point in result ]

def test_simplehf_covmat(filename,workspacename):
    from RooFitUtils.io import dict2mat
    return dict2mat(runfit(filename,workspacename,{"correlationMatrix":True})["min"]["cov"])

def test_eft_poivals(filename,workspacename):
    return runfit(filename,workspacename)["min"]["parameters"]
