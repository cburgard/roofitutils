# this is the implementation of the ATLAS statistics challenges
#
# Start by setting up your environment, for example using
#
#     setupATLAS
#     lsetup "root 6.22.06-x86_64-centos7-gcc8-opt"
#     lsetup cmake
#     python3 -m venv myenv
#     source myenv/bin/activate
#     python -m pip install pytest
#
# Next, check out and compile RooFitUtils
#
#     git clone ssh://git@gitlab.cern.ch:7999/cburgard/RooFitUtils.git
#     mkdir build
#     cd build
#     cmake ../RooFitUtils
#     make -j4
#     source setup.sh
#     cd ..
#
# Then, run the challenges
#
#     git clone ssh://git@gitlab.cern.ch:7999/cburgard/statchallenges.git
#     cd statchallenges
#     pytest --subjects=../RooFitUtils/python/statchallenges.py -k simplehf

def runfit(filename,workspacename,moreargs={}):
    from RooFitUtils.fit import setup,buildModel,buildMinimizer,fit
    from RooFitUtils.io import result2dict
    args = {"inFileName":filename,"wsName":workspacename,"loglevel":"ERROR"}
    args.update(moreargs)
    setup(args)
    model = buildModel(args)
    minimizer = buildMinimizer(args,model)
    result = fit(args,model,minimizer)
    d = result2dict(result)
    return d

def test_simplehf_poivals(filename,workspacename,dataname):
    return runfit(filename,workspacename,{"hesse":False,"dataName":dataname})["min"]["parameters"]

def test_hww_poivals(filename,workspacename,dataname):
    return runfit(filename,workspacename,{"hesse":False,"dataName":dataname,"poi":["mu_ggF","mu_VBF"]})["min"]["parameters"]

def test_simplehf_poiscan(filename,workspacename,dataname,scanpars):
    result = runfit(filename,workspacename,{"scan":scanpars,"hesse":False,"dataName":dataname})["scans"]["mu"]
    return [ point["nll"] for point in result ]

def test_simplehf_covmat(filename,workspacename,dataname):
    from RooFitUtils.io import dict2mat
    return dict2mat(runfit(filename,workspacename,{"hesse":True,"findSigma":False,"dataName":dataname})["min"]["cov"])

# this helper function is only for debugging
if __name__ == "__main__":
    for elem in dir():
        if elem.startswith("test_simplehf_"):
            try:
#                globals()[elem]("input/simpleHistFactory/workspace.root","Test")
                test_simplehf_covmat("input/simpleHistFactory/workspace.root","Test")
            except err:
                print(err)
                pass

