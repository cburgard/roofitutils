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

def test_simplehf_poivals(filename,workspacename):
    return runfit(filename,workspacename,{"hesse":False})["min"]["parameters"]

def test_simplehf_poiscan(filename,workspacename,scanpars):
    result = runfit(filename,workspacename,{"scan":scanpars,"hesse":False})["scans"]["mu"]
    return [ point["nll"] for point in result ]

def test_simplehf_covmat(filename,workspacename):
    from RooFitUtils.io import dict2mat
#    constpars = ["gamma_stat_CR_bin_0","gamma_stat_SR_bin_0","gamma_stat_SR_bin_1","gamma_stat_SR_bin_10","gamma_stat_SR_bin_2","gamma_stat_SR_bin_3","gamma_stat_SR_bin_4","gamma_stat_SR_bin_5","gamma_stat_SR_bin_6","gamma_stat_SR_bin_7","gamma_stat_SR_bin_8","gamma_stat_SR_bin_9"]
    return dict2mat(runfit(filename,workspacename,{"hesse":True,"findSigma":False,#"fixParameters":constpars
                                                   })["min"]["cov"])

def test_eft_poivals(filename,workspacename):
    return runfit(filename,workspacename)["min"]["parameters"]

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

