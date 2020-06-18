# RooFitUtils

This package is a collection of tools for advanced workspace
manipulation and fitting using RooFit. Among others, it includes

 * ExtendedMinimizer and ExtendedModel, advanced fitting tools originally designed by Stefan Gadatsch
 * editws, workspace manipulation functions designed by Tim Adye.
 * guessCorrelations, a helper script able to intelligently identify nuisance parameters when combining workspaces.

## Setup

This package supports setup within RootCore and CMake based ASG
releases as well as standalone compilation with ROOT.

In order to compile with CMake, type

    mkdir build
    cd build
    cmake .. -DCMAKE_MODULE_PATH=$ROOTSYS/cmake/
    make -j4
    cd ..
    source setup.sh

Now, you are ready to use the package. Don't forget to 

    source setup.sh

every time you create a new shell.

It is recommended to have python 2.7 or higher.
scipy and scikit-image packages are required 
for plotting likelihood results. 

First make sure you have pip installed

    curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
    python get-pip.py

and then install the packages by typing

    pip install --user scipy
    pip install --user scikit-image

## Usage

A fitting font-end are provided in the form of a python script 

    scripts/fit.py

which provides extensive help with the `--help` command line option. The results of likelihood scans obtained with this script can be plotted using 

    scripts/plotscan.py


## Pruning

As an example, you can try the pruning code on a Higgs 5XS workspace. The directory `/afs/cern.ch/user/r/rabalasu/public/test_prune` contains the 5XS 80ifb workspace (pre-fit and post-fit) and the corresponding hesse and fitresults as root files. To test the pruning code, copy it's contents to `test`.
The first step would be to obtain the following inputs for the pruning code, the postfit workspace, hesse matrix, and the fitresult. You can skip this step for the example as they
are already available. 

    python scripts/fit.py --input test/WS-Comb-5XS_80ifb.root --output test/WS-Comb-5XS.txt --workspace combWS --data combData —-no-findSigma —-hesse --writeResult —-make-snapshots —-write-workspace prune_test/WS-Comb-5XS_80ifb_postFit.root

 To perform the prune run the `order_NPs.py` script. Running the order_NPs script with `--writeSubmit` options creates a txt file whose lines split the rank-finding into multiple jobs.

    python scripts/order_NPs.py --pois r_ggF,r_VBF,r_WH,r_ZH,r_ttH --hesse test/WS-Comb-5XS_80ifb_hesse.root --fitResult test/WS-Comb-5XS_80ifb_fitresult.root --writeSubmit test/joblines_5XS.txt --jobTime 10 --output test/order_NPs.txt

This creates the joblines file `test/joblines_5XS.txt`. You can run these lines in parallel which will write out the nuisance parameters and their ranks.
Once all of them have finished, you will see the ranks for the nuisance parameters witten in §`test/order_NPs_*.txt`. You can then run the 
pruning procedure with the following line,

    python scripts/prune_NPs.py --input test/WS-Comb-5XS_80ifb_postFit.root --hesse test/WS-Comb-5XS_80ifb_hesse.root --fitResult test/WS-Comb-5XS_80ifb_fitresult.root --pois r_ggF,r_VBF,r_WH,r_ZH,r_ttH --order test/order_*.txt --output test/WS-Comb-5XS_80ifb_postFit_prune.root 2>&1 | tee test/prune_NPs_WS-Comb-5XS_80ifb.log
    
This creates a workspace `WS-Comb-5XS_80ifb_postFit_prune.root` that contain the snapshots `prune_combData_1pct`,`prune_combData_5pct`, and `prune_combData_10pct` 
that contains the nuisance parameters identified for their respective thresholds set to constant. The snapshot can be loaded before subsequent fits to fix the pruned nuisance parameters.
