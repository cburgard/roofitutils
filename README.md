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

## Basic Usage

A fitting font-end are provided in the form of a python script 

    scripts/fit.py

which provides extensive help with the `--help` command line option. 

The results of likelihood scans obtained with this script can be plotted using 

    scripts/plotscan.py

The results of individual fits (parameter pulls) can be with this script can be plotted using 

    scripts/plotpulls.py

There resulting files are `.tex` files that encode the graphics using `TikZ` and can be rendered with `pdflatex` or `xelatex`.

You can investigate log files with `scripts/logdoc.py`.

    
## Advanced Usage

This section includes some examples on how this package can be used to perform some common tasks.

### Liklihood scan

Create job definitions with

    python scripts/fit.py --input workspace.root --fit --data asimovData --scan mu 100 0 5 --writeSubmit scan
    
Use some batch submission of your choice to submit the jobs defined in 

    scan/jobs.txt

After your jobs finish, plot a simple, 1D likelihood scan with

    scripts/plotscan.py -i red scans/*.txt -o my1dscan.tex --atlas Preliminary --poi mu


### Ranking

Create breakdown job definitions with

    python scripts/fit.py --input workspace.root --fit --data asimovData --breakdown "ATLAS_*" --breakdown "theo_*" --breakdown "gamma_*" --writeSubmit breakdown

This will _group_ parameters together, one group per `--breakdown`
argument, and evaluate their combined influence of the POI(s) by
calculating the difference in quadrature between a fit with these
parameters set constant and the full fit.

If you are instead interested in evaluating a ranking by fixing the
nuisance parameters to their +/-1 sigma values and re-fitting to
evaluate the influence on the central value of the POI, you can use
the `impact` mechanism. Similarly, you can create impact job
definitions with

    python scripts/fit.py --input workspace.root --fit --data asimovData --impacts "ATLAS_*" "theo_*"  "gamma_*" --writeSubmit impacts
    
If you want to derive postfit impacts, you need to run in two steps like this:

    python scripts/fit.py -i workspace.root --data asimovData --poi mu --fit --write-workspace postfit.root
    python scripts/fit.py -i postfit.root --poi mu --data asimovData --impacts 'norm_bkg' 'gamma_stat_*' --writeSubmit impacts 
    
Use some batch submission of your choice to submit the jobs defined in 

    breakdown/jobs.txt
    impacts/job.txt

After your jobs finish, plot a full breakdown and impact plot including parameter pulls with

    scripts/plotpulls.py -i red "impacts/*nominal*.json" --impacts mu "impacts/*.json" --range -2 2 --scaleimpacts 5 --numbers --atlas True --output impacts.tex

### Pruning

As an example, you can try the pruning code on a Higgs 5XS
workspace. The directory
`/afs/cern.ch/user/r/rabalasu/public/test_prune` contains the 5XS
80ifb workspace (pre-fit and post-fit) and the corresponding hesse and
fitresults as root files. To test the pruning code, copy its contents
to `test`.  The first step would be to obtain the following inputs for
the pruning code, the postfit workspace, and the fitresult. If the
workspace contains more POIs than those that should be profiled, make sure
to specify only those that should be profiled. 
You can skip this step for the example as they are already
available.

    python scripts/fit.py --poi r_ggF r_VBF r_WH r_ZH r_ttH --input test/WS-Comb-5XS_80ifb.root --output test/WS-Comb-5XS.txt --workspace combWS --data combData --no-findSigma --hesse --writeResult --make-snapshots --write-workspace prune_test/WS-Comb-5XS_80ifb_postFit.root

 To perform the prune run the `order_NPs.py` script. Running the
 order_NPs script with `--writeSubmit` options creates a txt file
 whose lines split the rank-finding into multiple jobs.

    python scripts/order_NPs.py --pois r_ggF,r_VBF,r_WH,r_ZH,r_ttH --fitResult test/WS-Comb-5XS_80ifb_fitresult.root --writeSubmit test/joblines_5XS.txt --jobTime 10 --output test/order_NPs.txt

This creates the joblines file `test/joblines_5XS.txt`. You can run
these lines in parallel which will write out the nuisance parameters
and their ranks.  Once all of them have finished, you will see the
ranks for the nuisance parameters witten in
`test/order_NPs_*.txt`. You can then run the pruning procedure with
the following line,

    python scripts/prune_NPs.py --input test/WS-Comb-5XS_80ifb_postFit.root --fitResult test/WS-Comb-5XS_80ifb_fitresult.root --pois r_ggF,r_VBF,r_WH,r_ZH,r_ttH --snapshot-name prune_combData_1pct --order test/order_*.txt --output test/WS-Comb-5XS_80ifb_postFit_prune.root --write-list test/pruned_NPs.txt

This creates a workspace `WS-Comb-5XS_80ifb_postFit_prune.root` that
contain the snapshots `prune_combData_1pct` that contains the nuisance
parameters identified for their respective thresholds set to 1% of max.
change in the hesse error of any POI. The snapshot can be loaded before
subsequent fits to fix the pruned nuisance parameters.
