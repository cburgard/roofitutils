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




