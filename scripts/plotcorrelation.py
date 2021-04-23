#!/usr/bin/env python


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser("plot a correlation matrix")
    parser.add_argument('-i','--input',nargs="+",help="text files with the input information",required=True)
    parser.add_argument("--atlas",type=str,help="ATLAS plot label, will enable ATLAS style if used",required=False,default=None)
    parser.add_argument("-o", "--output",type=str,help="output file name",default="scan.tex",required=True)
    parser.add_argument("--filter",type=str,help="filter for parameters in correlation matrix",default=".*",required=False)
    args = parser.parse_args()

    results  = {}

    from RooFitUtils.io import collectresults
    collectresults(results,args.input,label="default")

    d = results["cov"]["default"]
    import re
    parfilter = re.compile(args.filter)
    parnames = [ k for k in d.keys() if parfilter.match(k) ]
    from RooFitUtils.io import dict2mat,cov2corr
    matrix = cov2corr(dict2mat(d,parnames))

    with open(args.output,"w") as outfile:
        from RooFitUtils.pgfplotter import writecorrmatrix
        writecorrmatrix(args.atlas,parnames,matrix,args.output)
