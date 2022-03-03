#!/bin/env python

def main(args):
    from RooFitUtils.io import collectresults,writeResultJSON
    results = {}
    collectresults(results,args.input,args.label,filterScans=args.pois)
    with open(args.output,"wt") as out:
        writeResultJSON(out,results,True)

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser("convert results to JSON format")
    parser.add_argument('--label',type=str,default="converted")
    parser.add_argument('-i','--input',type=str,nargs="+",help="files with input information",required=True)
    parser.add_argument("-o","--output",type=str,help="output file name",default="result.json",required=True)
    parser.add_argument('--pois',type=str,nargs="+",help="list of POIs to use")    

    args = parser.parse_args()
    main(args)
