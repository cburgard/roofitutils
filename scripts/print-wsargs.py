#!/bin/env python

def main(args):
    import ROOT

    for inputfile in args.input:
        infile = ROOT.TFile.Open(inputfile,"READ")
    
        for k in infile.GetListOfKeys():
            if k.GetClassName()=="RooWorkspace":
                thisws = k.ReadObj()
                break
        if not thisws:
            exit(0)
    
        if args.variables:
            from re import compile
            expressions = [ compile(exp) for exp in args.variables ]
            varnames = [ v.GetName() for v in thisws.allVars() ]
            for var in varnames:
                for exp in expressions:
                    if exp.match(var):
                        print(var)
            

        

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="print some arguments from a workspace")
    parser.add_argument("input",nargs="+",help="workspace")
    parser.add_argument("--variables",nargs="+",help="list of regular expressions")
    args = parser.parse_args()
    main(args)
