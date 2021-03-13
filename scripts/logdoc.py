#!/bin/env python
    
def parselines(elements,lines):
    from RooFitUtils.logfile_diagnostics import make_parser
    parser = make_parser()
    
    if args.debugcategory:
        import re
        for k,v in parser.parsers.items():
            if re.match(args.debugcategory,k):
                v.set_debug(True)
    
    from more_itertools import peekable
    logfile = peekable(lines)
    
    parser.parse(elements,logfile)

            
def main(args):
    elements = {}
    with open(args.logfile,"rt") as logfile:
        parselines(elements,logfile)
    for category in elements.keys():
        print("{:s}: {:d}".format(category,len(elements[category])))

    for category in args.printcategories:
        print("the following messages were encountered in category "+category+":")
        for element in elements[category]:
            print(element)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("log file doctor")
    parser.add_argument("logfile")
    parser.add_argument("--print",dest="printcategories",nargs="+",type=str,default=[])
    parser.add_argument("--debug",dest="debugcategory",type=str,default=None)        
    args = parser.parse_args()
    main(args)
