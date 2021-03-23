#!/bin/env python
    
def parse(messages,lines):
    from RooFitUtils.logfile_diagnostics import make_parser
    parser = make_parser()
    
    if args.debugcategory:
        import re
        for k,v in parser.parsers.items():
            if re.match(args.debugcategory,k):
                v.set_debug(True)
    
    from more_itertools import peekable
    logfile = peekable(lines)
    
    parser.parse(messages,logfile)

def printdict(d,indent=0):
    for k,v in d.items():
        if isinstance(v,dict):
            print(" "*4*indent + k)
            printdict(v,indent+1)            
        if isinstance(v,list):
            for e in v:
                print(" "*4*indent + k)            
                printdict(e,indent+1)
        else:
            print(" "*4*indent + k + " " + v)

def printsummary(messages):
    print("Summary of messages")
    for category in messages.keys():
        print("  {:s}: {:d}".format(category,len(messages[category])))

def println():
    print()
        
def printdetails(messages,printcategories):
    # print a summary
    printed = False
    for category in printcategories:
        print("the following messages were encountered in category "+category+":")
        for message in messages[category]:
            printed = True
            print(category)
            printdict(message,1)
    return printed
            
def main(args):
    messages = {}
    with open(args.logfile,"rt") as logfile:
        parse(messages,logfile)

    printsummary(messages)

    println()
    
    if printdetails(messages,args.printcategories):
        println()
    
    # diagnose messages
    from RooFitUtils.logfile_diagnostics import diagnose
    print("detailed diagnosis")
    diagnose(messages)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="the log file doctor")
    parser.add_argument("logfile")
    parser.add_argument("--print",dest="printcategories",nargs="+",type=str,default=[])
    parser.add_argument("--debug",dest="debugcategory",type=str,default=None)        
    args = parser.parse_args()
    main(args)