#!/usr/bin/env python

def sortparameters(parset,rankings):
    sorting = []
    for par in parset:
        score = 0.
        count = 0
        for ranking in rankings.values():
            score += ranking[par][0]*ranking[par][0] + ranking[par][1]*ranking[par][1]
            count += 1
        sorting.append((par,score))
    sortedlist = sorted(sorting,key=lambda x:x[1])
    sortedparams = [ e[0] for e in sortedlist ]
    for par in sorted(parset):
        if not par in sortedparams:
            sortedparams.append(par)
    return list(reversed(sortedparams))

def filterparameters(parlist,blacklist):
    newlist = []
    for par in parlist:
        remove = False
        for pat in blacklist:
            import re
            regex = re.compile(pat)
            if regex.match(par):
                remove = True
        if not remove:
            newlist.append(par)
    return newlist

def specifyparameters(parlist,whitelist):
    newlist = []
    for parin in whitelist:
        exist = False
        newlist.append(parin)
        for par in parlist:
            if parin == par:
                exist = True
        if not exist:
            print('ERROR ',parin,'cannot be found in the workspace')
    return newlist

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser("plot pulls of parameters")
    parser.add_argument('--range',nargs=2,default=[-5,5],help="range to plot",type=float)    
    parser.add_argument('--scaleimpacts',default=1,help="Scale impacts for better visualization",type=float)
    parser.add_argument('-i','--input',action='append',nargs="+",metavar=('drawoptions','file.txt'),help="text files with the input information",default=[])
    parser.add_argument('--impacts',action='append',nargs="+",metavar=('POI','file.txt'),help="text files with the input information for ranking the nuisance parameters",default=[])
    parser.add_argument('--breakdown',action='append',nargs="+",metavar=('POI','file.txt'),help="text files with the input information for ranking the nuisance parameters",default=[])
    parser.add_argument('--blacklist',nargs="+",metavar="NP",help="parameters and groups to avoid showing",default=[])    
    parser.add_argument('--whitelist',nargs="+",metavar="NPw",help="parameters to show",default=[])
    parser.add_argument("--atlas",type=str,help="ATLAS plot label, will enable ATLAS style if used",required=False,default=None)
    parser.add_argument("--numbers",action="store_true",help="show numbers",default=False)    
    parser.add_argument("--labels",help="position of the labels",choices=["l","r"],default="r")    
    
    parser.add_argument('-o',"--output",type=str,help="output file name",default="pulls.tex",required=True)
    args = parser.parse_args()

    scans = {}
    pullresults = {}
    from RooFitUtils.io import collectresults
    for inset in args.input:
        label = inset[0]
        files = inset[1:]
        collectresults(scans,pullresults,files,label)
        
    parset = set(pullresults.keys())
    rankings = {}
    from RooFitUtils.io import collectimpacts
    for inset in args.impacts:
        poiname = inset[0]
        files = inset[1:]
        parset = parset.union(set(collectimpacts(rankings,files,poiname)))
    from RooFitUtils.io import collectbreakdowns
    for inset in args.breakdown:
        poiname = inset[0]
        files = inset[1:]
        parset = parset.union(set(collectbreakdowns(rankings,files,poiname)))        

    # whitelist used only if the blacklist is not specified
    filtallpars = []
    if len(args.blacklist) < 1 and len(args.whitelist) > 0:
        filtallpars = specifyparameters(parset,args.whitelist)
    else:
        filtallpars = filterparameters(parset,args.blacklist)
    allpars = sortparameters(filtallpars,rankings)
    npar = len(allpars)
        
    with open(args.output,"wt") as outfile:
        from RooFitUtils.pgfplotter import writehead,writefoot,writepullshead,writepullsfoot,writepull,writepullnumbers,writeranking,writerankinghead
        writehead(outfile,args.atlas)
        writepullshead(args,outfile,allpars,args.labels)
        writerankinghead(args,outfile,allpars,args.scaleimpacts)
        for np in range(0,npar):
            for poi,ranking in rankings.items():
                if allpars[np] in ranking.keys():
                    npname = allpars[np]
                    writeranking(args,rankings[poiname][npname],outfile,np-npar,args.scaleimpacts)
            if allpars[np] in pullresults.keys():
                writepull(args,pullresults,outfile,allpars[np],np-npar,True)
                writepullnumbers(args,pullresults,outfile,allpars[np],np-npar,args.labels,args.numbers)
        writepullsfoot(args,outfile,allpars)        
        writefoot(outfile)
