#!/usr/bin/env python
import re
import sys

global ignoreComponents
global namemap
ignoreComponents = []
namemap = []
blacklist = []

_files = []

tolerance = pow(10,-4)


colors = {
    "blue":'\033[94m',
    "green":'\033[92m',
    "yellow":'\033[93m',
    "red":'\033[91m',
    "bold":'\033[1m',
    "underline":'\033[4m',
}

def colored(s,c):
    return colors[c] + s + '\033[0m'


def isclose(a, b, rel_tol):
    return abs(a-b) <= rel_tol * max(abs(a), abs(b))

def notin(a,l):
    for x in l:
        if x.match(a):
            return False
    return True

def groupnames(m):
    return tuple( m.groupdict().keys() )

def namelist(args):
    retval = []
    if hasattr(args,"createIterator"):
      iter = args.createIterator()
    else:
      iter = args
    try:
        var = iter.Next()
        while var :
            retval.append(var.GetName())
            var = iter.Next()
    except AttributeError:
        for var in args: retval.append(var.GetName())
    return retval

def list_and(inlist):
    retval = True
    for x in inlist:
        retval = x and retval
    return retval

def list_or(inlist):
    retval = False
    for x in inlist:
        retval = x or retval
    return retval

def makelist(args,type=None):
    retval = []

    if hasattr(args,"createIterator"):
        iter = args.createIterator()
    else:
        iter = args
    try:
        var = iter.Next()
        while var :
            if type and not var.InheritsFrom(type.Class()): continue
            retval.append(var)
            var = iter.Next()
    except AttributeError:
        for var in args:
            if type and not var.InheritsFrom(type.Class()): continue
            retval.append(var)
    return retval

def filter(inlist,whitelist,blacklist=[]):
    return [ x for x in inlist if (list_or([x.InheritsFrom(y.Class()) for y in whitelist]) and not list_or([x.InheritsFrom(y.Class()) for y in blacklist])) ]

def removeAll(s,l):
    mys = str(s)
    for x in l:
        mys = x.sub("",mys)
    return mys


def binValues(dhist):
    values = []
    import ROOT
    for obs in makelist(dhist.get()):
        if obs.InheritsFrom(ROOT.RooRealVar.Class()):
            for i in range(0,obs.getBins()):
                x = dhist.get(i)
                values.append(dhist.weight())
        if obs.InheritsFrom(ROOT.RooAbsCategory.Class()):
            continue
    return values

def compare(first,second):
    #returns true if either one of the two names of the TObjects provided are equal
    #if one of the two arguments is None False is returned, if both are None True is returned
    global ignoreComponents
    f = removeAll(first,ignoreComponents)
    s = removeAll(second,ignoreComponents)
    if f == s: return True
    if f.startswith(s) and f[len(s)]=="_": return True
    if s.startswith(f) and s[len(f)]=="_": return True
    return False

def compareNames(first,second):
    if (not first) and (not second):
        return True
    if (not first) or (not second):
        return False
    return compare(first.GetName(),second.GetName())


def findByName(somelist,name):
    for obj in somelist:
        ma,mb = None,None
        for a,b in namemap:
            ma = a.match(obj.GetName())
            mb = b.match(name)
            if ma and mb: break
        if ma and mb:
            names = groupnames(ma)
            ok = True
            for x in names:
                if ma.group(x) != mb.group(x):
                    ok = False
            if ok:
                return obj
        if compare(obj.GetName(),name):
            return obj
    return None

def diffObjects(objtype,label,firstvars,secondvars):
    missing_in_second = []
    for firstvar in firstvars:
        secondvar = findByName(secondvars,firstvar.GetName())
        if not secondvar and notin(firstvar.GetName(),blacklist):
            missing_in_second.append(firstvar.GetName())

    if len(missing_in_second)>0:
        print("the following "+objtype+" are missing in "+label+":")
        print(",".join(sorted(missing_in_second)))

def diffVarVariations(firstlabel,firstvars,firstpdf,firstnset,secondlabel,secondvars,secondpdf,secondnset):
    print("name "+firstlabel+" "+secondlabel)
    pairs = []
    for firstvar in firstvars:
        secondvar = findByName(secondvars,firstvar.GetName())
        if secondvar:
            if firstvar.isConstant() and secondvar.isConstant(): continue
            pairs.append((firstvar,secondvar))
    for a,b in pairs:
        aval = a.getVal()
        bval = b.getVal()
        na = firstpdf.expectedEvents(firstnset)
        nb = secondpdf.expectedEvents(secondnset)
        a.setVal(aval+abs(a.getErrorHi()))
        b.setVal(bval+abs(b.getErrorHi()))
        nahi = firstpdf.expectedEvents(firstnset)
        nbhi = secondpdf.expectedEvents(secondnset)
        a.setVal(aval-abs(a.getErrorLo()))
        b.setVal(bval-abs(b.getErrorLo()))
        nalo = firstpdf.expectedEvents(firstnset)
        nblo = secondpdf.expectedEvents(secondnset)
        a.setVal(aval)
        b.setVal(bval)
        if not isclose(na,nb,tolerance) or not isclose(nahi,nbhi,tolerance) or not isclose(nalo,nblo,tolerance):
            print("{:s}/{:s}:{:s}/{:s} {:f} {:f} {:f} : {:f} {:f} {:f}".format(firstpdf.GetName(),secondpdf.GetName(),a.GetName(),b.GetName(),na,nahi-na,nalo-na,nb,nbhi-nb,nblo-nb))


def getHistFunc(datahist,ws):
    import ROOT
    obs = datahist.get().first()
    wsobs = ws.var(obs.GetName())
    if not wsobs or not wsobs.InheritsFrom(ROOT.RooRealVar.Class()): return None
    for obj in makelist(wsobs.clientIterator()):
        if obj.InheritsFrom(ROOT.RooHistFunc.Class()):
            if obj.dataHist() == datahist: return obj
    return None

def getInterpolationVariable(hf):
    if not hf: return None
    for c in makelist(hf.clientIterator()):
        if c.InheritsFrom(ROOT.PiecewiseInterpolation.Class()):
            for high,low,var in zip(makelist(c.highList()),makelist(c.lowList()),makelist(c.paramList())):
                if hf == high or hf == low:
                    return var
    return None

def longestCommonSubstring(s1, s2):
    m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(s1)):
        for y in range(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest: x_longest]


def getSisterFunction(hf,varname):
    useHigh = False
    useLow = False
    highFunc = None
    lowFunc = None
    for c in makelist(hf.clientIterator()):
        if c.InheritsFrom(ROOT.PiecewiseInterpolation.Class()):
            for high,low,var in zip(makelist(c.highList()),makelist(c.lowList()),makelist(c.paramList())):
                if hf == high: useHigh=True
                elif hf == low: useLow=True
                if compare(var.GetName(),varname):
                    highFunc = high
                    lowFunc = low
    if useHigh: return highFunc
    elif useLow: return lowFunc
    return None



def diffFuncVals(firstlabel,firstfuncs,secondlabel,secondfuncs):
    print("name "+firstlabel+" "+secondlabel)
    for firstfunc in firstfuncs:
        secondfunc = findByName(secondfuncs,firstfunc.GetName())
        if secondfunc and not firstfunc.getVal() == secondfunc.getVal():
            print("{:s} {:f} {:f}".format(firstfunc.GetName(),firstfunc.getVal(),secondfunc.getVal()))

def diffDataVals(firstlabel,firstdata,secondlabel,seconddata,firstws,secondws):
    print("name "+firstlabel+" "+secondlabel)
    for firstd in firstdata:
        firsthf = getHistFunc(firstd,firstws)
        firstvar = getInterpolationVariable(firsthf)

        secondd = findByName(seconddata,firstd.GetName())
        if not secondd: continue

        secondhf = getHistFunc(secondd,secondws)
        secondvar = getInterpolationVariable(secondhf)
        if firstvar and secondvar and not compareNames(firstvar,secondvar):
            secondhf = getSisterFunction(secondhf,firstvar.GetName())
            secondd = secondhf.dataHist()

        firstvals = binValues(firstd)
        secondvals = binValues(secondd)

        equal = True
        for a,b in zip(firstvals,secondvals):
            if not a == b:
                equal = False
        if not equal:

            if firstvar:
                print(firstd.GetName()+"/"+secondd.GetName() + " ("+firstvar.GetName()+"/"+secondvar.GetName()+") " + " ".join(["{:f}/{:f}".format(a,b) for a,b in zip(firstvals,secondvals)]))
            else:
                print(firstd.GetName() + " " + " ".join(["{:f}/{:f}".format(a,b) for a,b in zip(firstvals,secondvals)]))


def getWS(filename):
    import ROOT
    rootfile = ROOT.TFile.Open(filename,"READ")
    _files.append(rootfile)
    for k in rootfile.GetListOfKeys():
        if k.GetClassName() == "RooWorkspace":
            return rootfile.Get(k.GetName())

def sortobjects(mylist):
    retval = {}
    for obj in mylist:
        if obj.ClassName() not in retval.keys(): retval[obj.ClassName()] = []
        retval[obj.ClassName()].append(obj)
        for k in retval.keys():
            retval[k] = sorted(retval[k],key = lambda x: str(x.GetName()).lower()+"_")
    return retval

def checkCompleteness(orig,sync):
    for o in orig:
        if not o: continue
        if not o in sync:
            raise RuntimeError("object '{:s}' is not in output list".format(o.GetName()))

def synchronizeLists(first,second, recurse=True):
    #returns tuple of two lists corresponding to 'first' and 'second' with additionally inserted None entries in order to align the two lists based on the names of the contained objets (using compareNames)
    #do not externally provide the 'recurse' argument!
    firstOut = []
    secondOut = []
    offset = 0
    for index in range(0,len(first)):
        f = first[index]
        localOffset = 0
        s = None
        while ( (index+offset+localOffset) < len(second)):
            s = second[index+offset+localOffset]
            if compareNames(f,s):
                #we got a match stop iterating
                break
            s = None
            localOffset += 1 #for next try

        if (index+offset+localOffset) < len(second): #no matching element in second list
            #we got a match
            firstOut += [None]*localOffset #we skipped this many elements from the second list until the match. Hence, fill first list with the same number of dummy (None) entries
            secondOut += second[index+offset:index+offset+localOffset] #add the corresponding (skipped) elements from the second input to the second output
            firstOut.append(f) #add the matched elements
            secondOut.append(s)
            offset += localOffset
        else: #no match, add f to firstOut and a dummy to secondOut
            if f:
                firstOut.append(f)
                secondOut.append(None)
                offset -= 1 #original first list advances by one more step compared to the second one.
            elif (index+offset<len(second)) and second[index+offset]: #we might have 'None' in first but a non-None in the second list. We therefore need to add the one in the second input to the second output (and None to the first output)
                firstOut.append(None)
                secondOut.append(second[index+offset])
    #we might not have reached the end of 'second' yet, so add the remaining items from second input to second output (and Nones to first output)
    if (len(first)+offset < len(second)):
        secondOut += second[len(first)+offset-len(second) : ]
        firstOut += [None]*(len(second)-len(first)-offset) #flipped sign w.r.t. previous line -> ok!
    #sanity check for the logic above, uncomment for debugging
    #checkCompleteness(first,firstOut)
    #checkCompleteness(second,secondOut)
    return (firstOut,secondOut)

def diffPdfTrees(firstpdf,secondpdf,indent=0):
    myindent = indent+1
    firstcomps  = makelist(firstpdf.serverIterator())
    secondcomps = makelist(secondpdf.serverIterator())
    firstfuncs  = sortobjects(filter(firstcomps, [ROOT.RooAbsReal],[ROOT.RooConstVar]))
    secondfuncs = sortobjects(filter(secondcomps,[ROOT.RooAbsReal],[ROOT.RooConstVar]))

    classes = set(firstfuncs.keys())
    classes.update(secondfuncs.keys())
    firstfuncsSync,secondfuncsSync = ({},{})
    for c in classes:
        if (c in firstfuncs) and (c in secondfuncs):
            firstfuncsSync[c],secondfuncsSync[c] = synchronizeLists(firstfuncs[c],secondfuncs[c])
        elif c in secondfuncs:
            secondfuncsSync[c] = secondfuncs[c]
            firstfuncsSync[c] = [None]*len(secondfuncs[c])
        elif c in firstfuncs:
            firstfuncsSync[c] = firstfuncs[c]
            secondfuncsSync[c] = [None]*len(firstfuncs[c])

    for c in classes:
        if not c in firstfuncsSync.keys():
            for elem in secondfuncsSync[c]:
                print("  "*(myindent)+colored(elem.GetName()+"/(None)","red")+ " ("+elem.ClassName()+")")
        elif not c in secondfuncsSync.keys():
            for elem in firstfuncsSync[c]:
                print("  "*(myindent)+colored("(None)/"+elem.GetName(),"red")+ " ("+elem.ClassName()+")")
        else:
            for i in range(max(len(firstfuncsSync[c]),len(secondfuncsSync[c]))):
                if i < len(firstfuncsSync[c]): first = firstfuncsSync[c][i]
                else: first = None
                if i < len(secondfuncsSync[c]): second = secondfuncsSync[c][i]
                else: second = None
                if first and second:
                    #if not first.GetName().startswith(second.GetName()) and not second.GetName().startswith(first.GetName()):
                    if not compareNames(first,second):
                        print("  "*(myindent) + colored("{:s}/{:s} ".format(first.GetName(),second.GetName()),"red") + " ("+c+")")
                        diffPdfTrees(first,second,myindent)
                    else:
                        print("  "*(myindent) + first.GetName() + " (" + c + ")")
                        diffPdfTrees(first,second,myindent)
                elif first:
                    print("  "*(myindent) + colored("{:s}/(None) ".format(first.GetName()),"red") + " ("+c+")")
                elif second:
                    print("  "*(myindent) + colored("(None)/{:s} ".format(second.GetName()),"red") + " ("+c+")")




if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("wsdiff",description="investigate differences between workspaces")
    parser.add_argument("first",type=str,help="first workspace")
    parser.add_argument("second",type=str,help="second workspace")
    parser.add_argument("--ignore",nargs="+",type=str,default=["auto","alpha","Constraint$","constraint$","_"],help="ignore certain name components")
    parser.add_argument("--blacklist",nargs="+",type=str,default=["nom_.*"],help="do not show certain objects")
    parser.add_argument("--map",nargs="+",type=str,help="map certain names",default=[])
    parser.add_argument('--models',                           dest='models', action='store_true'   , help="show models", default=True )
    parser.add_argument('--no-models',                        dest='models', action='store_false'  , help="do not show models", default=True )
    parser.add_argument('--variables',                           dest='vars', action='store_true'   , help="show variables", default=True )
    parser.add_argument('--no-variables',                        dest='vars', action='store_false'  , help="do not show variables", default=True )
    parser.add_argument('--variations',                           dest='variations', action='store_true'   , help="show variations", default=True )
    parser.add_argument('--no-variations',                        dest='variations', action='store_false'  , help="do not show variations", default=True )
    parser.add_argument('--functions',                           dest='funcs', action='store_true'   , help="show functions", default=True )
    parser.add_argument('--no-functions',                        dest='funcs', action='store_false'  , help="do not show functions", default=True )
    parser.add_argument('--data',                           dest='data', action='store_true'   , help="show data", default=True )
    parser.add_argument('--no-data',                        dest='data', action='store_false'  , help="do not show data", default=True )
    parser.add_argument('--all',                           dest='all', action='store_true'   , help="show list of all objects", default=True )
    parser.add_argument('--no-all',                        dest='all', action='store_false'  , help="do not show list of all objects", default=True )
    parser.add_argument('--tolerance',                      dest='tolerance', type=int, default=4  , help="number of digits in precision to compare" )
    args = parser.parse_args()
    first = getWS(args.first)
    first.SetTitle(args.first)
    tolerance = pow(10,-args.tolerance)

    second = getWS(args.second)
    second.SetTitle(args.second)

    import re
    ignoreComponents = [ re.compile(x) for x in args.ignore ]
    blacklist = [ re.compile(x) for x in args.blacklist ]

    for mapfile in args.map:
        with open(mapfile,"r") as infile:
            for line in infile:
                if len(line.strip()) == 0: continue
                if line.startswith("#"):
                    continue
                bits = line.split()
                if not len(bits) == 2:
                    print("unable to read line :"+line)
                    continue
                a = re.compile(bits[0])
                b = re.compile(bits[1])
                namemap.append((a,b))
                namemap.append((b,a))

    import ROOT
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

    firstmodels = makelist(first.allGenericObjects(),ROOT.RooStats.ModelConfig)
    secondmodels = makelist(second.allGenericObjects(),ROOT.RooStats.ModelConfig)

    if args.models:
        print("== models ==")
        diffObjects("Models",first.GetTitle(),secondmodels,firstmodels)
        print("")
        diffObjects("Models",second.GetTitle(),firstmodels,secondmodels)
        print("")
        for firstmc,secondmc in zip(firstmodels,secondmodels):
            diffPdfTrees(firstmc.GetPdf(),secondmc.GetPdf())

            print("POIs:")
            firstpois = makelist(firstmc.GetParametersOfInterest())
            secondpois = makelist(secondmc.GetParametersOfInterest())
            diffObjects("POIs",first.GetTitle(),firstpois,secondpois)
            diffObjects("POIs",second.GetTitle(),secondpois,firstpois)
            print("NPs:")
            firstnps = makelist(firstmc.GetNuisanceParameters())
            secondnps = makelist(secondmc.GetNuisanceParameters())
            diffObjects("NPs",first.GetTitle(),firstnps,secondnps)
            diffObjects("NPs",second.GetTitle(),secondnps,firstnps)
            print("Obs.:")
            firstobs = makelist(firstmc.GetObservables())
            secondobs = makelist(secondmc.GetObservables())
            for obs in firstobs:
                if obs.InheritsFrom(ROOT.RooRealVar.Class()):
                    obs.setVal(obs.getMin())
                    obs.setConstant(True)
                if obs.InheritsFrom(ROOT.RooAbsCategory.Class()):
                    obs.setIndex(0)
            for obs in secondobs:
                if obs.InheritsFrom(ROOT.RooRealVar.Class()):
                    obs.setVal(obs.getMin())
                    obs.setConstant(True)
                if obs.InheritsFrom(ROOT.RooAbsCategory.Class()):
                    obs.setIndex(0)
            diffObjects("Obs",first.GetTitle(),firstobs,secondobs)
            diffObjects("Obs",second.GetTitle(),secondobs,firstobs)

    if args.vars:
        print("\n\n== variables ==")
        firstvars = makelist(first.allVars())
        secondvars = makelist(second.allVars())
        diffObjects("Variables",first.GetTitle(),secondvars,firstvars)
        print("")
        diffObjects("Variables",second.GetTitle(),firstvars,secondvars)
        print("")
        if args.variations:
            for firstmc,secondmc in zip(firstmodels,secondmodels):
                diffVarVariations(first.GetTitle(),firstvars,firstmc.GetPdf(),firstmc.GetObservables(),second.GetTitle(),secondvars,secondmc.GetPdf(),secondmc.GetObservables())

    if args.funcs:
        print("\n\n== functions ==")
        firstfuncs = makelist(first.allFunctions())
        secondfuncs = makelist(second.allFunctions())
        diffObjects("Functions",second.GetTitle(),firstfuncs,secondfuncs)
        print("")
        diffObjects("Functions",first.GetTitle(),secondfuncs,firstfuncs)
        print("")
        diffFuncVals(first.GetTitle(),firstfuncs,second.GetTitle(),secondfuncs)

    if args.data:
        print("\n\n== datasets ==")
        firstdata = first.allData()
        seconddata = second.allData()
        diffObjects("Datasets",second.GetTitle(),firstdata,seconddata)
        print("")
        diffObjects("Datasets",first.GetTitle(),seconddata,firstdata)
        print("")
        diffDataVals(first.GetTitle(),firstdata,second.GetTitle(),seconddata,first,second)

    if args.all:
        print("\n\n== generic objects ==")
        firstfuncs = makelist(first.allGenericObjects())
        secondfuncs = makelist(second.allGenericObjects())
        diffObjects("Objects",second.GetTitle(),firstfuncs,secondfuncs)
        print("")
        diffObjects("Objects",first.GetTitle(),secondfuncs,firstfuncs)

