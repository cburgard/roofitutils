#!/bin/evn python

_nodel = []
def nodel(something):
    """mark an object as undeletable for the garbage collector"""
    _nodel.append(something)

try:
    isinstance("", basestring)
    def isstr(s):
        """test if an object is a string"""
        return isinstance(s, basestring)
except NameError:
    def isstr(s):
        """test if an object is a string"""
        return isinstance(s, str)

def concat(pieces,delim=" "):
    """concatenate a list of strings to a single string"""
    if isstr(pieces):
        return pieces
    else:
        return delim.join(pieces)

def shellexec(command,inputs=[],verbose=False,allowErrors=True):
    """execute a command by invoking it in a shell"""
    if verbose:
        print(concat(command))
        for inputline in inputs:
            if len(inputline) > 0: print("> "+concat(inputline))
            pass
        pass
    command = concat(command).split()
    stdinlines = [ concat(pieces) for pieces in inputs if len(pieces) > 0]
    s = concat(stdinlines,"\n")

    from subprocess import Popen, PIPE
    p = Popen(command, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    b = s.encode('ascii')
    stdout_data, stderr_data = p.communicate(input=b)
    if len(stderr_data) > 0:
        if allowErrors and verbose:
            QFramework.ERROR("Error from binary '"+concat(command)+"': "+concat(stderr_data))
        if not allowErrors:
            QFramework.BREAK("Error from binary '"+concat(command)+"': "+concat(stderr_data))
    return stdout_data

def mkdir(path):
    """call 'mkdir -p' on the given path"""
    shellexec(["mkdir","-p",path])

def mod(a,b):
    """calculate modulus"""
    return a%b

def linspace(vmin,vmax,npoints):
    """create a linear spacing of points in the given range"""
    step = (vmax - vmin)/(npoints)
    vals = [ float(vmin + i*step) for i in range(0,npoints+1) ]
    return vals

def loadRooFitUtils():
    """load the RooFitUtils shared library"""
    # retrieve the root core dir environment variable
    from ROOT import gSystem
    if gSystem.Load("libRooFitUtils"):
        raise ImportError("unable to load standalone libRooFitUtils.so!")

def timestamp(seconds):
    """get the current timestamp"""
    s = "{0:.3f}s".format(mod(seconds,60))
    if seconds < 60:
        return s
    m = "{0:d}m".format(int(mod(seconds/60,60)))
    if seconds < 3600:
        return m+" "+s
    h = "{0:d}h".format(int(seconds/3600))
    return h+" "+m+":"+s

def makelist(coll):
    """turn a RooFit collection into a list"""
    itr = coll.createIterator()
    var = itr.Next()
    retval = []
    while var :
        retval.append(var)
        var = itr.Next()
    return retval

def union(listoflists):
    """create the union of a list of lists"""
    s = set()
    for l in listoflists:
        s.update(l)
    return list(s)

def vec(l,t):
    """turn a python list into a c++ vector of the given type"""
    import ROOT
    v = ROOT.vector(t)()
    for e in l:
        v.push_back(e)
    return v

def graphoffset(d,nllmin):
    """transform the values of a curve by sorting them in x and subtracting a fixed value in y"""
    xvals = sorted(d.keys())
    yvals = [ max(d[k] - nllmin,0) for k in xvals ]
    return list(zip(xvals,yvals))
        

def parsedict(s):
    """parse a string of the format "a=b,c=d" into a dictionary"""
    d = {}
    for kv in s.split(","):
        k,v = kv.split("=")
        d[k] = v
    return d

def product(pools):
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

def generateCoordsDict(scan):
    """generate a dictionary of coordinates"""
    nvar = len(scan)
    names = []
    ranges = []
    for dim in scan:
        parname = str(dim[0])
        n = int(dim[1])
        minx = float(dim[2])
        maxx = float(dim[3])
        names.append(parname)
        ranges.append(linspace(minx,maxx,n))
    return [ {names[i]:point[i] for i in range(len(names))} for point in product(ranges) ]

def parsePoint(line):
    """parse a point into a dictionary"""
    return { n.strip():float(v) for (n,v) in [x.split("=") for x in line.split(",") ]}

def reconstructCall(args,arglist,blacklist):
    options = {}
    from argparse import _StoreTrueAction,_StoreFalseAction,_StoreAction
    for arg in args.__dict__:
        actions = []
        for action in arglist:
            if arg in blacklist: continue
            if arg == action.dest:
                argval = args.__dict__[arg]
                argopt = action.option_strings[0]
                if argval == True and isinstance(action,_StoreTrueAction):
                    options[argopt]=None
                if argval == False and isinstance(action,_StoreFalseAction):
                    options[argopt]=None
                if isinstance(action,_StoreAction) and not argval == action.default:
                    options[argopt]=argval
    return options

def stringify(s):
    if s == None: return ""
    if isinstance(s,list): return " ".join(['"'+v+'"' for v in s])
    return str(s)

def makepoint(coord):
    return ",".join([k+"="+str(v) for k,v in coord.items()])

def clearfile(filename):
    import os
    if os.path.exists(filename):
        os.remove(filename)

