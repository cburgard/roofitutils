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

def flipped(l,flip=True):
    if flip:
        return reversed(l)
    else:
        return [i for i in l]

def flattened(l):
    """flatten a list of lists to a single list"""    
    for i in l:
        for j in i:
            yield j

def python_version():
    import sys
    return sys.version_info[0]

def println():
    if python_version() <= 2:
        print
    else:
        print()

def isdict(d):
    return isinstance(d,dict)

def islist(l):
    return isinstance(l,list)

def islist(l):
    return isinstance(l,tuple)
            
def isclose(a,b,rel_tol=1e-9,abs_tol=1e-9):
    return abs(a-b) <= max( rel_tol * max(abs(a), abs(b)), abs_tol )
            
def keys(d):
    if type(d) == dict:
        return d.keys()
    if type(d) == str:
        return [d]

def strip(s):
    return s.strip()

def striplatex(s):
    return s.replace("\n","").replace("\\toprule","").replace("\\bottomrule","").replace("\\hline","").replace("\\textbf","").replace("\\textit","").replace("{","").replace("}","").replace("$","").strip()
    
def concat(pieces,delim=" "):
    """concatenate a list of strings to a single string"""
    if isstr(pieces):
        return pieces
    else:
        return delim.join([p for p in pieces if p != None and p != ""])

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
    gSystem.Load("libRooFitUtils")

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
    if hasattr(coll,"createIterator"):
        itr = coll.createIterator()
        var = itr.Next()
        retval = []
        while var :
            retval.append(var)
            var = itr.Next()
        return retval
    else:
        return coll


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

def graphrescale(d,nllmin,scale=2):
    """transform the values of a curve by sorting them in x and subtracting a fixed value in y"""
    xvals = sorted(d.keys())
    yvals = [ max(d[k] - nllmin,0)*scale for k in xvals ]
    return list(zip(xvals,yvals))


def parsedict(s,cast=str):
    """parse a string of the format "a=b,c=d" into a dictionary"""
    d = {}
    for kv in s.split(","):
        try:
            k,v = kv.split("=")
        except ValueError:
            print("unable to parse key-value-pair '"+kv+"', need syntax 'key=value'")
            exit(0)
        d[k] = cast(v)
    return d

def product(pools):
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

def names(objects):
    return [ obj.GetName() for obj in makelist(objects) ]
        
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

def allkeys(d):
    elems = list(d.keys())
    for v in d.values():
        if type(v) == dict:
            elems += allkeys(v)
    return sorted(list(set(elems)))

def reconstructCall(args,arglist,blacklist):
    options = {}
    from argparse import _StoreTrueAction,_StoreFalseAction,_StoreAction
    for arg in args.items():
        actions = []
        for action in arglist:
            if arg[0] in blacklist: continue
            if arg[0] == action.dest:
                argval = arg[1]
                argopt = action.option_strings[0]
                if argval == True and isinstance(action,_StoreTrueAction):
                    options[argopt]=None
                if argval == False and isinstance(action,_StoreFalseAction):
                    options[argopt]=None
                if isinstance(action,_StoreAction) and not argval == action.default:
                    options[argopt]=argval
    return options

def lineSegmentLengths(points):
    segments = []
    from math import sqrt
    lastpoint = None
    for point in points:
        if not lastpoint:
            lastpoint = point
            continue
        else:
            l = length(subtract(point,lastpoint))
            segments.append(l)
            lastpoint = point
    return segments

def weighted_choice(weights):
    if len(weights) == 0:
        raise RuntimeError("unable to choose from empty list")
    import random
    totals = []
    running_total = 0

    for w in weights:
        running_total += w
        totals.append(running_total)

    rnd = random.random() * running_total
    for i, total in enumerate(totals):
        if rnd < total:
            return i
    raise RuntimeError("unable to choose with value "+str(rnd))

def interpolate(a,b,x):
    """interpolate between to vectors with a given interpolation in the range 0-1"""
    return add(scale(a,1-x),scale(b,x))

def length(a):
    """return the length of a vector"""
    from math import sqrt
    return sqrt(norm(a))

def norm(a):
    """return the norm of a vector"""
    return sum([ x*x for x in a ])

def subtract(a,b):
    """return the difference between to vectors"""
    return [ b[i] - a[i] for i in range(0,len(a)) ]

def add(a,b):
    """return the sum of two vectors"""
    return [ b[i] + a[i] for i in range(0,len(a)) ]

def normalized(a):
    """return a normalized copy of a vector"""
    if length(a) == 0:
        return scale(a,1.)
    else:
        return scale(a,1./length(a))

def scale(a,x):
    """scale a vector by a number"""
    return [ e*x for e in a ]

def equal(a,b):
    """test if two vectors are equal"""
    for i in range(0,len(a)):
        if a[i] != b[i]: return False
    return True

def orthogonal(a):
    """return a vector orthogonal to the one given"""
    from math import atan2,sin,cos,pi
    angle = atan2(a[1],a[0])
    # by deliberately swapping sin and cos, we introduce a shift by 90 deg.
    return [ sin(angle),cos(angle) ]

def distributePointsAroundLine(dimlabels,coords,contour,npoints,spread=0.005):
    """distribute points randomly around a line"""
    # closed contours have the first and the last point equal
    closed = equal(contour[0],contour[-1])
    lengths = lineSegmentLengths(contour)
    distpar = spread*sum(lengths)
    for p in range(0,npoints):
        # pick a segment randomly by length
        segment = weighted_choice(lengths)
        # for the vector calculations, we need the vectors of the previous and next segment
        # for closed-contours, this is almost trivial, except for the cases handled in the following if-then-else
        # for non-closed contours, we simply construct an orthogonal vector later on (and mark the segments as 'None' here)
        start = segment
        end   = segment+1
        if segment == 0:
            # we picked the first segment
            # if this is a closed contour, we use -2 instead of -1 in order to avoid a length of 0
            if closed: before = -2
            else:      before = None
        else:
            before = start-1
        if segment == len(lengths)-1:
            # we picked the last segment
            # if this is a closed contour, we use 1 instead of 0 in order to avoid a length of 0
            if closed:  after = 1
            else:       after = None
        else:
            after = end+1
        # randomly select the position along the segment
        import random
        pos = random.uniform(0, 1)
        linepoint = interpolate(contour[start],contour[end],pos)
        # calculate the trapezoid vectors at the beginning and end of the segment
        # by adding the vectors of the previous and current (current and next) segments, we get the "direction" in which the corner is pointing
        # by interpolating between the directions before and after the current segment, we cover the full area around the segment and don't skip-over/double count areas with edges
        if before != None:
            prev_segment_forward = normalized(subtract(contour[before],contour[start]))
            this_segment_reverse = normalized(subtract(contour[end],contour[start]))
            prevector   = normalized(add(prev_segment_forward,this_segment_reverse))
            if length(prevector) == 0:
                prevector = normalized(orthogonal(subtract(contour[end],contour[start])))                
        else:
            prevector = normalized(orthogonal(subtract(contour[end],contour[start])))
        if after != None:
            next_segment_reverse = normalized(subtract(contour[after],contour[end]))
            this_segment_forward = normalized(subtract(contour[start],contour[end]))
            postvector  = normalized(add(this_segment_forward,next_segment_reverse))
            if length(postvector) == 0:
                postvector = normalized(orthogonal(subtract(contour[end],contour[start])))
        else:
            postvector = normalized(orthogonal(subtract(contour[end],contour[start])))
        # calculate the shift by interpolating between the pre- and post-vector and scaling it randomly
        shift = scale(normalized(interpolate(prevector,postvector,pos)),distpar*random.gauss(0,1))
        # put everything together
        coord = add(linepoint,shift)
        coords.append( {dimlabels[i]:coord[i] for i in range(0,len(dimlabels)) } )

def distributePointsAroundPoint(dimlabels,coords,minimum,npoints,distpar):
    import numpy as np
    xy = np.random.multivariate_normal(minimum, distpar * np.identity(len(dimlabels)), npoints)
    for j in range(0,npoints):
        coords.append({dimlabels[i]:xy[j][i] for i in range(0,len(dimlabels))})

def stringify(s):
    if s == None: return ""
    if isinstance(s,list): return " ".join(['"'+v+'"' for v in s])
    return str(s)

def makepoint(coord):
    return ",".join([k+"="+str(v) for k,v in coord.items()])

def clearfile(filename):
    import os
    if os.path.isfile(filename):
        os.remove(filename)

def getThreshold(pcent,ndim):
    from scipy.stats import chi2
    return chi2.ppf(pcent,ndim)

def getPercent(nsigma):
    from scipy.stats import norm
    return 1 - (1-norm.cdf(nsigma))*2

def make1dscan(scan):
    return {k[0]:v for k,v in scan.items()}

def mergescan1d(scanvals1, scanvals2):
    import numpy as np

    # prepare input scan points for 1d interpolation
    # sort the scan points in x 
    x1,x2 = sorted(scanvals1.keys()), sorted(scanvals2.keys()) 
    y1 = [ scanvals1[val] for val in x1]  
    y2 = [ scanvals2[val] for val in x2]  
    # don't extrapolate,
    xmin = max(min(x1),min(x2))
    xmax = min(max(x1),max(x2))
    # grid granularity is finer than the inputs
    nx = 10*len(x1+x2)

    x_inter = np.linspace(xmin,xmax,nx)
    mergevals = {}
    
    y1_inter = np.interp(x_inter,x1,y1)
    y2_inter = np.interp(x_inter,x2,y2)

    # find minima of the two curves in the 1d-grid
    for i in range(len(x_inter)):
        mergevals[(x_inter[i],)] = float(min(y1_inter[i],y2_inter[i]))

    return mergevals

def mergescans1d(scan1,scan2):

    allvals1 = [ make1dscan(scan) for pnamelist,curves in scan1.items() for options,scan in curves.items()]
    allvals2 = [ make1dscan(scan) for pnamelist,curves in scan2.items() for options,scan in curves.items()]
    mergevals = mergescan1d(allvals1[0], allvals2[0])
    mergescans = {}
    for param in scan1:
        mergescans[param] = {}
        for option in scan1[param]:
            mergescans[param][option] = mergevals
    return mergescans
  
def createAsimov(ws,mc,asmName):
    import ROOT
    allParams = ROOT.RooArgSet()
    allParams.add(mc.GetGlobalObservables())
    allParams.add(mc.GetObservables())
    allParams.add(mc.GetNuisanceParameters())
    allParams.add(mc.GetParametersOfInterest())
    globs = mc.GetGlobalObservables().snapshot()
    asimovData = ROOT.RooStats.AsymptoticCalculator.MakeAsimovData(mc,allParams,globs)
    asimovData.SetName(asmName)
    getattr(ws,"import")(asimovData)

def makepctstring(pois,pct):
    pctstring, npois = "", len(pois)
    for i in range(0,npois): pctstring = pctstring + str(pct)+","
    return pctstring[0:len(pctstring)-1]

def getjobdims(dim,time):
    import math
    dimred = dim/1000.
    t = 25*dimred*dimred*dimred*dimred
    # split into n jobs such that 
    njobs = int(math.ceil(t/time))
    dimsize = int(math.ceil(dim/njobs))
    k = 0
    dims = []
    for i in range(0,njobs):
        lo = k 
        hi = k+dimsize
        if lo >=dim: break
        if hi >=dim: hi = dim
        dims.append([lo,hi])
        k = k+dimsize+1
    return dims

def countNPs(hesse,pois):
    N = hesse.GetNcols()
    npois = len(pois)
    return N - npois

def sgnstr(x):
    if x<0: return str(x)
    return "+"+str(x)

def parsegroups(pois):
    if isinstance(pois,dict):
        return None
    if len(pois) > 0 and isinstance(pois[0],tuple):
        return [ (t[0],len(t[-1])) for t in pois]
    return None


def parsepois(pois,options=False):
    if isdict(pois):
        poinames = list(pois.keys())
        poiopts = pois
    else:
        if isstr(pois[0]):
            poinames = pois
            poiopts = {}
        elif istuple(pois[0]):
            poinames = list(flattened([t[-1] for t in pois]))
            poiopts = {}
        else:
            poinames = [ p[0] for p in pois ]
            poiopts = { p[0]:p[1] for p in pois }
    if options:
        return poinames,poiopts
    else:
        return poinames

def first_n_nonzero_digits(l, n):
    from math import isnan
    if isnan(l): return 0
    from itertools import islice
    return int(''.join(islice((i for i in str(abs(l)) if i not in {'0', '.'}), n)).ljust(n, '0'))

def formatNumber(x,ndigits):
    from math import log10,isnan
    if isnan(x): return "NaN"
    if log10(abs(x)) < -ndigits: return "0"
    return "{:f}".format(round(x,ndigits)).rstrip("0").rstrip(".")

def findSignificantDigits(x):
    if x == 0: return 0
    from math import floor,log10,isnan
    if isnan(x): return "NaN"
    digits = first_n_nonzero_digits(x,3)
    scale = floor(log10(abs(x)))
    if digits < 354:
        return max(0,-scale+1)
    elif digits < 949:
        return max(0,-scale)
    else:
        return max(0,-scale-1)

def roundAutoUp(m):
    from math import log10
    unit = int(log10(m))
    return pow(10,unit) * (int(m / pow(10,unit))+1)
    
def formatPDG(x,xup,xdn,lap=False,show="v+-"):
    from math import isnan
    if isnan(x): return "NaN"
    elif isnan(xup) and isnan(xdn): return str(x)
    elif isnan(xup): digits = findSignificantDigits(abs(xdn))
    elif isnan(xdn): digits = findSignificantDigits(abs(xup))
    else:
        digits = max(findSignificantDigits(abs(xup)),findSignificantDigits(abs(xdn)))
    fmt = "{:."+str(int(digits))+"f}"
    retval = "\\ensuremath{"
    if "v" in show:
        num = fmt.format(x)
        if lap and "." in num:
            parts = num.split(".")
            num = "\\llap{\\ensuremath{"+parts[0]+"}}."+parts[1]
        retval += num
    else:
        retval += "{}"
    if "+" in show:
        retval += "^{+" + fmt.format(abs(xup)) +"}"
    if "-" in show:
        retval += "_{" + fmt.format(-abs(xdn))+"}"
    retval += "}"
    return retval
    
    
def formatNumberPDG(x,forceSign=False):
    digits = findSignificantDigits(x)
    fmt = "{:."+str(int(digits))+"f}"
    s = fmt.format(x)
    if not forceSign or x < 0:
        return s
    else:
        return "+"+s

def redirect_c_output(filename,keep_python_output=True):
    import sys
    import os
    sys.stdout.flush() # <--- important when redirecting to files

    # Duplicate stdout (file descriptor 1)
    # to a different file descriptor number
    newstdout = os.dup(1)
    newstderr = os.dup(2)    

    # open the output
    output = os.open(filename, os.O_WRONLY|os.O_CREAT)

    # Duplicate the file descriptor for output
    # and overwrite the value for stdout (file descriptor 1)
    os.dup2(output, 1)
    os.dup2(output, 2)    

    # Close output after duplication (no longer needed)
    os.close(output)

    # Use the original stdout to still be able
    # to print to stdout within python
    if keep_python_output:
        sys.stdout = os.fdopen(newstdout, 'w')
        sys.stderr = os.fdopen(newstderr, 'w')    

    return (newstdout,newstderr)

def restore_c_output(tup):
    import os
    import sys
    os.dup2(tup[0], 1)
    os.dup2(tup[1], 2)
    newstdout = os.dup(1)
    newstderr = os.dup(2)
    sys.stdout = os.fdopen(newstdout, 'w')
    sys.stderr = os.fdopen(newstderr, 'w')        

def extend(d,otherd={}):
    d.update(otherd)
    return d

def superset(thelist):
    return sorted(list(set(flattened(thelist))))

def typecast(d):
    if not isdict(d):
        raise TypeError("typecast only works on dictionaries!")
    retval = {}
    for k,v in d.items():
        if isdict(v):
            retval[k] = typecast(v)
            continue
        if islist(v):
            retval[k] = [typecast(e) for e in v]
            continue
        try:
            retval[k] = int(v)
            continue
        except ValueError:
            pass
        try:
            retval[k] = float(v)
            continue
        except ValueError:
            pass
        retval[k] = str(v).strip()
    return retval
        
    

def filterdict(data, types=[], keys=[], allow_empty=False):
    if isdict(data):
        result = {}  # dict-like, use dict as a base
        for k, v in data.items():
            if k in keys or any([isinstance(v, t) for t in types]):  # skip key/type
                continue
            try:
                result[k] = filterdict(v,types=types,keys=keys,allow_empty=allow_empty)
            except ValueError:
                pass
        if result or allow_empty:
            return result
    elif islist(data):
        result = []  # a sequence, use list as a base
        for v in data:
            if any([isinstance(v, t) for t in types]):  # skip type
                continue
            try:
                result.append(filterdict(v,types=types,keys=keys,allow_empty=allow_empty))
            except ValueError:
                pass
        if result or allow_empty:
            return result
    else:  # we don't know how to traverse this structure...
        return data  # return it as-is, hope for the best...
    raise ValueError

def parsecsv(infile,delimiter=',',quotechar='|'):
    """ requires a csv form where the first column serves as a label """
    import csv
    with open(infile) as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter, quotechar=quotechar)
        header = next(reader, None)
        tabledict = {}
        for row in reader:
            tabledict[row[0]] = {}
            for icol in range(1,len(header)):
                tabledict[row[0]][header[icol]] = float(row[icol])

    return tabledict

def dicttoarr(indict,names=False):
    """ convert a dict of dict of a matrix to an array
        pass names as True if the names of the parameters
        is to be returned """
    matarr = []
    for ipar1 in sorted(indict.keys()):
        tmparr = []
        for ipar2 in sorted(indict.keys()):
            tmparr.append(indict[ipar1][ipar2])
        matarr.append(tmparr)
    if names:
        return sorted(indict.keys()), matarr
    else: return matarr

def makecovmat(corarr,uncer):
    """ make array of a covariance matrix given a
        array of correlation matrix and uncertainities"""
    import math
    covarr, npar = [], len(uncer)
    for ipar1 in range(0,npar):
        val1 = uncer[ipar1]
        tmparr = []
        for ipar2 in range(0,npar):
            val2 = uncer[ipar2]
            tmparr.append(corarr[ipar1][ipar2]*val1*val2)
        covarr.append(tmparr)
    return covarr

def arrtotmat(arr):
    """ return array to a TMatrix object """
    import ROOT

    npar = len(arr)
    tmat = ROOT.TMatrixDSym(npar)
    for ipar1 in range(0,npar):
        for ipar2 in range(0,npar):
            tmat[ipar1][ipar2] = arr[ipar1][ipar2]
    return tmat

def getunc(indict,name="asm up"):
    import math
    uncer = []
    for par in sorted(indict.keys()):
        uncer.append(0.5*(indict[par][name+"up"]+abs(indict[par][name+"down"])))
    return uncer

