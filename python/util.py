
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

def distributePointsAroundLine(dimlabels,coords,contour,npoints):
    """distribute points randomly around a line"""
    # closed contours have the first and the last point equal
    closed = equal(contour[0],contour[-1])
    lengths = lineSegmentLengths(contour)
    distpar = 0.005*sum(lengths)
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
        else:
            prevector = normalized(orthogonal(subtract(contour[end],contour[start])))
        if after != None:
            next_segment_reverse = normalized(subtract(contour[after],contour[end]))
            this_segment_forward = normalized(subtract(contour[start],contour[end]))
            postvector  = normalized(add(this_segment_forward,next_segment_reverse))
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

def makePOIstring(infilename):
  from ROOT import TFile, RooStats
  f = TFile.Open(infilename)
  for key in f.GetListOfKeys():
    if "ProcessID" not in key.GetName():
      wkspc = f.Get(key.GetName())
      modelconfig = wkspc.obj("ModelConfig")
      pois = modelconfig.GetParametersOfInterest()
      itr = pois.createIterator()
      pois, var = "", itr.Next()
      while var:
        pois = pois + var.GetName() + ","
        var = itr.Next()
      pois = pois[0:len(pois)-1]
  return pois

def retriveObj(filename):
  import ROOT
  file0 = ROOT.TFile.Open(filename)
  for x in file0.GetListOfKeys():
    if "ProcessID" not in x.GetName():
      obj = file0.Get(x.GetName())
  return obj

def makepctstring(pois,pct):
  pctstring, npois = "", len(pois.split(","))
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

def getnofNPs(hesse,pois):
  N = hesse.GetNcols()
  npois = len(pois.strip(","))
  return N - npois
