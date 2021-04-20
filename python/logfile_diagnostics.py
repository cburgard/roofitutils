class MetaParser:
    # a MetaParser is a class intended to parse a block of logfile lines
    # a MetaParser owns a colleciton of subparsers
    # for every line, the MetaParser attempts to pass the line to every one of its parsers
    # if a line is encountered that does not fit any of the parsers this MetaParser is managing, one of two things will happen
    #   if the MetaParser has a parent parser, it will return there, leaving the line in-place
    #   if the MetaParser does not have a parent parser, it will flag the line as "unknown" and move on
    
    def __init__(self,parsers):
        self.parsers = parsers
        for k,v in parsers.items():
            v.label = k
        self.debug = False
        self.parent = None
        
    def set_debug(self,debug):
        self.debug = debug        
        for parser in self.parsers.values():
            parser.set_debug(debug)

    def is_parsed(self,elem):
        return elem != False
            
    def parse(self,elements,lines):
        while lines:
            parsed = False
            for label,parser in self.parsers.items():
                parsed = parser.parse(lines)
                if self.is_parsed(parsed):
                    if not label in elements.keys():
                        elements[label] = []
                    parsed[".item"] = label
                    elements[label].append(parsed)
                    break
            if self.is_parsed(parsed):
                continue
            if not self.parent:
                line = next(lines).strip()
                if not "unknown" in elements.keys():
                    elements["unknown"] = []
                elements["unknown"].append(line)
            else:
                line = lines.peek().strip("\n")
                if self.debug:
                    print(self.parent.label + " encountered line: '"+line+"' unknown to any of: " + " ".join(self.parsers.keys()))
                return False
        return True

class Parser:
    # a Parser is a class intended to parse a single line of logfile code
    # there is exactly one "trigger" regex that causes this parser to recognize and claim ownership of a line
    # this line is then read, processed, and the information retruned
    # if the line is not recongized, the parser will just return False
    #
    # a Parser can (but does not need to) have an associated MetaParser
    # that is activated right after a line has been claimed successfully
    
    def __init__(self,trigger_regex,metaparser=None):
        import re
        self.trigger = re.compile(trigger_regex)
        self.subparser = metaparser
        if self.subparser:
            self.subparser.parent = self
        self.debug = False
        self.label = None

    def set_debug(self,debug):
        self.debug = debug
        if self.subparser:
            self.subparser.set_debug(debug)
        
    def parse(self,lines):
        line = lines.peek().strip()
        while not line:
            next(lines)
            line = lines.peek().strip()            
        match = self.trigger.match(line)
        if not match:
            if self.debug:
                print(self.label + " unable to parse line: '"+line.strip("\n")+"'")
            return False
        if self.debug:
            print(self.label + " parsed line '"+line.strip()+"'")
        keys = match.groupdict()
        next(lines)
        if self.subparser:
            self.subparser.parse(keys,lines)
        return keys

def make_parser():
    # this is the gory, catch-all function to create the parsing tree for all known lines in a logfile
    # this function will never be complete and needs to be extended continuously until all possible messages in a logfile can be identified
    # this will of course never happen, but the completeness will improve over time
    
    NUM=r'[+-]?\d+[.]?\d*[e]?[+-]?\d*'
    TIME=r"\d+:\d+:\d+"

    minimum = MetaParser({
            "NLL":      Parser(r"FVAL\s+= (?P<VAL>"+NUM+")"),
            "Edm":      Parser(r"Edm\s+= (?P<VAL>"+NUM+")"),
            "Nfcn":     Parser(r"Nfcn\s+= (?P<VAL>"+NUM+")"),
            "Timing":             Parser(r"Real time (?P<RealTime>"+TIME+"), CP time (?P<CPTime>"+NUM+")"),
            "Timing+Slices":      Parser(r"Real time (?P<RealTime>"+TIME+"), CP time (?P<CPTime>"+NUM+"), (?P<Slices>\d+) slices"),            
            "FloatParameter":Parser(r"(?P<name>\w+)\s+=\s+(?P<val>"+NUM+")\s+\+/-\s+(?P<err>"+NUM+")"),
            "ConstParameter":Parser(r"(?P<name>\w+)\s+=\s+(?P<val>"+NUM+")\s+\(fixed\)"),            
        })

    
    minimization = MetaParser({
        "DefaultOptionChange":Parser(r"Minuit2Minimizer::Minuit\s+- Changing default options",MetaParser({
            "StorageLevel":       Parser(r"StorageLevel\s+:\s+(?P<StorageLevel>\d+)"),
        })),
        "InitialState1":       Parser(r"MnSeedGenerator: for initial parameters FCN\s*=\s*(?P<initFCN>"+NUM+")"),
        "InitialState1a":      Parser(r"(?P<source>.*): Initial state[:]?\s*-\s*FCN\s*=\s*(?P<FCN>"+NUM+")\s*Edm\s*=\s*(?P<Edm>"+NUM+")\s*NCalls\s*=\s*(?P<NCalls>"+NUM+")"),
        "InitialState2":       Parser(r"Info in <Minuit2>: MnSeedGenerator Initial state:"),# FCN =\s+(?P<FCN>"+NUM+")"),#\s+Edm =\s+(?P<Edm>"+NUM+")\s+NCalls =\s+(?P<NCalls>"+NUM+")"),
        "StartIteration":     Parser(r".*VariableMetric.*\s*[Ss]tart iterating until Edm is < (?P<Eps>"+NUM+")[ with call limit = (?P<CallLimit>\d+)]?"),
        "lowTolerance":       Parser(r".*VariableMetric.*\s*Tolerance is not sufficient, continue the minimization"),
        "edminfo":            Parser(r"Info in\s*(?P<Label>\w+)\s*Edm\s*is\s*:\s*edm[val]*\s*=\s*(?P<edm>"+NUM+")"),
        "Iteration":          Parser(r".*VariableMetric.*\s*(?P<it>\d+)\s+-\s+FCN\s+=\s+(?P<FCN>"+NUM+")\s+Edm\s+=\s+(?P<Edm>"+NUM+")\s+NCalls\s+=\s+(?P<NCalls>\d+)"),
        "AfterHessian":       Parser(r".*VariableMetric.*\s*After Hessian"),
        "FlatLH":             Parser(r"Minuit2:0: RuntimeWarning: VariableMetricBuilder No improvement in line search"),
        "runHesse":           Parser(r".*MnSeedGenerator[:]? run Hesse\s*-\s*new state:\s*-\s*FCN\s*=\s*(?P<FCN>"+NUM+")\s*Edm\s*=\s*(?P<Edm>"+NUM+")\s*NCalls\s*=\s*(?P<NCalls>"+NUM+")"),
        "hesseCalls":         Parser(r".*Hesse using max-calls (?P<maxcalls>\d+)"),
        "hesseInfo":          Parser(r".*::Hesse : Hesse is (?P<hesseStatus>\w+) - matrix is (?P<matrixStatus>\w+)",minimum),
        "NegativeG2":         Parser(r".*MnSeedGenerator[:]? Negative G2 found - new state:",MetaParser({
            "MinVal":Parser(r"\s+Minimum value\s+: (?P<val>"+NUM+")"),
            "Edm":Parser(r"\s+Edm\s+: (?P<val>"+NUM+")"),
            "InternalParameters":Parser(r"\s*Internal parameters\s*:",MetaParser({
                "value":Parser(r"\s*(?P<value>"+NUM+")")
            })),
            "InternalGradient":Parser(r"\s*Internal gradient\s*:",MetaParser({
                "value":Parser(r"\s*(?P<value>"+NUM+")")
            })),
            "InternalCovMat":Parser(r"\s*Internal covariance matrix\s*:",MetaParser({
                "value":Parser(r"\s*(?P<value>"+NUM+"\s*)+")
            })),
            "posdef":Parser(r"Info in matrix forced pos-def by adding to diagonal : padd = (?P<padd>"+NUM+")"),
            "hesseposdef":Parser(r"Info: MnHesse: matrix was forced pos. def.")
        }))})
        
    parser = MetaParser({
        "intro":Parser(r"RooFit v(?P<Version>[\d.]+) -- Developed by Wouter Verkerke and David Kirkby",MetaParser({
            "copyright":Parser(r"Copyright \(C\) (?P<Years>[\d-]+) NIKHEF, University of California & Stanford University"),
            "rights":Parser(r"All rights reserved, please read http://roofit.sourceforge.net/license.txt")
        })),
        "minimum":Parser(r"Minuit2Minimizer : (?P<State>\w+) minimum - status = (?P<Status>\d+)",minimum),
        "roofitresult":Parser(r"RooFitResult: minimized FCN value: (?P<FCN>"+NUM+"), estimated distance to minimum: (?P<edm>"+NUM+")",MetaParser({
            "covqual":Parser(r"covariance matrix quality: (?P<status>.*)"),
            "status":Parser(r"Status : (?P<tags>.*)"),
            "floatpars":Parser(r"Floating Parameter\s+Final\s*Value [+]/[-]\s*Error",MetaParser({
                "dashes":Parser(r"\-+\s+\-+"),
                "floatPar":Parser(r"(?P<parname>.*)\s+(?P<val>"+NUM+")\s*[+]/[-]\s*(?P<err>"+NUM+")")
            }))
        })),
        "minos":      Parser(r".*GetMinosError for parameter (?P<parno>\d+) (?P<parname>.+) using max-calls (?P<maxcalls>\d+), tolerance (?P<tolerance>\d+)",MetaParser({
            "mncross":Parser(r"MnFunctionCross: parameter 0 set to (?P<val>"+NUM+")",minimization),
            "result":Parser(r"Minos: (?P<direction>\w+) error for parameter (?P<parname>.+)\s+:\s+(?P<val>"+NUM+")"),
            "Timing":             Parser(r"Real time (?P<RealTime>"+TIME+"), CP time (?P<CPTime>"+NUM+")"),
            "Timing+Slices":      Parser(r"Real time (?P<RealTime>"+TIME+"), CP time (?P<CPTime>"+NUM+"), (?P<Slices>\d+) slices"),                        
        })),
        "minimization":Parser(r"Minuit2Minimizer: Minimize with max-calls (?P<MaxCalls>\d+) convergence for edm < (?P<Edm>"+NUM+") strategy (?P<Strategy>\d)",minimization),
        "quickfit-snapshot":Parser(r"REGTEST: Loading snapshot (?P<Snapshot>.*)"),        
        "quickfit-preparing":Parser(r"Preparing parameters of interest\s*:\s*(?P<poiset>.*)",MetaParser({
            "firstpoi":Parser(r"REGTEST: Set first POI to (?P<POI>.*)"),
            "summary":Parser(r"Summary of POI setup in the likelihood model",MetaParser({
                "poi":Parser(r".*RooRealVar::\s*(?P<poiname>.*)\s*=\s*(?P<poival>.*)\s*L\((?P<poirange>.*)\)")
            }))
        })),
        "quickfit-startminos":Parser(r"Evaluating MINOS errors for all POIs...",minimization),        
        "quickfit-starthesse":Parser(r"Starting fit with HESSE...",minimization),
        "quickfit-saveresults":Parser(r"Saving results to (?P<outpath>.*)"),
        "quickfit-done":Parser(r"All fits done in (?P<cputime>"+NUM+") min \(cpu\), (?P<realtime>"+NUM+") min \(real\)",MetaParser({
            "poisummary":Parser(r"Fit Summary of POIs \( STATUS (?P<status>\w+) \)",MetaParser({
                "dashes":Parser(r"\-+"),
                "FloatParameter":Parser(r"RooRealVar::(?P<name>.+)\s+=\s+(?P<val>"+NUM+")\s+\+/\-\s+(?P<err>.*)\s+L\((?P<range>.*)\)"),
                "ConstParameter":Parser(r"RooRealVar::(?P<name>.+)\s+=\s+(?P<val>"+NUM+")\s+\(fixed\)")
            }))
        })),
        "quickfit-start":Parser(r"Starting fit...",MetaParser({
            "quickfit-binned":Parser(r"set binned likelihood for:\s*(?P<region>.*)"),        
            "quickfit-buildnll":Parser(r"Building NLL...",MetaParser({
                "time":Parser("NLL built in (?P<cputime>.*) min \(cpu\), (?P<realtime>.*) min \(real\)")
            }))
        })),
        "robustminimizer-start":Parser(r".*ExtendedMinimizer::robustMinimize\(minimizer\): starting minimization with strategy (?P<Strategy>\d+)"),
        "robustminimizer-end"  :Parser(r".*ExtendedMinimizer::robustMinimize\(minimizer\) fit succeeded with status (?P<Status>\d+)"),
        "roofitutils-close":Parser(r"Fitting time: (?P<Time>"+NUM+")s",MetaParser({
            "nll":Parser(r"NLL after minimisation: (?P<Time>"+NUM+")"),
            "noop":Parser(r"no (?P<obj>\w+) requested")
        }))
    })
    return parser

def diagnose_minimizations(minimizations):
    # diagnose all the messages from the minimization branch for any severe issues
    i = 0
    nOk = 0
    nNegG2 = 0
    nLowTol = 0
    for i in range(0,len(minimizations)):
        minimization = minimizations[i]
        ok = 1
        if "NegativeG2" in minimization.keys():
            nNegG2 += 1
            ok = 0
        if "lowTolerance" in minimization.keys():
            nLowTol += 1
            ok = 0
        nOk += ok
    if nOk == i+1:
        print("all {:d} minimizations OK".format(nOk))
    else:
        print("some minimizations encountered issues ({:d}/{:d}):".format(nOk,i+1))
        if nNegG2 > 0:
            print("  {:d} instances of negative G2".format(nNegG2))
        if nLowTol > 0:
            print("  {:d} instances of insufficient tolerance".format(nLowTol))
    
def diagnose_minima(minima):
    # diagnose all the messages from the minima branch for any severe issues    
    i = 0
    nOk = 0
    for i in range(0,len(minima)):
        minimum = minima[i]
        if int(minimum["Status"]) != 0:
            print("minimum {:d} has status {:d} ".format(i,int(minimum["Status"])))
        else:
            nOk += 1
    if nOk == i+1:
        print("all {:d} minima OK".format(nOk))
    else:
        print("some minima were not OK ({:d}/{:d})".format(nOk,i+1))        
        

def diagnose(messages):
    # evaluate all known diagnostics
    diagnose_minima(messages["minimum"])
    diagnose_minimizations(messages["minimization"])    

    
        
