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
                    elements[label].append(parsed)
                    break
            if self.is_parsed(parsed):
                continue
            if not self.parent:
                line = next(lines).strip()
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
        line = lines.peek()
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
    parser = MetaParser({
        "minimum":Parser(r"Minuit2Minimizer : (?P<State>\w+) minimum - status = (?P<Status>\d+)",MetaParser({
            "NLL":      Parser(r"FVAL\s+= (?P<VAL>"+NUM+")"),
            "Edm":      Parser(r"Edm\s+= (?P<VAL>"+NUM+")"),
            "Nfcn":     Parser(r"Nfcn\s+= (?P<VAL>"+NUM+")"),
            "Timing":             Parser(r"Real time (?P<RealTime>"+TIME+"), CP time (?P<CPTime>"+NUM+")"),
            "Timing+Slices":      Parser(r"Real time (?P<RealTime>"+TIME+"), CP time (?P<CPTime>"+NUM+"), (?P<Slices>\d+) slices"),            
            "FloatParameter":Parser(r"(?P<name>\w+)\s+=\s+(?P<val>"+NUM+")\s+\+/-\s+(?P<err>"+NUM+")"),
            "ConstParameter":Parser(r"(?P<name>\w+)\s+=\s+(?P<val>"+NUM+")\s+\(fixed\)"),            
        })),
        "minimization":Parser(r"Minuit2Minimizer: Minimize with max-calls (?P<MaxCalls>\d+) convergence for edm < (?P<Edm>"+NUM+") strategy (?P<Strategy>\d)",MetaParser({
            "DefaultOptionChange":Parser(r"Minuit2Minimizer::Minuit\s+- Changing default options",MetaParser({
                "StorageLevel":       Parser(r"\s+StorageLevel\s+:\s+(?P<StorageLevel>\d+)"),
            })),
            "InitialState":       Parser(r"Info in <Minuit2>: MnSeedGenerator Initial state:"),# FCN =\s+(?P<FCN>"+NUM+")"),#\s+Edm =\s+(?P<Edm>"+NUM+")\s+NCalls =\s+(?P<NCalls>"+NUM+")"),
            "StartIteration":     Parser(r"Info in <Minuit2>: VariableMetricBuilder Start iterating until Edm is < (?P<Eps>"+NUM+") with call limit = (?P<CallLimit>\d+)"),
            "Iteration":          Parser(r"Info in <Minuit2>: VariableMetricBuilder\s+(?P<Iteration>\d+)\s+-\s+FCN\s+=\s+(?P<FCN>"+NUM+")\s+Edm\s+=\s+(?P<Edm>"+NUM+")\s+NCalls\s+=\s+(?P<NCalls>\d+)"),
            "AfterHessian":       Parser(r"Info in <Minuit2>: VariableMetricBuilder After Hessian"),
            "FlatLH":             Parser(r"Minuit2:0: RuntimeWarning: VariableMetricBuilder No improvement in line search"),
            "NegativeG2":         Parser(r"Info in <Minuit2>: MnSeedGenerator Negative G2 found - new state:",MetaParser({
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
    for i in range(0,len(minimizations)):
        minimization = minimizations[i]
        if "NegativeG2" in minimization.keys():
            nNegG2 += 1
        else:
            nOk += 1
    if nOk == i+1:
        print("all {:d} minimizations OK".format(nOk))
    else:
        print("some minimizations encountered issues ({:d}/{:d}):".format(nOk,i+1))
        if nNegG2 > 0:
            print("  {:d} instances of negative G2".format(nNegG2))
    
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

    
        
