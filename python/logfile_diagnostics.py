class MetaParser:
    # a MetaParser is a class intended to parse a block of logfile lines
    # a MetaParser owns a colleciton of subparsers
    # for every line, the MetaParser attempts to pass the line to every one of its parsers
    # if a line is encountered that does not fit any of the parsers this MetaParser is managing, one of two things will happen
    #   if the MetaParser has a parent parser, it will return there, leaving the line in-place
    #   if the MetaParser does not have a parent parser, it will flag the line as "unknown" and move on
    
    def __init__(self,parsers):
        if not parsers:
            raise RuntimeError("empty MetaParser dictionary passed!")
        self.parsers = parsers
        for k,v in parsers.items():
            try:
                v.label = k
            except AttributeError:
                raise RuntimeError("object passed along as a parser '"+k+"' is not one, intead "+str(type(v)))
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
            try:
                self.subparser.parent = self
            except AttributeError:
                raise RuntimeError("object passed as metaparser with regex "+str(trigger_regex)+" is not one, instead '"+str(type(metaparser))+"'")
        self.debug = False
        self.label = None

    def set_debug(self,debug):
        self.debug = debug
        if self.subparser:
            self.subparser.set_debug(debug)
        
    def parse(self,lines):
        from RooFitUtils.util import typecast
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
        keys = typecast(match.groupdict())
        next(lines)
        if self.subparser:
            self.subparser.parse(keys,lines)
        return keys

class Severity:
    INFO = 0
    WARNING = 1
    ERROR = 5
    FATAL = 10
    
class IssueHandler:
    def __init__(self,key,severity,info):
        self.key = key
        self.severity = severity
        self.info = info

    def handle(self,instance):
        raise NotImplementedError()
        
    def handle_all(self,messages):
        summary = []
        if self.key in messages.keys():
            for instance in messages[self.key]:
                summary.append(self.handle(instance))
        return summary

class Handler(IssueHandler):
    def __init__(self,key,message="",severity=Severity.ERROR,info=""):
        super().__init__(key=key,severity=severity,info=info)
        self.message = message
    
    def handle(self,instance):
        try:
            return self.message.format_map(instance)
        except KeyError as e:
            raise RuntimeError("cannot fine key "+str(e)+" in "+str(instance))
    
def diagnose_minimizations(minimizations):
    # diagnose all the messages from the minimization branch for any severe issues
    i = 0
    nOk = 0

    issue_handlers = [
        Handler("NegativeG2",  severity = Severity.ERROR, info="some of your parameters might be interdependent",
                message = "the second derivative matrix was made invertible by adding {posdef[0][padd]} to the diagonal."),
        Handler("lowTolerance",severity = Severity.WARNING, info="the minimization needed to be prolonged as the minimum was not yet reached",
                message = "the estimated distances to minimum were recorded as {edminfo[0][edm]} ({edminfo[0][Label]}) and {edminfo[1][edm]} ({edminfo[1][Label]})."),        
    ]
    
    errors = {}
    warnings = {}    
    for i in range(0,len(minimizations)):
        minimization = minimizations[i]
        ok = True
        for handler in issue_handlers:
            summary = handler.handle_all(minimization)
            if summary:
                if handler.severity >= Severity.ERROR:
                    errors[handler.key] = {"info":handler.info,"messages":summary}
                    ok = False
                else:
                    warnings[handler.key] = {"info":handler.info,"messages":summary}
        if ok:
            nOk+=1
    if nOk == i+1:
        print("all {:d} minimizations OK".format(nOk))
    else:
        print("some minimizations encountered errors ({:d}/{:d}):".format(i+1-nOk,i+1))
        for key,summary in errors.items():
            print("  {:d} instances of {:s}. {:s}".format(len(summary["messages"]),key,summary["info"]))
            for elem in summary["messages"]:
                if elem:
                    print("    "+elem)
        if len(warnings) > 0:
            print("some minimizations encountered warnings:")
            for key,summary in warnings.items():
                print("  {:d} instances of {:s}. {:s}".format(len(summary["messages"]),key,summary["info"]))
                for elem in summary["messages"]:
                    if elem:
                        print("    "+elem)                    
    
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

    
        
