def make_parser():
    # this is the gory, catch-all function to create the parsing tree for all known lines in a logfile
    # this function will never be complete and needs to be extended continuously until all possible messages in a logfile can be identified
    # this will of course never happen, but the completeness will improve over time

    from RooFitUtils.logfile_diagnostics import Parser,MetaParser
    from RooFitUtils.util import extend
    
    NUM=r'[+-]?\d+[.]?\d*[e]?[+-]?\d*'
    TIME=r"\d+:\d+:\d+"

    minimum = {
            "NLL":      Parser(r"FVAL\s+= (?P<VAL>"+NUM+")"),
            "Edm":      Parser(r"Edm\s+= (?P<VAL>"+NUM+")"),
            "Nfcn":     Parser(r"Nfcn\s+= (?P<VAL>"+NUM+")"),
            "Timing":             Parser(r"Real time (?P<RealTime>"+TIME+"), CP time (?P<CPTime>"+NUM+")"),
            "Timing+Slices":      Parser(r"Real time (?P<RealTime>"+TIME+"), CP time (?P<CPTime>"+NUM+"), (?P<Slices>\d+) slices"),            
            "FloatParameter":Parser(r"(?P<name>\w+)\s+=\s+(?P<val>"+NUM+")\s+\+/-\s+(?P<err>"+NUM+")"),
            "FloatParameterNoErr":Parser(r"(?P<name>\w+)\s+=\s+(?P<val>"+NUM+")"),        
            "ConstParameter":Parser(r"(?P<name>\w+)\s+=\s+(?P<val>"+NUM+")\s+\(fixed\)"),            
        }


    hessestate = {
        "MinVal":Parser(r"\s*Minimum value\s*:\s*(?P<val>"+NUM+")"),
        "Edm":Parser(r"\s*Edm\s*:\s*(?P<val>"+NUM+")"),
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
    }
    
    minimization = {
        "DefaultOptionChange":Parser(r"Minuit2Minimizer::Minuit\s+-\s+Changing default.*options",MetaParser({
            "StorageLevel":       Parser(r"StorageLevel\s+:\s+(?P<StorageLevel>\d+)"),
        })),
        "InitialState1":       Parser(r"MnSeedGenerator: for initial parameters FCN\s*=\s*(?P<initFCN>"+NUM+")"),
        "InitialState1a":      Parser(r"(?P<source>.*): Initial state[:]?\s*-\s*FCN\s*=\s*(?P<FCN>"+NUM+")\s*Edm\s*=\s*(?P<Edm>"+NUM+")\s*NCalls\s*=\s*(?P<NCalls>"+NUM+")"),
        "InitialState2":       Parser(r"Info in <Minuit2>: MnSeedGenerator Initial state:"),# FCN =\s+(?P<FCN>"+NUM+")"),#\s+Edm =\s+(?P<Edm>"+NUM+")\s+NCalls =\s+(?P<NCalls>"+NUM+")"),
        "StartIteration":     Parser(r".*VariableMetric.*\s*[Ss]tart iterating until Edm is < (?P<Eps>"+NUM+")[ with call limit = (?P<CallLimit>\d+)]?"),
        "lowTolerance":       Parser(r".*VariableMetric.*\s*Tolerance is not sufficient, continue the minimization",MetaParser({
            "edminfo":            Parser(r"Info in\s*(?P<Label>\w+)\s*Edm\s*is\s*:\s*edm[val]*\s*=\s*(?P<edm>"+NUM+")")
        })),
        "lowTolerance2":  Parser(r".*VariableMetric.*Tolerance not sufficient, continue minimization; Edm (?P<edm>"+NUM+") Required (?P<required>"+NUM+")"),
        "iterationsFinishNoConvergence": Parser(r".*VariableMetric.*\s*[Ii]terations finish without convergence.*",MetaParser({
            "edminfo":            Parser(r"Info in\s*(?P<Label>\w+)\s*: \s*edm[val]*\s*=\s*(?P<edm>"+NUM+")"),
            "edmreq":            Parser(r"Info in\s*requested\s*: edmval\s*=\s*(?P<edmreq>"+NUM+")"),
            "invalid":            Parser(r"Info in\s*(?P<Label>\w+)\s*: INVALID function minimum - edm is above tolerance.*:\s*edm[val]*\s*=\s*(?P<edm>"+NUM+")"),
            "reqtol":            Parser(r"Info in\s*(?P<Label>\w+)\s*: Required tolerance\s*is 10 x edmval\s*:\s*edm[val]*\s*=\s*(?P<edm>"+NUM+")"),            
        })),
        "MinimumInvalidTry2": Parser(r".*VariableMetricBuilder FunctionMinimum is invalid after second try"),
        "DavidonErrorUpdator": Parser(".*DavidonErrorUpdator[:]?\s*delgam\s*<\s*0\s*:\s*first derivatives increasing along search line.*"),
        "NoImprovementInLineSearch": Parser(r".*VariableMetricBuilder: no improvement in line search"),
        "Iteration":          Parser(r".*VariableMetric.*\s*(?P<it>\d+)\s+-\s+FCN\s+=\s+(?P<FCN>"+NUM+")\s+Edm\s+=\s+(?P<Edm>"+NUM+")\s+NCalls\s+=\s+(?P<NCalls>\d+)"),
        "AfterHessian":       Parser(r".*VariableMetric.*\s*After Hessian"),
        "NotConverged" :      Parser(r".*Minimization did NOT converge, Edm is above max"),
        "FlatLH":             Parser(r"Minuit2:0: RuntimeWarning: VariableMetricBuilder No improvement in line search"),
        "runHesseVerbose":           Parser(r".*MnSeedGenerator[:]? run Hesse\s*-\s*new state:\s*-\s*FCN\s*=\s*(?P<FCN>"+NUM+")\s*Edm\s*=\s*(?P<Edm>"+NUM+")\s*NCalls\s*=\s*(?P<NCalls>"+NUM+")"),
        "runHesse":           Parser(r".*MnSeedGenerator[:]? run Hesse\s*-\s*new state:",MetaParser(hessestate)),
        "hesseCalls":         Parser(r".*Hesse [Uu]sing max-calls (?P<maxcalls>\d+)"),
        "hesseInfo":          Parser(r".*::Hesse[ :]*Hesse is (?P<hesseStatus>\w+) - matrix is (?P<matrixStatus>\w+)",MetaParser(minimum)),
        "hesseNonPosDef": Parser(r".*VariableMetricBuilder: matrix not pos.def. : edm is < 0. Make pos def...",MetaParser({        
            "posdef":Parser(r"Info in matrix forced pos-def by adding to diagonal : padd = (?P<padd>"+NUM+")")
        })),
        "hesseNonPosDefGDel": Parser(r".*VariableMetricBuilder: matrix not pos.def, gdel > 0",MetaParser({        
            "gdel":Parser(r"Info: gdel = (?P<gdel>"+NUM+")"),
            "posdef":Parser(r"Info in matrix forced pos-def by adding to diagonal : padd = (?P<padd>"+NUM+")")            
        })),
        "NegativeG2":         Parser(r".*MnSeedGenerator[:]? Negative G2 found - new state:",MetaParser(hessestate))}

    minos_migrad = extend({
        "parameter":Parser(r"Pos (?P<pos>\d+): (?P<parname>.*) = (?P<val>"+NUM+")"),
        "result":Parser(r".*MnFunctionCross Result after [()23rndMigrad: ]* FCN =\s+(?P<FCN>"+NUM+")\s+Edm =\s+(?P<Edm>"+NUM+")\s+NCalls =\s+(?P<NCalls>\d+)",MetaParser({
            "header":Parser(r"Pos\s+\|\s+Name\s+\|\s+type\s+\|\s+Value\s+\|\s+Error\s+[+]/[-]"),
            "parameter":Parser(r"(?P<Pos>\d+)\s+\|\s+(?P<Name>.*)\s+\|\s+(?P<type>\w+)\s+\|\s+(?P<Value>"+NUM+")\s+\|\s+(?P<Error>"+NUM+")")
        })),
    },minimization)
    
    minos = extend({
        "mnminos":Parser(r".*MnMinos Determination of (?P<direction>\w+) Minos error for parameter (?P<parno>\d+)"),
        "mncross":       Parser(r".*MnFunctionCross: parameter 0 set to (?P<val>"+NUM+")",MetaParser(minimization)),
        "end":           Parser(r".*MnMinos end of Minos scan for (?P<direction>\w+) interval for parameter (?P<parname>.*)"),
        "error":         Parser(r".*MnMinos could not find (?P<direction>\w+) Value for Parameter\s*:\s*par_name = (?P<parname>.*)"),
        "invalid":       Parser(r".*Minos:\s*Invalid (?P<direction>\w+) error for parameter (?P<parname>.*)"),
        "matinvfail":    Parser(r".*MnGlobalCorrelationCoeff: inversion of matrix fails."),
        "newmin":        Parser(r".*Minos[:]?\s*[Nn]ew Minimum found while looking for Parameter\s*: par_name = (?P<parname>.*)"),
        "newmin2":       Parser(r".*Minos[:]?\s*[Nn]ew Minimum found while running Minos for (?P<direction>\w+) error for parameter (?P<parname>.*)"),
        "rerun":         Parser(r".*Minuit2Minimizer::GetMinosError : Found a new minimum: run again the Minimization\s*starting from the new\s*point",MetaParser({
            "newmin":      Parser(r"New minimum point found by Minos\s*:",MetaParser(minimum)),
        })),
        "limit":       Parser(r".*Minos\s*Parameter is at (?P<direction>\w+) limit. : par_name = (?P<parname>.*)"),
        "limit2":       Parser(r".*Minos:\s*Parameter\s*:\s*(?P<parname>.*)\s*is at (?P<direction>\w+) limit"),                        
        "mncross-migrad":Parser(r".*MnFunctionCross[: ]+Run Migrad [again (23rnd)]*with fixed parameters:",MetaParser(minos_migrad)),
        "result":Parser(r"Minos: (?P<direction>\w+) error for parameter (?P<parname>.+)\s+:\s+(?P<val>"+NUM+")"),
        "Timing":             Parser(r"Real time (?P<RealTime>"+TIME+"), CP time (?P<CPTime>"+NUM+")"),
        "Timing+Slices":      Parser(r"Real time (?P<RealTime>"+TIME+"), CP time (?P<CPTime>"+NUM+"), (?P<Slices>\d+) slices"),
        "stars":Parser("[*]+"),        
    },minimization)
    
    parser = MetaParser({
        "intro":Parser(r"RooFit v(?P<Version>[\d.]+) -- Developed by Wouter Verkerke and David Kirkby",MetaParser({
            "copyright":Parser(r"Copyright \(C\) (?P<Years>[\d-]+) NIKHEF, University of California & Stanford University"),
            "rights":Parser(r"All rights reserved, please read http://roofit.sourceforge.net/license.txt")
        })),
        "minimum":Parser(r".*Minuit2Minimizer\s*:\s*(?P<State>\w+) [Mm]inimum - status = (?P<Status>\d+)",MetaParser(minimum)),
        "roofitresult":Parser(r"RooFitResult: minimized FCN value: (?P<FCN>"+NUM+"), estimated distance to minimum: (?P<edm>"+NUM+")",MetaParser({
            "covqual":Parser(r"covariance matrix quality: (?P<status>.*)"),
            "status":Parser(r"Status : (?P<tags>.*)"),
            "floatpars":Parser(r"Floating Parameter\s+Final\s*Value [+]/[-]\s*Error",MetaParser({
                "dashes":Parser(r"\-+\s+\-+"),
                "floatPar":Parser(r"(?P<parname>.*)\s+(?P<val>"+NUM+")\s*[+]/[-]\s*(?P<err>"+NUM+")")
            }))
        })),
        "minos":      Parser(r".*GetMinosError for parameter (?P<parno>\d+) (?P<parname>.+) using max-calls (?P<maxcalls>\d+), tolerance (?P<tolerance>\d+)",MetaParser(minos)),
        "minos2":     Parser(r".*GetMinosError - Run MINOS (?P<direction>\w+) error for parameter #(?P<parno>\d+) : (?P<parname>.+) using max-calls (?P<maxcalls>\d+), tolerance (?P<tolerance>\d+)",MetaParser(minos)),
        "minos2bidir": Parser(r".*GetMinosError - Run MINOS (?P<direction1>\w+) and (?P<direction2>\w+) error for parameter #(?P<parno>\d+) : (?P<parname>.+) using max-calls (?P<maxcalls>\d+), tolerance (?P<tolerance>\d+)",MetaParser(minos)),
        "minosrerun":Parser(r".*Minuit2Minimizer::GetMinosError : Run now again Minos from the new found Minimum",MetaParser(minos)),
        "minimization":Parser(r"Minuit2Minimizer: Minimize with max-calls (?P<MaxCalls>\d+) convergence for edm < (?P<Edm>"+NUM+") strategy (?P<Strategy>\d)",MetaParser(minimization)),
        "quickfit-snapshot":Parser(r"REGTEST: Loading snapshot (?P<Snapshot>.*)"),        
        "quickfit-preparing":Parser(r"Preparing parameters of interest\s*:\s*(?P<poiset>.*)",MetaParser({
            "firstpoi":Parser(r"REGTEST: Set first POI to (?P<POI>.*)"),
            "summary":Parser(r"Summary of POI setup in the likelihood model",MetaParser({
                "poi":Parser(r".*RooRealVar::\s*(?P<poiname>.*)\s*=\s*(?P<poival>.*)\s*L\((?P<poirange>.*)\)")
            }))
        })),
        "quickfit-startminos":Parser(r"Evaluating MINOS errors for all POIs...",MetaParser(minimization)),        
        "quickfit-starthesse":Parser(r"Starting fit with HESSE...",MetaParser(minimization)),
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
        "sframework-initialize-parameter":Parser(r".*SFramework/(?P<appname>.*): initializing parameter '(?P<name>.*)' to '(?P<value>"+NUM+")'"),
        "sframework-skip-constant-parameter":Parser(r".*WARNING: skipping constant parameter '(?P<name>.*)'"),
        "sframework-float-parameters":Parser(r".*SFramework/(?P<appname>.*): parameters are: (?P<info>.*)"),                
        "sframework-initialize-parameter-error":Parser(r".*unable to initialize parameter '(?P<name>.*)'.*"),                
        "sframework-loadsnapshot":Parser(r".*SFramework/(?P<appname>.*): loaded snapshot '(?P<snapshot>.*)' for Nll construction, (?P<npar>\d+) parameters floating"),
        "sframework-createnll":Parser(r".*SFramework/(?P<appname>.*): constructed Nll with (?P<npar>\d+) floating parameters and starting value (?P<nllval>"+NUM+")"),        
        "robustminimizer-start":Parser(r".*ExtendedMinimizer::robustMinimize\((?P<appname>.)*\): starting minimization with strategy (?P<Strategy>\d+)"),
        "robustminimizer-starthesse":Parser(r".*ExtendedMinimizer::runHesse\((?P<appname>.*)\) running after minimization \(this might take a while\) ...",MetaParser(minimization)),
        "robustminimizer-startminos":Parser(r".*ExtendedMinimizer::minimize\((?P<appname>.*)\): Running Minos",MetaParser({
            "stars":Parser("[*]+"),
        })),
        "robustminimizer-covqual":Parser(r".*ExtendedMinimizer::minimize\((?P<appname>.*)\): Covariance quality is (?P<covqual>"+NUM+")"),
        "robustminimizer-invert":Parser(r".*ExtendedMinimizer::minimize\((?P<appname>.*)\): attempting to invert covariance matrix..."),        
        "robustminimizer-endhesse":Parser(r".*ExtendedMinimizer::runHesse\((?P<appname>.*)\) finished with status (?P<status>\d) \(covqual=(?P<covqual>\d)\)"),
        "robustminimizer-end"  :Parser(r".*ExtendedMinimizer::robustMinimize\((?P<appname>.*)\) fit succeeded with status (?P<Status>\d+)"),
        "robustminimizer-retry"  :Parser(r".*ExtendedMinimizer::robustMinimize\((?P<appname>.*)\) fit failed with status (?P<Status>[-]?\d+). [Rr]etrying with strategy (?P<Strategy>\d+)",MetaParser(minimization)),
        "robustminimizer-giveup"  :Parser(r".*ExtendedMinimizer::robustMinimize\((?P<appname>.*)\) fit failed with status (?P<Status>[-]?\d+)[, giving up.]?"),
        "roofitutils-close":Parser(r"Fitting time: (?P<Time>"+NUM+")s",MetaParser({
            "nll":Parser(r"NLL after minimisation: (?P<Time>"+NUM+")"),
            "poi":Parser(r"(?P<poiname>.*) = (?P<val>"+NUM+") (?P<uperr>"+NUM+") (?P<dnerr>"+NUM+")"),
            "noop":Parser(r"no (?P<obj>\w+) requested")
        }))
    })
    return parser

