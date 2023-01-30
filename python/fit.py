verbosity = {"DEBUG":4,"INFO":3,"WARNING":2,"ERROR":1}

def simplefit(args):
    """this is a simple drop-in replacement running a fit that is useful for debugging"""
    from ROOT import TFile
    infile = TFile.Open(args["inFileName"],"READ")
    if not infile or not infile.IsOpen(): raise RuntimeError("unable to open file "+args["inFileName"])    
    workspace = infile.Get(args["wsName"])
    if not workspace: raise RuntimeError("unable to find workspace "+args["wsName"])
    model = workspace.obj(args["modelConfigName"])
    if not model: raise RuntimeError("unable to obtain "+args["modelConfigName"])        
    pdf = model.GetPdf()
    if not pdf: raise RuntimeError("unable to obtain pdf from model!")
    data = workspace.data(args["dataName"])
    if not data: raise RuntimeError("no dataset '"+args["dataName"])
    pdf.fitTo(data)


def setup(args):
    # general setup and loading of modules
    import ROOT

    from RooFitUtils.util import loadRooFitUtils
    # load libraries
    loadRooFitUtils()

    # setup verbosity
    ROOT.RooFitUtils.Log.SetReportingLevel(ROOT.RooFitUtils.Log.FromString(args.get("loglevel","DEBUG")))
    if args.get("loglevel","DEBUG") == "DEBUG":
        ROOT.gErrorIgnoreLevel = 0
        try:
            ROOT.Minuit2.MnPrint.ClearFilter()
        except AttributeError:
            pass
        ROOT.Math.MinimizerOptions.SetDefaultPrintLevel(2)
    else:
        ROOT.Math.MinimizerOptions.SetDefaultPrintLevel(-1)

    if ROOT.Math.MinimizerOptions.DefaultPrintLevel() < 0: ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

    # Configuration of minimizer
    ROOT.Math.MinimizerOptions.SetDefaultMinimizer(args.get("minimizerType","Minuit2"), args.get("minimizerAlgo","Migrad"))
    ROOT.Math.MinimizerOptions.SetDefaultStrategy(args.get("defaultStrategy",1))

    # patches for HCombRoot
    ROOT.RooFitUtils.RooStarMomentMorphFix = args.get("fixCache",True)
    ROOT.RooFitUtils.RooMultiPdfFix = args.get("fixMulti",True)

def buildModel(args):
    # Load the model
    import ROOT

    model = ROOT.RooFitUtils.ExtendedModel("model", args["inFileName"], args["wsName"],
                                           args.get("modelConfigName","ModelConfig"), args.get("dataName","asimovData"), args.get("snapshot",""), args.get("pdfName","simPdf"),
                                           args.get("binnedLikelihood",True))

    pois = args.get("poi",None)
    if pois and model.GetModelConfig(): model.GetModelConfig().SetParametersOfInterest(",".join(pois))

    if args.get("penalty",False):
        npens = len(args["penalty"])
        ws = model.GetWorkspace()
        if npens > 0:
            for ipens in range(0, npens):
                pars = ROOT.RooArgList()
                allpars = args["penalty"][0][1].split(",")
                for par in allpars[0:len(allpars)]:
                    par = par.strip(" ") 
                    par = par.strip("\"") 
                    pars.add(ws.obj(par))
        
                name = "penalty_"+str(ipens)
                penaltyform = ROOT.RooFormulaVar(name, args["penalty"][ipens][0], pars)
                model.addPenalty(penaltyform)

    if args.get("penaltyfile",None):
        ws = model.GetWorkspace()
        ipens = 0
        for penline in open(args["penaltyfile"],"r"):
            ipens = ipens + 1
            lineparts = penline.strip("\n").split("|")
            lineparts = (lineparts[0], lineparts[1])
            pars = ROOT.RooArgList()
            allpars = lineparts[1].split(",")
            for par in allpars[0:len(allpars)]:
                par = par.strip(" ")
                pars.add(ws.obj(par))
            
            name = "penalty_"+str(ipens)
            penaltyform = ROOT.RooFormulaVar(name, lineparts[0], pars)
            model.addPenalty(penaltyform)

    if args.get("fixAllNP",False):          model.fixNuisanceParameters()
    if args.get("setInitialError",False):   model.setInitialErrors()
    if args.get("fixParameters",None):      model.fixParameters(",".join(args["fixParameters"]))
    if args.get("floatParameters",None):    model.floatParameters(",".join(args["floatParameters"]))
    if args.get("randomizeParameters",None):model.randomizeParameters(",".join(args["randomizeParameters"]))        

    model.fixParametersOfInterest()
    model.profileParameters(",".join(args.get("profile",[])))

    return model

def buildMinimizer(args,model):
    import ROOT
    from RooFitUtils.util import makelist

    ws = model.GetWorkspace()
    allparams = ROOT.RooArgSet()
    nuis = model.GetNuisanceParameters()
    pois = model.GetParametersOfInterest()
    globs = model.GetGlobalObservables()
    if nuis:
        ROOT.RooFitUtils.addArgSet(allparams, nuis)
    if pois:
        ROOT.RooFitUtils.addArgSet(allparams, pois)
    if globs:
        ROOT.RooFitUtils.addArgSet(allparams, globs)
    obs = model.GetObservables()

    if args.get("makeParameterSnapshots",False):
       # Save the snapshots of nominal parameters
        print("Saving nominal snapshots.")
        ws.saveSnapshot("nominalGlobs", globs)
        ws.saveSnapshot("nominalNuis", nuis)
        ws.saveSnapshot("nominalPois", pois)
        ws.saveSnapshot("nominalAll", allparams)

    # Collect POIs
    poiset = ROOT.RooArgSet()
    if args.get("pois",None):
        poinames = args["pois"]
    else:
        poinames = [ p.GetName() for p in makelist(pois) ]
    for poi in poinames:
        p = model.configureParameter(poi)
        if not p:
            raise(RuntimeError("unable to find parameter '{0:s}'".format(poi)))
        p.setConstant(False)
        poiset.add(p)
    pdf = model.GetPdf()

    findSigma = args.get("findSigma",False)
    constOpt = args.get("constOpt",2)
    if findSigma and constOpt > 0:
        print("deactivating constant term optimization: incompatible with findSigma")
        constOpt = 0

    argelems = [ROOT.RooFit.Minimizer(args.get("minimizerType","Minuit2"), args.get("minimizerAlgo","Migrad")),
                ROOT.RooFit.Strategy(args.get("defaultStrategy",2)),
                ROOT.RooFitUtils.ExtendedMinimizer.Eps(args.get("eps",1e-3)),
                ROOT.RooFitUtils.ExtendedMinimizer.ReuseMinimizer(args.get("reuseMinimizer",False)),
                ROOT.RooFitUtils.ExtendedMinimizer.ReuseNLL(args.get("reuseNll",False)),
                ROOT.RooFitUtils.ExtendedMinimizer.MaxIterations(50*pdf.getVariables().getSize()),
                ROOT.RooFitUtils.ExtendedMinimizer.MaxCalls(5000*pdf.getVariables().getSize()),                
                ROOT.RooFitUtils.ExtendedMinimizer.Verbose(verbosity[args.get("loglevel","DEBUG")]),
                ROOT.RooFit.NumCPU(args.get("numCPU",1), args.get("mpStrategy",3)),
                ROOT.RooFit.Offset(args.get("offsetting",True)),
                ROOT.RooFit.Optimize(constOpt),
                ROOT.RooFit.PrintLevel(args.get("printLevel",ROOT.Math.MinimizerOptions.DefaultPrintLevel())),
                ROOT.RooFit.Hesse(args.get("hesse",True)),
                ROOT.RooFit.Save()]

    if globs and globs.getSize() > 0:
        argelems.append(ROOT.RooFit.GlobalObservables(globs))
    elif nuis:
        argelems.append(ROOT.RooFit.Constrain(nuis))
    if args.get("batchMode"):
        argelems.append(ROOT.RooFit.BatchMode(args.get("batchMode")))
    if args.get("findSigma",False):
        argelems.append(ROOT.RooFitUtils.ExtendedMinimizer.FindSigma(poiset))
    elif args.get("minos",True):
        argelems.append(ROOT.RooFit.Minos(poiset))        

    from RooFitUtils.util import nodel
    nodel(poiset)
    nodel(argelems)
    arglist = ROOT.RooLinkedList()
    for arg in argelems: arglist.Add(arg)

    minimizer = ROOT.RooFitUtils.ExtendedMinimizer("minimizer", model, arglist)
    return minimizer


def fit(args,model,minimizer):
    from time import time
    from RooFitUtils.util import parsePoint,makelist,timestamp,linspace,vec,mkdir,union,names
    from RooFitUtils.io import writeResultJSON as writeResult
    import ROOT
    opt = ROOT.Math.MinimizerOptions.Default("Minuit2")
    opt.SetValue("StorageLevel",0)
    ROOT.Math.MinimizerOptions.SetDefaultExtraOptions(opt)


    # Collect POIs
    if args.get("pois",None):
        poinames = args["pois"]
        pois = [ model.GetWorkspace().var(p) for p in poinames ]
    else:
        pois = model.GetParametersOfInterest()
        poinames = names(pois)
    for poi in poinames:
        p = model.configureParameter(poi)

        if not p:
            raise(RuntimeError("unable to find parameter '{0:s}'".format(poi)))
        p.setConstant(False)

    if args.get("fit",True):
        start = time()
        if not args.get("dummy",False):
            minimizer.minimize()

        end = time()
        print("Fitting time: " + timestamp(end-start))
        minNll = minimizer.GetMinNll()
        print("NLL after minimisation: "+str(minNll))

        if args.get("makeParameterSnapshots",False):
            ws = model.GetWorkspace()
            allparams = ROOT.RooArgSet()
            nuis = model.GetNuisanceParameters()
            allparams.add(nuis)
            globs = model.GetGlobalObservables()
            allparams.add(globs)
            allparams.add(pois)
            obs = model.GetObservables()
            # Save the snapshots of nominal parameters
            print("Saving minimum snapshots.")
            ws.saveSnapshot("minimumGlobs", globs)
            ws.saveSnapshot("minimumNuis", nuis)
            ws.saveSnapshot("minimumPois", pois)
            ws.saveSnapshot("minimumAll", allparams)

    else:
        print("no fit requested")

    parnames = None
    coords = None
    scan = args.get("scan",False)
    if scan and len(scan)>0:
        from RooFitUtils.util import isstr
        if isstr(scan[0]): scan = [scan]
        from RooFitUtils.util import generateCoordsDict
        coordsdict = generateCoordsDict(scan)
        parnames = vec(sorted(coordsdict[0].keys()),"string")
        coords = vec([ vec([d[str(k)] for k in parnames],"double") for d in coordsdict],"vector<double>")

    if args.get("points",None) != None:
        with open(args["points"]) as infile:
            points = [ parsePoint(line) for line in infile if len(line)>0 ]
            parnames = vec(sorted(union([p.keys() for p in points])),"string")
            coords = vec( [ vec( [ point[str(p)] for p in parnames ] , "double") for point in points ], "vector<double>")

    if args.get("point",None) != None:
        point = parsePoint(args["point"])
        parnames = vec(sorted(point.keys()),"string")
        coords = vec( [ vec( [ point[str(p)] for p in parnames ] , "double") ], "vector<double>")

    if parnames and coords and not args.get("dummy",False):
        minimizer.scan(parnames,coords)
    else:
        print("no scan requested")

    result = minimizer.getResult()

    if result != None:
        outFileName = args.get("outFileName","output.root")
        if outFileName:
            import os
            outpath,outfile = os.path.split(outFileName)
            mkdir(outpath)
            with open(outFileName,'w') as out:
                writeResult(out,result,args.get("correlationMatrix",True))
            print("wrote output to "+outFileName)
            if result.min.fit:
                if args.get("writeResult",False):
                    import os
                    outpath,outfile = os.path.split(outFileName)
                    fitresfile = os.path.join(outpath,os.path.splitext(outfile)[0]+".root")
                    result.min.fit.SaveAs(fitresfile)
                if args.get("printResult",True):
                    for poi in makelist(pois):
                        p = result.min.fit.floatParsFinal().find(poi)
                        if p:
                            print("{:s} = {:g} +{:g} -{:g}".format(poi.GetName(),poi.getVal(),abs(poi.getErrorHi()),abs(poi.getErrorLo())))
            else:
                print("invalid fit result!")
        if not args.get("writeResult",False):
            print("no output requested")
    else:
        print("received invalid result")

    if args.get("outWsName",None):
        ws = model.GetWorkspace()
        ws.writeToFile(args["outWsName"])

    return result

def createScanJobs(args,arglist):
    from os.path import join as pjoin
    from RooFitUtils.util import stringify,makepoint,reconstructCall,generateCoordsDict,mkdir,concat
    from RooFitUtils.util import distributePointsAroundPoint,distributePointsAroundLine
    options = reconstructCall(args,arglist,["scan","findSigma","writeSubmit","writeSubmitPoints","refineScan","refineScanThresholds"])
    import sys
    submitCommand = concat([args["submitCommand"],sys.argv[0]]," ")    
    if args["refineScan"]:
        from RooFitUtils.io import collectresults
        prescans = {}
        if "scan" in args.keys() and args["scan"]:
            parnamelist = [ s[0] for s in args["scan"] ]
            ranges = [(float(s[2]),float(s[3])) for s in args["scan"]]
        else:
            parnamelist = None
            ranges = None
        collectresults(prescans,args["refineScan"],"dummy",filterScans=parnamelist)
        from RooFitUtils.interpolate import findcontours
        coords = []
        for label,scan in prescans["scans"].items():
            if not parnamelist: parnamelist = label
            for labels,points in scan.items():
                # for now, use as many points for the new scan as for the old one
                npoints = args["refineScanPoints"]
                if len(parnamelist) == 2:
                    percent_thresholds = args["refineScanThresholds"]
                    from RooFitUtils.util import getThreshold
                    thresholds = [ getThreshold(p/100,2)*0.5 for p in percent_thresholds ]
                    contours,minimum = findcontours(points,values=thresholds,smooth=False,npoints=100,ranges=ranges)
                    # for now, assign 10% of the points to the minimum, divide the rest evenly among the contours
                    if len(contours) == 0:
                        raise RuntimeError("cannot refine: no contours found!")
                    nEach = int(1 * npoints / len(contours))
                    for contour in contours:
                        for graph in contour:
                            distributePointsAroundLine(parnamelist,coords,graph,nEach,args["refineScanSpread"])
                    # the distpar argument needs to be tuned to fit the coodinate sytem, TODO: come up with a smart way of guessing it#
                    from RooFitUtils.util import inarea
                    if not ranges or inarea(minimum,ranges):
                        distributePointsAroundPoint(parnamelist,coords,minimum,int(0.1*npoints),0.05*args["refineScanSpread"])
                else:
                    percent_thresholds = args["refineScanThresholds"]
                    thresholds = [ getThreshold(p/100,1)*0.5 for p in percent_thresholds ]                    
                    for t in thresholds:
                        from RooFitUtils.interpolate import findcrossings
                        cv,down,up = findcrossings(points,t)
                        distributePointsAroundPoint(parnamelist,coords,down,npoints/4,0.1)
                        distributePointsAroundPoint(parnamelist,coords,up,npoints/4,0.1)
    elif args["scan"]:
        coords = generateCoordsDict(args["scan"])
    idx = 0
    import os
    outpath = args["writeSubmit"]
    mkdir(outpath)
    outfile = "jobs.txt"

    pointspath = pjoin(outpath,"coords_0.txt")
    from RooFitUtils.util import clearfile
    clearfile(pjoin(outpath,outfile))
    clearfile(pointspath)

    idx = 0
    if not args.get("outFileName",None):
        print("output file name mandatory for use of batch scanning!")
        exit(0)
    with open(pjoin(outpath,outfile),"w") as jobs:
        for coord in coords:
            if  idx % args["writeSubmitPoints"] == 0:
                pointspath =outpath+"/coords" +"_"+str(idx)+".txt"
                clearfile(pointspath)
                from os.path import splitext,exists
                outfilename,outfileext = splitext(args["outFileName"])
                ofname = outfilename + ".part"+str(idx)+outfileext
                options[" --no-findSigma --points"]=pointspath
                if "-o" in options.keys(): del options["-o"]
                options["--output"]=ofname
                cmd = " ".join([k+" "+stringify(v) for k,v in options.items()])
                jobs.write(submitCommand+" "+cmd+"\n")
            idx = idx + 1
            with open(pointspath,"a") as coordlist:
                point = makepoint(coord)
                coordlist.write(point+"\n")

    print("wrote "+args["writeSubmit"])


def createImpactJobs(args,arglist):
    from os.path import join as pjoin
    from RooFitUtils.util import stringify,makepoint,reconstructCall,mkdir,names,concat,makelist
    options = reconstructCall(args,arglist,["impacts","writeSubmit","writeSubmitPoints","refineScan","refineScanThresholds","findSigma"])
    model = buildModel(args)    
    import sys
    submitCommand = concat([args["submitCommand"],sys.argv[0]]," ")    
    pars = []
    for group in args["impacts"]:
        import re
        regex = re.compile(group)
        for par in makelist(model.GetNuisanceParameters()):
            if regex.match(par.GetName()):
                if abs(par.getErrorHi()) > 0 and abs(par.getErrorLo()) > 0:
                    pars.append(par)
                else:
                    print("warning: cannot produce impact of parameter '"+par.GetName()+"' without valid error range set!")
    outpath = args["writeSubmit"]
    mkdir(outpath)
    outfile = "jobs.txt"
    from os.path import join as pjoin

    options["--no-findSigma"] = ""
    
    with open(pjoin(outpath,outfile),"w") as jobs:
        if "-o" in options.keys(): del options["-o"]        
        options["--output"]=pjoin(outpath,"impact.nominal.json")
        cmd = " ".join([k+" "+stringify(v) for k,v in options.items()])
        jobs.write(submitCommand+" "+cmd+"\n")
        for par in pars:
            options["--fix"]="{:s}={:f}".format(par.GetName(),par.getVal()+abs(par.getErrorHi()))
            options["--output"]=pjoin(outpath,"impact."+par.GetName()+".up.json")
            cmd = " ".join([k+" "+stringify(v) for k,v in options.items()])                    
            jobs.write(submitCommand+" "+cmd+"\n")
            options["--fix"]="{:s}={:f}".format(par.GetName(),par.getVal()-abs(par.getErrorLo()))
            options["--output"]=pjoin(outpath,"impact."+par.GetName()+".dn.json")
            cmd = " ".join([k+" "+stringify(v) for k,v in options.items()])                    
            jobs.write(submitCommand+" "+cmd+"\n")            
    print("wrote "+args["writeSubmit"])


def createBreakdownJobs(args,arglist):
    from os.path import join as pjoin
    from RooFitUtils.util import names,reconstructCall,stringify,mkdir,concat
    options = reconstructCall(args,arglist,["breakdown","findSigma","writeSubmit","writeSubmitPoints"])
    model = buildModel(args)
    import sys
    submitCommand = concat([args["submitCommand"],sys.argv[0]]," ")
    import re
    groups = {}
    for group in args["breakdown"]:
        pars = []
        for ex in group:
            regex = re.compile(ex)
            if model.GetNuisanceParameters():
                for parameter in model.GetNuisanceParameters():
                    if regex.match(parameter.GetName()):
                        pars.append(parameter)
        groupname = "_".join([ ex.replace("_","").replace("*","") for ex in group])
        groups[groupname] = names(pars)
    outpath = args["writeSubmit"]
    mkdir(outpath)
    outfile = "jobs.txt"
    with open(pjoin(outpath,outfile),"w") as jobs:
        options["--findSigma"] = ""
        if "-o" in options.keys(): del options["-o"]        
        options["--output"]=pjoin(outpath,"breakdown.nominal.json")
        cmd = " ".join([k+" "+stringify(v) for k,v in options.items()])
        jobs.write(submitCommand+" "+cmd+"\n")
        for groupname,parnames in groups.items():
            options["--fix"] = " ".join(parnames)
            options["--output"]=pjoin(outpath,"breakdown."+groupname+".json")
            cmd = " ".join([k+" "+stringify(v) for k,v in options.items()])                    
            jobs.write(submitCommand+" "+cmd+"\n")
    print("wrote "+args["writeSubmit"])


