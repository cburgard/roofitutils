def setup(args):
    # general setup and loading of modules
    import ROOT

    from RooFitUtils.util import loadRooFitUtils
    # load libraries
    loadRooFitUtils()

    # setup verbosity
    ROOT.RooFitUtils.Log.SetReportingLevel(ROOT.RooFitUtils.Log.FromString(args.loglevel))
    if args.loglevel == "DEBUG":
        ROOT.Math.MinimizerOptions.SetDefaultPrintLevel(2)
    else:
        ROOT.Math.MinimizerOptions.SetDefaultPrintLevel(-1)

    if ROOT.Math.MinimizerOptions.DefaultPrintLevel() < 0: ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

    # Configuration of minimizer
    ROOT.Math.MinimizerOptions.SetDefaultMinimizer(args.minimizerType, args.minimizerAlgo)
    ROOT.Math.MinimizerOptions.SetDefaultStrategy(args.defaultStrategy)

    # patches for HCombRoot
    ROOT.RooFitUtils.RooStarMomentMorphFix = args.fixCache
    ROOT.RooFitUtils.RooMultiPdfFix = args.fixMulti

def buildModel(args):
    # Load the model
    import ROOT


    model = ROOT.RooFitUtils.ExtendedModel("model", args.inFileName, args.wsName,
                                           args.modelConfigName, args.dataName, args.snapshot,
                                           args.binnedLikelihood, ROOT.RooArgSet(), "pdf_")

    if args.penalty:
        npens = len(args.penalty)
        ws = model.GetWorkspace()
        if npens > 0:
            for ipens in range(0, npens):
                pars = ROOT.RooArgList()
                allpars = args.penalty[0][1].split(",")
                for par in allpars[0:len(allpars)]:
                    par = par.strip(" ") 
                    par = par.strip("\"") 
                    ws.obj(par).Print()
                    pars.add(ws.obj(par))
        
                name = "penalty_"+str(ipens)
                penaltyform = ROOT.RooFormulaVar(name, args.penalty[ipens][0], pars)
                model.addPenalty(penaltyform)

    if args.penaltyfile:
        ws = model.GetWorkspace()
        ipens = 0
        for penline in open(args.penaltyfile,"r"):
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

    if args.fixAllNP:          model.fixNuisanceParameters()
    if args.setInitialError:   model.setInitialErrors()
    if args.fixParameters:     model.fixParameters(",".join(args.fixParameters))
    if args.floatParameters:   model.floatParameters(",".join(args.fixParameters))    

    model.fixParametersOfInterest()
    model.profileParameters(",".join(args.profile))
 
    return model

def buildMinimizer(args,model):
    import ROOT
    from RooFitUtils.util import makelist
    ws = model.GetWorkspace()
    mc = model.GetModelConfig()
    allparams = ROOT.RooArgSet()
    nuis = model.GetNuisanceParameters()
    pois = model.GetParametersOfInterest()
    globs = model.GetGlobalObservables()
    ROOT.RooFitUtils.addArgSet(allparams, nuis)
    ROOT.RooFitUtils.addArgSet(allparams, pois)
    ROOT.RooFitUtils.addArgSet(allparams, globs)
    obs = model.GetObservables()
    if args.makeParameterSnapshots:
       # Save the snapshots of nominal parameters
        print("Saving nominal snapshots.")
        ws.saveSnapshot("nominalGlobs", globs)
        ws.saveSnapshot("nominalNuis", nuis)
        ws.saveSnapshot("nominalPois", pois)
        ws.saveSnapshot("nominalAll", allparams)


    # Collect POIs
    pois = model.GetParametersOfInterest()
    poiset = ROOT.RooArgSet()
    if args.pois:
        poinames = args.pois
    else:
        poinames = [ p.GetName() for p in makelist(pois) ]
    for poi in poinames:
        p = model.configureParameter(poi)
        if not p:
            raise(RuntimeError("unable to find parameter '{0:s}'".format(poi)))
        p.setConstant(False)
        poiset.add(p)

    pdf = model.GetPdf()

    argelems = [ROOT.RooFit.Minimizer(args.minimizerType, args.minimizerAlgo),
                ROOT.RooFit.Strategy(args.defaultStrategy),
                ROOT.RooFitUtils.ExtendedMinimizer.Eps(args.eps),
                ROOT.RooFitUtils.ExtendedMinimizer.ReuseMinimizer(args.reuseMinimizer),
                ROOT.RooFitUtils.ExtendedMinimizer.ReuseNLL(args.reuseNll),
                ROOT.RooFitUtils.ExtendedMinimizer.MaxCalls(5000*pdf.getVariables().getSize()),
                ROOT.RooFit.Constrain(nuis),
                ROOT.RooFit.GlobalObservables(globs),
                ROOT.RooFit.NumCPU(args.numCPU, args.mpStrategy),
                ROOT.RooFit.Offset(args.offsetting),
                ROOT.RooFit.Optimize(args.constOpt),
                ROOT.RooFit.PrintLevel(args.printLevel),
                ROOT.RooFit.Precision(args.precision),
                ROOT.RooFit.Hesse(args.hesse),
                ROOT.RooFit.Save()]
    if args.findSigma:
        argelems.append(ROOT.RooFitUtils.ExtendedMinimizer.Scan(poiset))

    from RooFitUtils.util import nodel
    nodel(poiset)
    nodel(argelems)
    arglist = ROOT.RooLinkedList()
    for arg in argelems: arglist.Add(arg)

    minimizer = ROOT.RooFitUtils.ExtendedMinimizer("minimizer", model,arglist)
    return minimizer


def fit(args,model,minimizer):
    from time import time
    from RooFitUtils.util import parsePoint,makelist,timestamp,linspace,vec,mkdir,union
    from RooFitUtils.util import generateCoordsDict
    from RooFitUtils.io import writeResultJSON as writeResult
    import ROOT
    opt = ROOT.Math.MinimizerOptions.Default("Minuit2")
    opt.SetValue("StorageLevel",0)
    ROOT.Math.MinimizerOptions.SetDefaultExtraOptions(opt)


    # Collect POIs
    pois = model.GetParametersOfInterest()
    if args.pois:
        poinames = args.pois
    else:
        poinames = [ p.GetName() for p in makelist(pois) ]
    for poi in poinames:
        p = model.configureParameter(poi)

        if not p:
            raise(RuntimeError("unable to find parameter '{0:s}'".format(poi)))
        p.setConstant(False)

    if args.fit:
        start = time()
        if not args.dummy:
            minimizer.minimize()

        end = time()
        print("Fitting time: " + timestamp(end-start))
        minNll = minimizer.GetMinNll()
        print("NLL after minimisation: "+str(minNll))

        if args.makeParameterSnapshots:
            ws = model.GetWorkspace()
            mc = model.GetModelConfig()
            allparams = ROOT.RooArgSet()
            nuis = model.GetNuisanceParameters()
            allparams.add(nuis)
            globs = model.GetGlobalObservables()
            allparams.add(globs)
            pois = model.GetParametersOfInterest()
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
    if args.scan:
        val = args.scan[0]
        parname = val[0]
        parrange = linspace(float(val[2]),float(val[3]),int(val[1]))
        parnames = vec([parname],"string")
        coords = vec([ vec([val],"double") for val in parrange],"vector<double>")
        coordsdict = generateCoordsDict(args.scan)
        parnames = vec(sorted(coordsdict[0].keys()),"string")
        coords = vec([ vec([d[str(k)] for k in parnames],"double") for d in coordsdict],"vector<double>")

    if args.points != None:
        with open(args.points) as infile:
            points = [ parsePoint(line) for line in infile if len(line)>0 ]
            parnames = vec(sorted(union([p.keys() for p in points])),"string")
            coords = vec( [ vec( [ point[str(p)] for p in parnames ] , "double") for point in points ], "vector<double>")

    if args.point != None:
        point = parsePoint(args.point)
        parnames = vec(sorted(point.keys()),"string")
        coords = vec( [ vec( [ point[str(p)] for p in parnames ] , "double") ], "vector<double>")

    if parnames and coords and not args.dummy:
        minimizer.scan(parnames,coords)
    else:
        print("no scan requested")

    result = minimizer.getResult()

    if result:
        if args.outFileName:
            import os
            outpath,outfile = os.path.split(args.outFileName)
            mkdir(outpath)
            with open(args.outFileName,'w') as out:
                writeResult(out,result,args.correlationMatrix)
            print("wrote output to "+args.outFileName)
        if args.writeResult:
            outpath,outfile = os.path.split(args.outFileName.replace(".txt",""))
            result.fit.SaveAs(outfile+"_fitresult.root")
        if not args.writeResult and not args.outFileName:
            print("no output requested")
    else:
        print("received invalid result")

    if args.outWsName:
        ws = model.GetWorkspace()
        ws.writeToFile(args.outWsName)


def createScanJobs(args,arglist):
    from os.path import join as pjoin
    from RooFitUtils.util import stringify,makepoint,reconstructCall,generateCoordsDict,mkdir,concat
    from RooFitUtils.util import distributePointsAroundPoint,distributePointsAroundLine
    options = reconstructCall(args,arglist,["scan","findSigma","writeSubmit","writeSubmitPoints","refineScan","refineScanThresholds"])
    import sys
    submitCommand = concat([args.submitCommand,sys.argv[0]]," ")    
    if args.refineScan:
        from RooFitUtils.io import collectresults
        prescans = {}
        preresults = {}
        collectresults(prescans,preresults,args.refineScan,"dummy")
        from RooFitUtils.interpolate import findcontours
        coords = []
        for parnamelist,scan in prescans.items():
            for labels,points in scan.items():
                # for now, use as many points for the new scan as for the old one
                npoints = 1000
                if len(parnamelist) == 2:
                    # 1 sigma (=68.26895% CL):  2.296
                    # 2 sigma (=95.44997% CL):  6.180
#                    thresholds = [0.5*2.296,0.5*6.180]
#                    thresholds = [0.5*2.28]
                    if args.refineScanThresholds:
                        thresholds = args.refineScanThresholds
                    else:
                        thresholds = [0.5*2.28,0.5*5.99]
                    contours,minimum = findcontours(points,thresholds,False)
                    # for now, assign 10% of the points to the minimum, divide the rest evenly among the contours
                    nEach = int(1 * npoints / len(contours))
                    for contour in contours:
                        for graph in contour:
                            distributePointsAroundLine(parnamelist,coords,graph,nEach)
                    # the distpar argument needs to be tuned to fit the coodinate sytem, TODO: come up with a smart way of guessing it
                    #distributePointsAroundPoint(parnamelist,coords,minimum,int(0.1*npoints),0.001)
                else:
                    if args.refineScanThresholds:
                        thresholds = args.refineScanThresholds
                    else:
                        thresholds = [0.5,2]
                    for t in thresholds:
                        cv,down,up = findcrossings(points,t)
                        distributePointsAroundPoint(parnamelist,coords,down,npoints/4,0.1)
                        distributePointsAroundPoint(parnamelist,coords,up,npoints/4,0.1)
    elif args.scan:
        coords = generateCoordsDict(args.scan)
    idx = 0
    import os
    outpath = args.writeSubmit
    mkdir(outpath)
    outfile = "jobs.txt"

    pointspath = pjoin(outpath,"coords_0.txt")
    from RooFitUtils.util import clearfile
    clearfile(pjoin(outpath,outfile))
    clearfile(pointspath)

    idx = 0
    if not args.outFileName:
        print("output file name mandatory for use of batch scanning!")
        exit(0)
    with open(pjoin(outpath,outfile),"w") as jobs:
        for coord in coords:
            if  idx % args.writeSubmitPoints == 0:
                pointspath =outpath+"/coords" +"_"+str(idx)+".txt"
                clearfile(pointspath)
                options[" --no-findSigma --points"]=pointspath
                options["--output"]=args.outFileName+".part"+str(idx)
                cmd = " ".join([k+" "+stringify(v) for k,v in options.items()])
                if not os.path.exists(args.outFileName+".part"+str(idx)):
                    jobs.write(submitCommand+" "+cmd+"\n")
            idx = idx + 1
            with open(pointspath,"a") as coordlist:
                point = makepoint(coord)
                coordlist.write(point+"\n")

    print("wrote "+args.writeSubmit)


def createImpactJobs(args,arglist):
    from os.path import join as pjoin
    from RooFitUtils.util import stringify,makepoint,reconstructCall,mkdir,names,concat
    options = reconstructCall(args,arglist,["impacts","writeSubmit","writeSubmitPoints","refineScan","refineScanThresholds","findSigma"])
    model = buildModel(args)    
    import sys
    submitCommand = concat([args.submitCommand,sys.argv[0]]," ")    
    pars = []
    for group in args.impacts:
        import re
        regex = re.compile(group)
        for par in model.GetNuisanceParameters():
            if regex.match(par.GetName()):
                if abs(par.getErrorHi()) > 0 and abs(par.getErrorLo()) > 0:
                    pars.append(par)
                else:
                    print("warning: cannot produce impact of parameter '"+par.GetName()+"' without valid error range set!")
    outpath = args.writeSubmit
    mkdir(outpath)
    outfile = "jobs.txt"
    from os.path import join as pjoin

    options["--no-findSigma"] = ""
    
    with open(pjoin(outpath,outfile),"w") as jobs:
        options["--output"]=pjoin(outpath,"impact.nominal.txt")
        cmd = " ".join([k+" "+stringify(v) for k,v in options.items()])
        jobs.write(submitCommand+" "+cmd+"\n")
        for par in pars:
            options["--fix"]="{:s}={:f}".format(par.GetName(),par.getVal()+abs(par.getErrorHi()))
            options["--output"]=pjoin(outpath,"impact."+par.GetName()+".up.txt")
            cmd = " ".join([k+" "+stringify(v) for k,v in options.items()])                    
            jobs.write(submitCommand+" "+cmd+"\n")
            options["--fix"]="{:s}={:f}".format(par.GetName(),par.getVal()-abs(par.getErrorLo()))
            options["--output"]=pjoin(outpath,"impact."+par.GetName()+".dn.txt")
            cmd = " ".join([k+" "+stringify(v) for k,v in options.items()])                    
            jobs.write(submitCommand+" "+cmd+"\n")            
    print("wrote "+args.writeSubmit)


def createBreakdownJobs(args,arglist):
    from os.path import join as pjoin
    from RooFitUtils.util import names,reconstructCall,stringify,mkdir,concat
    options = reconstructCall(args,arglist,["breadkown","findSigma","writeSubmit","writeSubmitPoints"])
    model = buildModel(args)
    import sys
    submitCommand = concat([args.submitCommand,sys.argv[0]]," ")
    import re
    groups = {}
    for group in args.breakdown:
        pars = []
        for ex in group:
            regex = re.compile(ex)
            for parameter in model.GetNuisanceParameters():
                if regex.match(parameter.GetName()):
                    pars.append(parameter)
        groupname = "_".join([ ex.replace("_","").replace("*","") for ex in group])
        groups[groupname] = names(pars)
    outpath = args.writeSubmit
    mkdir(outpath)
    outfile = "jobs.txt"
    with open(pjoin(outpath,outfile),"w") as jobs:
        options["--findSigma"] = ""
        options["--output"]=pjoin(outpath,"breakdown.nominal.txt")
        cmd = " ".join([k+" "+stringify(v) for k,v in options.items()])
        jobs.write(submitCommand+" "+cmd+"\n")
        for groupname,parnames in groups.items():
            options["--fix"] = " ".join(parnames)
            options["--output"]=pjoin(outpath,"breakdown."+groupname+".txt")
            cmd = " ".join([k+" "+stringify(v) for k,v in options.items()])                    
            jobs.write(submitCommand+" "+cmd+"\n")
    print("wrote "+args.writeSubmit)


