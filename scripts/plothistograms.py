def main(args):
    from RooFitUtils.io import texdict,collecthistograms
    from RooFitUtils.util import flattened,keys
    from RooFitUtils.pgfplotter import plotBars
    parnames = []
    if args.plot:
        import yaml        
        with open(args.plot,"rt") as infile:
            plot = yaml.load(infile,Loader=yaml.FullLoader)
            parnames = sorted(list(set(flattened([keys(panel["data"]) for panel in plot.values()]))))
            catnames = sorted(list(set(flattened([keys(panel["bins"]) for panel in plot.values()]))))           
    else:
        if args.parameters:
            parameters = readcsv2dict(args.parameters)
            texdict({p["name"]:p["label"] for p in parameters.values()})
            parnames = parameters.keys()
        if args.categories:
            categories = readcsv2dict(args.categories)
            texdict({c["name"]:c["label"] for c in categories.values()})
            catnames = categories.keys()
        plot = {"main":{"bins":categories,"data":parameters}}
    histograms = {}    
    if args.infile.endswith(".root"):
        collecthistograms(histograms,{"input":args.infile,"labels":catnames},parnames)
    elif args.infile.endswith(".yml"):
        import yaml
        with open(args.infile,"rt") as infile:
            inputs = yaml.load(infile,Loader=yaml.FullLoader)
            for key,cfg in inputs.items():
                collecthistograms(histograms,cfg,parnames)
    plotBars(args.output,histograms,plot,xwidth=args.width)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="plot impacts")
    parser.add_argument("--output",help="output file name",required=True)
    parser.add_argument("--input",dest="infile",help="input file name (root or yml)",required=True)
    parser.add_argument("--plot",help="yml file configuring the plot layout")    
    parser.add_argument("--parameters",help="styles of parameters")
    parser.add_argument("--categories",help="labels of categories")
    parser.add_argument("--width",help="bar width",default="1em")    
    args = parser.parse_args()
    main(args)
    
