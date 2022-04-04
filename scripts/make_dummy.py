#!/bin/env python

def pseudo_poisson(v):
    from random import gauss
    from math import sqrt
    return gauss(v, sqrt(v))

def randomize_data(d):
    from numpy.random import uniform    
    if "counts" in d.keys():
        d["counts"] = [ pseudo_poisson(abs(v)) for v in d["counts"] ]
    if "weights" in d.keys():
        d["weights"] = [ pseudo_poisson(abs(v)) for v in d["weights"] ]
    if "errors" in d.keys():
        d["errors"] = [ uniform(0.5*v,1.5*v) for v in d["errors"] ]        

        
def randomize_embedded(d):
    from numpy.random import uniform
    if isinstance(d, list):
        for v in d:
            randomize_embedded(v)
    elif isinstance(d,dict):
        for k, v in d.items():
            if k == "data" or k == "dataHigh" or k=="dataLow":
                randomize_data(v)
            elif isinstance(v,dict) or isinstance(v,list):
                randomize_embedded(v)
            elif k == "low" or k == "high" or k == "value":
                if v > 1:
                    d[k] = uniform(v,1+2*(v-1))
                else:
                    d[k] = uniform(1-2*(1-v),1)

            
def remove_keys(d,keys):
    if not isinstance(d, (dict, list)):
        return d
    if isinstance(d, list):
        return [remove_keys(v,keys) for v in d if not v  in keys]
    return {k: remove_keys(v,keys) for k, v in d.items()
            if k not in keys }

def main(args):
    import json
    with open(args.input,"rt") as infile:
        ws = json.load(infile)
        ws = remove_keys(ws,["factory_tag","origName","SnapShot_ExtRefClone"])
        if "data" in ws.keys():
            for d in ws["data"].values():
                for region in d.values():
                    if not type(region) == dict:
                        continue
                    randomize_data(region)
        if args.randomize_embedded:
            if "pdfs" in ws.keys():
                for pdf in ws["pdfs"].values():
                    randomize_embedded(pdf)
            if "functions" in ws.keys():
                for pdf in ws["functions"].values():
                    randomize_embedded(pdf)                

    with open(args.output,"wt") as outfile:
        json.dump(ws,outfile)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="convert a JSON workspace to a dummy version of it")
    parser.add_argument("--randomize-embedded",default=False,action="store_true",help="randomize embedded data")
    parser.add_argument("input",help="input workspace in JSON format")
    parser.add_argument("output",help="output workspace in JSON format")    
    main(parser.parse_args())
