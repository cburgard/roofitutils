#!/bin/env python

def stripname(name,stripcomponents):
    while True:
        anymatch = False
        for c in stripcomponents:
            if name.endswith(c):
                name = name[:-len(c)].strip("_")
                anymatch = True
                break
            if name.startswith(c):
                name = name[len(c):].strip("_")
                anymatch = True
                break            
        if not anymatch:
            return name

def main(args):
    import re
    import json
    with open(args.input,"rt") as input:
        input_json = json.load(input)

    measurement = {"name":"meas","config":{"parameters":[],"poi":args.poi}}
    output_json = {"channels":[],"measurements":[measurement],"observations":[],"version":"1.0.0"}

    nps = set()
    nfs = set()    
    
    simpdf = [ pdf for pdf in input_json["pdfs"].values() if pdf["type"] == "simultaneous" ][0]
    for key,channel in simpdf["channels"].items():
        out_channel = {"name":stripname(key,args.strip_name_components),"samples":[]}
        output_json["channels"].append(out_channel)
        for key,sample in channel["samples"].items():
            out_sample = {"name":stripname(key,args.strip_name_components),"modifiers":[]}
            out_channel["samples"].append(out_sample)                    
            out_sample["data"] = sample["data"]["counts"]
            if "normFactors" in sample.keys():
                for nf in sample["normFactors"]:
                    if any([ re.match(b,nf) for b in args.nf_blacklist]): continue                    
                    nfs.add(nf)                    
                    out_sample["modifiers"].append({"type":"normfactor","name":nf,"data":None})
                    
            modifiernames = set()
            if "histogramSystematics" in sample.keys():
                for key,histoSys in sample["histogramSystematics"].items():
                    nps.add(key)
                    if args.defactorize and "overallSystematics" in sample.keys() and key in sample["overallSystematics"].keys():
                        overallSys = sample["overallSystematics"][key]
                        modifiernames.add(key)
                        out_sample["modifiers"].append({"type":"histosys","name":key,"data":{
                            "hi_data":[ v * overallSys["high"] for v in histoSys["dataHigh"]["counts"]],
                            "lo_data":[ v * overallSys["low"]  for v in histoSys["dataLow"]["counts"]]
                        }})
                    else:
                        out_sample["modifiers"].append({"type":"histosys","name":key,"data":{"hi_data":histoSys["dataHigh"]["counts"],"lo_data":histoSys["dataLow"]["counts"]}})
            if "overallSystematics" in sample.keys():
                for key,overallSys in sample["overallSystematics"].items():
                    if key in modifiernames: continue
                    nps.add(key)
                    out_sample["modifiers"].append({"type":"normsys","name":key,"data":{"hi":overallSys["high"],"lo":overallSys["low"]}})

    for par in nps:
        measurement["config"]["parameters"].append({"bounds":[[-1,1]],"fixed":False,"inits":[0.],"name":par})
    for par in nfs:
        if any([ re.match(b,par) for b in args.nf_blacklist]): continue
        measurement["config"]["parameters"].append({"bounds":[[0,50]],"fixed":False,"inits":[1.],"name":par})        
                    
    data = input_json["data"]["obsData"]
    for key,channel in data.items():
        if type(channel) != dict:
            continue
        if "counts" in channel.keys():
            output_json["observations"].append({"data":channel["counts"],"name":stripname(key,args.strip_name_components)})
        else:
            output_json["observations"].append({"data":channel["weights"],"name":stripname(key,args.strip_name_components)})            
                    
    with open(args.output,"wt") as output:
        json.dump(output_json,output,sort_keys=True)
                
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="format converter between ROOT JSON and pyhf JSON")
    parser.add_argument("--poi",default=None)
    parser.add_argument("--nf-blacklist",nargs="+",default=["binWidth_.*"],dest="nf_blacklist")
    parser.add_argument("--strip-name-components",nargs="+",default=[])
    parser.add_argument("--defactorize",action="store_true",help="merge OverallSys and HisotSys of the same name",default=False)
    parser.add_argument("input",help="RooFit JSON file")
    parser.add_argument("output",help="pyhf JSON output")
    args = parser.parse_args()
    main(args)
    
    
