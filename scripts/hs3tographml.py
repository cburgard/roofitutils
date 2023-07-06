def collect_strings(d,skipName):
    """
    Recursively collects unique strings from a dictionary.

    Args:
        d (dict): The input dictionary.
        skipName (bool): Whether to skip collecting strings from the "name" key.

    Returns:
        set: A set of unique strings collected from the dictionary.
    """
    values = set()
    for k,v in d.items():
        if skipName and k == "name": continue        
        elif type(v) == list:
            for e in v:
                if type(e) == str:
                    values.add(e)
                if type(e) == dict:
                    values = values.union(collect_strings(e,False))
        elif type(v) == dict:
            values = values.union(collect_strings(v,False))
        else:
            values.add(str(v))
    return values
    
def fill_graph(model,element,elements):
    """
    Fills a graph model recursively with the provided element and all its dependents

    Args:
        model (dict): The graph model.
        element (dict): The current element.
        elements (dict): Dictionary of elements.

    Returns:
        None
    """
    name = element["name"]
    if name in model.keys():
        return
    model[name] = set()
    values = collect_strings(element,True)

    for k in elements.keys():
        for v in values:
            if k in v:
                model[element["name"]].add(k)
                fill_graph(model,elements[k],elements)

def main(args):
    """
    The main entry point of the program.

    Args:
        args (argparse.Namespace): Command-line arguments.

    Returns:
        None
    """
    import json

    with open(args.infile, 'r') as f:
        data = json.load(f)

    likelihood = None
    if not args.likelihood:
        print("please select a likelihood with the -l/--likelihood option, available options for this model are:")
    for lh in data["likelihoods"]:
        if args.likelihood:
            if lh["name"] == args.likelihood:
                likelihood = lh
        else:
            print("  "+lh["name"])    
    if not likelihood:
        print("no likelihood selected, exiting")

    elements = {}
    for e in data["functions"]:
        elements[e["name"]] = e
    for e in data["distributions"]:
        elements[e["name"]] = e
    for d in data["domains"]:
        for v in d["axes"]:
            elements[v["name"]] = v

    model = {}

    dists = []
    for dist in likelihood["distributions"]:
        dists.append(dist)
        for d in data["distributions"]:
            if d["name"] == dist:
                fill_graph(model,d,elements)
    for dist in likelihood["aux_distributions"]:
        dists.append(dist)
        for d in data["distributions"]:
            if d["name"] == dist:
                fill_graph(model,d,elements)        

    import xml.etree.ElementTree as ET
    data = ET.Element("graphml")
    data.set("xmlns","http://graphml.graphdrawing.org/xmlns")
    tree = ET.ElementTree(data)

    graph = ET.SubElement(data,"graph")
    graph.set("id","model")
    graph.set("edgedefault","directed")

    # nodes
    lh = ET.SubElement(graph,"node")
    lh.set("id",likelihood["name"])
    for name in model.keys():
        node = ET.SubElement(graph,"node")
        node.set("id",str(name))

    # edges
    for client,serverlist in model.items():
        for server in serverlist:
            edge = ET.SubElement(graph,"edge")
            edge.set("source",str(server))
            edge.set("target",str(client))
    for dist in dists:
        edge = ET.SubElement(graph,"edge")
        edge.set("source",str(dist))
        edge.set("target",str(likelihood["name"]))        

    with open(args.outfile,"wb") as outfile:
        tree.write(outfile, encoding='utf-8', xml_declaration = True)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="convert a HS3 JSON file to GraphML")
    parser.add_argument("-i","--input",metavar="model.json",help="input JSON file",dest="infile",required=True)
    parser.add_argument("-o","--output",metavar="model.gml",help="GraphML output file",dest="outfile",required=True)
    parser.add_argument("--likelihood","-l",metavar="mylikelihood",help="name of the likleihood to use",dest="likelihood")    
    main(parser.parse_args())
    
