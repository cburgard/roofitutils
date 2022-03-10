#!/usr/bin/env python

prefix = "{http://www.w3.org/2000/svg}"

def collect(elem):
    elements = {}
    for c in elem:
        tag = c.tag[len(prefix):]
        if tag == "g":
            elements.update(collect(c))
        else:
            elements[c.attrib["id"]] = c.attrib
    return elements

def parsecoords(coordsstring):
    instructions = coordsstring.split()
    mode = None
    lastx = 0
    lasty = 0
    points = []
    for i in instructions:
        if i == "m" or i == "l" or i == "c":
            mode = "m"
        elif i == "M" or i == "L" or i == "C":
            mode = "M"
        elif i == "h":
            mode = "h"
        elif i == "H":
            mode = "H"
        elif i == "v":
            mode = "v"
        elif i == "V":
            mode = "V"                                    
        elif i == "z":
            mode = "z"
        else:
            if mode == "M":
                points.append(tuple(map(float,i.split(","))))
            if mode == "h":
                points.append((float(i)+lastx,lasty))
            if mode == "H":
                points.append((float(i),lasty))
            if mode == "v":
                points.append((lastx,float(i)+lasty))
            if mode == "V":
                points.append((lastx,float(i)))
            if mode == "m":
                x,y = map(float,i.split(","))
                x += lastx
                y += lasty
                lastx = x
                lasty = y
                points.append((x,y))
    if mode == "z":
        points.append(points[0])
    return points

def average(thelist):
    return (sum(map(lambda x:x[0],thelist))/len(thelist),sum(map(lambda x:x[1],thelist))/len(thelist))

def transform(plot):
    transformed = {}
    xunit = plot["xunit"]
    yunit = plot["yunit"]    
    for key,coords in plot["curves"].items():
        transformed[key] = [(x/xunit,y/yunit) for x,y in coords]
    for key,coords in plot["markers"].items():
        transformed[key] = (coords[0]/xunit,coords[1]/yunit)
    return transformed

def export(plot,outfile):
    import json
    with open(outfile, 'w') as fout:
        json_dumps_str = json.dumps(plot, indent=4)
        print(json_dumps_str, file=fout)

def identify(elems):
    plot = {}
    plot["curves"]={}
    plot["markers"]={}
    for item in elems.values():
        if "width" in item.keys() and "height" in item.keys():
            plot["xunit"] = float(item["width"])
            plot["yunit"] = float(item["height"])            
        elif "fill:none" in item["style"] and "d" in item.keys():
            plot["curves"][item["id"]] = parsecoords(item["d"])
        else:
            plot["markers"][item["id"]] = average(parsecoords(item["d"]))
    return plot

def load(infile):
    import xml.etree.ElementTree as ET
    tree = ET.parse(infile)
    root = tree.getroot()
    svg = root.find(prefix+"g")
    return svg
    
def main(args):
    svg = load(args.infile)
    
    elems = collect(svg)

    plot = identify(elems)

    coords = transform(plot)

    export(coords,args.outfile)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="convert an SVG into coordinates")
    parser.add_argument("infile",help="input file name")
    parser.add_argument("outfile",help="input file name")    
    main(parser.parse_args())
