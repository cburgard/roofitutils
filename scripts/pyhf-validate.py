#!/bin/env python

def main(args):
    import pyhf
    with open(args.input) as infile:
        ws = pyhf.Workspace(infile)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="validate a JSON file to see if it's valid pyhf JSON")
    parser.add_argument("input",help="the JSON file")
    args = parser.parse_args()
    main(args)


