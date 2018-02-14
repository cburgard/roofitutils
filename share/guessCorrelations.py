#!/bin/env python
try:
  from fuzzywuzzy import fuzz
  from fuzzywuzzy import process
except ImportError:
  print("this script requires the 'fuzzywuzzy' python package")
  print("to install this locally, execute the following commands:")
  print("  wget https://bootstrap.pypa.io/get-pip.py")
  print("  python get-pip.py --user && sed -i '1s|.*|#!/usr/bin/env python|' $HOME/.local/bin/pip")
  print("  $HOME/.local/bin/pip install python-Levenshtein --user")
  print("  $HOME/.local/bin/pip install fuzzywuzzy --user")
  exit(1)

import sys
if sys.version_info >= (3, 0):
  def iteritems(d):
    return d.items()
else:
  def iteritems(d):
    return d.iteritems()  
  
# argparse is for argument parsing
import argparse
# ROOT is needed to read the workspace
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
# re is needed to allow blacklisting parameters
import re
# json is needed to read the json/python file with the experimental parameters
import json
# xml.etree is needed to read the xml file with the theory parameters
from xml.etree import ElementTree

import itertools  

def convertCamelCaseToSnakeCase(name):
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()

nulltokens = []
  
def preprocess(instring):
  s = convertCamelCaseToSnakeCase(instring)
  for tok in nulltokens:
    s = s.replace(tok,"").strip("_")
  return sorted(re.split('_+',s))

def bestmatch(a,b,cutoff,func):
  best = 0
  xf,yf = None,None
  for x in a:
    for y in b:
      s = func(x,y)
      if s > best:
        best = s
        xf,yf = x,y
  if best < cutoff:
    return None,None,0
  return xf,yf,best

def comparesets(a,b,cutoff,func):
  ntok = max(len(a),len(b))
  if min(len(a),len(b)) < 1: return 0
  mya = [ x for x in a ]
  myb = [ y for y in b ]
  totalscore = 0.
  while len(mya)>0 and len(myb)>0:
    ea,eb,score = bestmatch(mya,myb,cutoff,func)
    if score < cutoff:
      return totalscore/ntok
    totalscore = totalscore + score
    mya.remove(ea)
    myb.remove(eb)
  return totalscore/ntok

def matchsets(a,b,cutoff):
  ntok = min(len(a),len(b))
  matches = []
  if ntok < 1: return (matches,0)
  mya = [ x for x in a ]
  myb = [ y for y in b ]
  totalscore = 0.
  while len(mya)>0 and len(myb)>0:
    ea,eb,score = bestmatch(mya,myb,cutoff,lambda x, y: comparesets(preprocess(x),preprocess(y),cutoff,fuzz.ratio))
    if score < cutoff:
      return (matches,totalscore)
    totalscore = totalscore + score
    mya.remove(ea)
    myb.remove(eb)
    matches.append((ea,eb,score))
  return (matches,totalscore)

def mkitr(orig):
  # convert anything into a list
  l = list()
  if not orig:
    return l
  try:
    itr = orig.fwdIterator()
    obj = itr.next()
    while obj:
      l.append(obj)
      obj = itr.next()
    return l
  except AttributeError:
    for x in orig:
      l.append(x)
      return l

def blacklisted(name,blacklist):
  # check if a string is blacklisted
  for i in blacklist:
    if i.match(name):
      return True
  return False


def getWorkspace(wsfilename):
  # open the file and get the workspace
  try:
    filename,wsname = wsfilename.split(":")
    if not filename or not wsname: raise ValueError()
  except ValueError as e:
    raise RuntimeError("workspace needs to be given in the format 'filename:workspacename', got '{0:s}' instead!".format(wsfilename))
  global infile
  infile = ROOT.TFile.Open(filename)
  if not infile or not infile.IsOpen():
    raise RuntimeError("unable to open file '{0:s}'!".format(filename))
  ws = infile.Get(wsname)
  if not ws:
    raise RuntimeError("unable to find workspace named '{0:s}' in file '{1:s}'!".format(wsname,filename))
  if not ws.InheritsFrom(ROOT.RooWorkspace.Class()):
    raise RuntimeError("object '{0:s}' in file '{1:s}' is not a workspace!".format(wsname,filename))
  return ws
  

def collectParameters(inlist,blacklist):
  params = []
  suppressed = []
  for var in mkitr(inlist):
    name = str(var.GetName())
    if blacklisted(name,blacklist):
      suppressed.append(name)
      continue
    params.append(name)
    pass
  return params,suppressed  

def collectParametersFromWorkspace(ws,blacklist):
  # crawl the workspace for all parameters
  return collectParameters(ws.allVars(),blacklist)

def collectParametersFromModelConfig(mc,blacklist):
  # crawl the workspace for all parameters
  np,npsupp = collectParameters(mc.GetParametersOfInterest(),blacklist)
  poi,poisupp = collectParameters(mc.GetNuisanceParameters(),blacklist)
  return np+poi,npsupp+poisupp

def printMatches(triplets):
  for orig,new,score in triplets:
    print("\t{:s} = {:s} (@{:.2f}%)".format(orig,new,score))

def printMissing(triplets,allpars):
  params = [ p for p in allpars ]
  for name,match,score in triplets:
    try:
      params.remove(name)
    except ValueError:
      pass
  for p in params:
    print("\t"+p)

def checkExperimentalParameters(params,infilenames,threshold=80):
  # identify experimental parameters based on the python/json file(s)
  print("Attempting to identify experimental NP Schemes...")
  results = []
  for pyjs in infilenames:
    with open(pyjs,"r") as f:
      for groupname,npgroup in iteritems(json.load(f)):
        highscore = 0.
        bestset = None
        mappings = []
        for setname,npset in iteritems(npgroup):
          matches,score = matchsets(npset,params,threshold)
          if score > highscore:
            mappings = matches
            highscore = score
            bestset = setname
        if bestset:
          bestgroup = npgroup[bestset]
          print("Identified '{:s}': '{:s}' (matched {:.2f}/{:d}):".format(groupname,bestset,float(highscore/100),int(len(bestgroup))))
          printMatches(mappings)
          print("  missing parameters are")
          printMissing(mappings,npgroup[bestset])
          results = results + mappings
          for new,old,score in mappings:
            params.remove(old)
            pass
          pass
        else:
          print("Unable to identify '{0:s}', tried '{1:s}'".format(groupname,"','".join(npgroup.keys())))
  return results

def checkTheoryParameters(params,infilenames,threshold=80):
  # identify theory parameters based on the XML file(s)
  print("\nAttempting to identify theory NPs")
  theoparams = []
  for xmlfilename in infilenames:
    with open(xmlfilename, 'rt') as xmlfile:
      xml = ElementTree.parse(xmlfile)
      for param in xml.iter('Param'):
        theoparams.append(str(param.attrib["Name"]))
  matches,totalscore = matchsets(theoparams,params,threshold)
  notfound = []
  for name,match,score in matches:
    params.remove(match)
    theoparams.remove(name)
  print("The following parameters were found and identified in the workspace: ({:d} total) ".format(len(matches)))
  printMatches(matches)
  print("The following parameters were not found in the workspace: ({:d} total)".format(len(theoparams)))
  printMissing(matches,theoparams)
  return matches
          
def main(args):
  # main steering function
  global nulltokens
  nulltokens = args.nulltokens
  
  # get the workspace
  ws = getWorkspace(args.workspace)

  # collect all the parameters
  blacklist = [ re.compile(i) for i in args.blacklist]
  if args.quick:
    mc = ws.obj("ModelConfig")
    params,suppressed = collectParametersFromModelConfig(mc,blacklist)
  else:
    params,suppressed = collectParametersFromWorkspace(ws,blacklist)
  mappings = [ ]

  # check the experimental ones
  mappings = mappings + checkExperimentalParameters(params,args.json,args.pyThresholds)

  # check the theory ones
  mappings = mappings + checkTheoryParameters(params,args.xml,args.xmlThresholds)

  # print everything to a correlation file
  if args.write:
    with open(args.write[0],"w") as outfile:
      for new,old,score in mappings:
        if old == new: continue
        outfile.write("{:s}::{:s}>>{:s}\n".format(args.write[1],old,new))

  # print a summary
  print("\nSummary:")
  print("The following parameters in the workspace could not be identified: {:s} ({:d} total)".format(",".join(sorted(params)),len(params)))
  print("The following expressions were blacklisted: {:s} ({:d} objects suppressed)".format(",".join(args.blacklist),len(suppressed)))

  #exit
  return True
  
if __name__ == "__main__":
  # create an argument parser
  parser = argparse.ArgumentParser("Workspace Checking Tool")
  parser.add_argument("workspace",metavar="file.root:workspacename",default=[],help="workspace to be checked")
  parser.add_argument("--xml",dest="xml",default=[],nargs="*",help="list of xml files to be used")
  parser.add_argument("--xmlThresholds",dest="xmlThresholds",type=int,help="threshold to be used for matching with XMLs",default=95)
  parser.add_argument("--py",dest="json",default=[],nargs="*",help="list of python/json files to be used")
  parser.add_argument("--pyThresholds",dest="pyThresholds",type=int,help="threshold to be used for matching with python/json",default=80)
  parser.add_argument("--blacklist",dest="blacklist",nargs="*",default=["obs_.*","binWidth_obs_.*",".*gamma_stat_.*","nom_.*",".*sampleNorm.*"],required=False,help="list of expressions to be ignored")
  parser.add_argument("--write",dest="write",default=None,metavar="FILE.TXT NAME",nargs=2,type=str,help="output file/name to write the correlations to")
  parser.add_argument('--quick', action='store_true',help="run in quick mode and only use the parameters present in the ModelConfig")
  parser.add_argument("--ignore",dest="nulltokens",nargs="*",default=["ATLAS","alpha"],required=False,help="list of name components to be ignored")
  # call the main function
  main(parser.parse_args())
