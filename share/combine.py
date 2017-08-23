#!/usr/bin/env python
import sys
import itertools
import re

import ConfigParser

def loadMeasurements(combined,filename):
  config = ConfigParser.ConfigParser()
  config.read(filename)
  for section in config.sections():
    if config.has_option(section,"Include"):
      include = config.get(section,"Include")
      config.read(include)
      config.remove_section(section)
  measurements = []
  for section in config.sections():
    print("adding "+section)
    measurement = ROOT.Measurement (section)
    measurement.SetBinnedLikelihood(config.getboolean(section,"BinnedLikelihood"))
    measurement.SetSnapshotName(config.get(section,"SnapshotName"))
    measurement.SetFileName(config.get(section,"FileName"))
    measurement.SetWorkspaceName(config.get(section,"WorkspaceName"))
    measurement.SetModelConfigName(config.get(section,"ModelConfigName"))
    measurement.SetDataName(config.get(section,"DataName"))
    combined.AddMeasurement(measurement)
    measurements.append(measurement)
  return measurements

def loadCorrelations(correlation,schemes):
  corrmap = {}
  for corrschemefile in schemes:
    with open(corrschemefile,'r') as corrscheme:
      for line in corrscheme:
        s = line.strip()
        if not s or s[0] == '#': continue
        source,target = s.split(">>")
        sources = source.split(",")
        if target in corrmap.keys():
          corrmap[target] += sources
        else:
          corrmap[target] = sources
  for target,sources in corrmap.iteritems():
    source = ",".join(sources)
    print("renaming '{0:s}' to '{1:s}'".format(source,target))
  return corrmap

def main(args):
  # Create a new combined measurement
  combined = ROOT.CombinedMeasurement("combined_master")
  
  measurements = loadMeasurements(combined,args.config)
  
  correlation = ROOT.CorrelationScheme("CorrelationScheme")
  correlation.SetParametersOfInterest("mu")
  corrmap = loadCorrelations(correlation,args.correlations)
  # Run the combination. First all measurements are regularised, i.e. the
  # structure of the PDF will be unified, parameters and constraint terms renamed
  # according to a common convention, etc. Then a new simultaneous PDF and
  # dataset is build.
  for target,source in corrmap.iteritems():
    correlation.CorrelateParameter(",".join(source),target)
  combined.SetCorrelationScheme(correlation)
  combined.CollectMeasurements()
  combined.CombineMeasurements()
  
  # Generate Asimov data (NP measured in unconditional fit, generated for mu = 1)
  if args.asimov:
      combined.MakeAsimovData(ROOT.kTRUE, ROOT.CombinedMeasurement.ucmles, ROOT.CombinedMeasurement.nominal)
  
  combined.writeToFile(args.output)
  
  # Print useful information like the correlation scheme, re-namings, etc.
  combined.Print()


# Read options from the command line, e.g. output filename, model, ...
if __name__ == "__main__":
  import argparse

  configexample = ""


  parser = argparse.ArgumentParser("NIKHEF Statistics Suite workspace combination script",formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("-i", "--input", nargs="+", dest="config", required=True,
                      help="""Config file(s) with channels that enter the combination. 
Each channel corresponds to a block like this:
  [BlockName]                         # e.g. 'HGAM'
  Measurement: measurementName        # e.g. 'hgam'
  BinnedLikelihood: [True|False]
  SnapshotName:  snapshotName         # e.g. 'nominalNuis'
  FileName: /path/to/workspace.root
  WorkspaceName: wsName               # e.g. 'combWS'
  ModelConfigName: ModelConfig
  DataName: datasName                 # e.g. 'combData' or 'AsimovSB'""")
  parser.add_argument("-o", "--output", default="combWS.root", help="Desired path to the output file")
  parser.add_argument("-a", "--asimov", help="Flag to generate Asimov data for the combined model",action="store_true")
  parser.add_argument("-c", "--correlations", nargs="+",default=["correlations.txt"],
                      help="""Config file(s) with correlation definitions. This file contains lines like
  meas1::muGGF>>mu
  meas1::muVBF>>mu
  meas2::mu>>mu
to correlate e.g. the parameters 'muGGF' and 'muVBF' from the measurement 'meas1'
and the parameter 'mu' from the measurement 'meas2'. These files can be auto-generated 
with the 'guessCorrelations.py' script.""")

  args = parser.parse_args()

  try:
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.ProcessLine(".x $ROOTCOREDIR/scripts/load_packages.C")
  except:
    print ("Could not load library. Make sure that it was compiled correctly.")

  main(args)
