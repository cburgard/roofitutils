#!/usr/bin/env python
import sys
import itertools
import re

def iteritems(d):
    try:
        return d.iteritems()
    except AttributeError:
        return d.items()

def loadRooFitUtils():
    # retrieve the root core dir environment variable
    from ROOT import gSystem
    if gSystem.Load("libRooFitUtils"):
        raise ImportError("unable to load standalone libRooFitUtils.so!")

def loadMeasurements(combined,filename):
  try:
    from ConfigParser import ConfigParser
  except:
    from configparser import ConfigParser    
  config = ConfigParser()
  config.read(filename)
  for section in config.sections():
    if config.has_option(section,"Include"):
      include = config.get(section,"Include")
      config.read(include)
      config.remove_section(section)
  measurements = []
  for section in config.sections():
    print("adding "+section)
    measurement = ROOT.RooFitUtils.Measurement (section)
    measurement.SetBinnedLikelihood(config.getboolean(section,"BinnedLikelihood"))
    if config.get(section,"SnapshotName"):
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
  for target,sources in iteritems(corrmap):
    source = ",".join(sources)
    print("renaming '{0:s}' to '{1:s}'".format(source,target))
  return corrmap

def main(args):
  # Create a new combined measurement
  combined = ROOT.RooFitUtils.CombinedMeasurement("combined_master")
  
  measurements = loadMeasurements(combined,args.config)
  
  correlation = ROOT.RooFitUtils.CorrelationScheme("CorrelationScheme")
  if args.autocorr: correlation.SetAutoCorrelation(True)
  correlation.SetParametersOfInterest(",".join(args.poi))
  corrmap = loadCorrelations(correlation,args.correlations)
  # Run the combination. First all measurements are regularised, i.e. the
  # structure of the PDF will be unified, parameters and constraint terms renamed
  # according to a common convention, etc. Then a new simultaneous PDF and
  # dataset is build.
  for target,source in iteritems(corrmap):
    correlation.CorrelateParameter(",".join(source),target)
  combined.SetCorrelationScheme(correlation)
  combined.CollectMeasurements()
  combined.CombineMeasurements()
  
  # Generate Asimov data (NP measured in unconditional fit, generated for mu = 1)
  if args.asimov:
      combined.MakeAsimovData(ROOT.kTRUE, ROOT.RooFitUtils.CombinedMeasurement.ucmles, ROOT.RooFitUtils.CombinedMeasurement.nominal)
  
  workspace = combined.GetWorkspace()
  import re
  for p in workspace.allVars():
      for constpar in args.const:
          if re.match(constpar,p.GetName()):
              p.setConstant(True)
      
  workspace.writeToFile(args.output)
  
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
  parser.add_argument("--autocorr", help="Flag to use autocorrelation", default=False, action="store_true")  
  parser.add_argument("-c", "--correlations", nargs="+",default=[],
                      help="""Config file(s) with correlation definitions. This file contains lines like

  meas1::muGGF>>mu
  meas1::muVBF>>mu
  meas2::mu>>mu
to correlate e.g. the parameters 'muGGF' and 'muVBF' from the measurement 'meas1'
and the parameter 'mu' from the measurement 'meas2'. These files can be auto-generated 
with the 'guessCorrelations.py' script.""")
  parser.add_argument("--const", default=[], nargs="+", help="Parameters to set constant")
  parser.add_argument("--poi", default=[], nargs="+", help="Parameters of Interest in the new workspace")

  args = parser.parse_args()

  try:
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    loadRooFitUtils()
  except:
    print ("Could not load library. Make sure that it was compiled correctly.")

  main(args)
