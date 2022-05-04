// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2016-03-17
// Description : Load models from ROOT file and prepare them for fits

#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"

#if ROOT_VERSION_CODE < ROOT_VERSION(6,25,0)
#include "RooNameSet.h"
#else
#include "RooFitLegacy/RooNameSet.h"
#endif

#include "RooFitUtils/ExtendedModel.h"
#include "RooFitUtils/Log.h"
#include "RooFitUtils/Utils.h"

#include <chrono>

#include "TMath.h"
#include "TROOT.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TTime.h"
#include "TRandom.h"

#include "RooArgSet.h"
#include "RooProdPdf.h"

#include "TKey.h"

#include <iostream>
#include <sstream>

ClassImp(RooFitUtils::ExtendedModel)

// _____________________________________________________________________________

RooFitUtils::ExtendedModel::ExtendedModel(
					  const std::string &ModelName, const std::string &FileName,
					  const std::string &WsName, const std::string &ModelConfigName,
					  const std::string &DataName, const std::string &SnapshotName,
					  bool binnedLikelihood, RooArgSet *penalty, const std::string &TagAsMeasurement)
: TNamed(ModelName.c_str(), ModelName.c_str()), fFileName(FileName),
  fWsName(WsName), fModelConfigName(ModelConfigName), fDataName(DataName),
  fSnapshotName(SnapshotName), fBinnedLikelihood(binnedLikelihood),
  fTagAsMeasurement(TagAsMeasurement), fPenalty(penalty) {
  // Constructor
  initialise();
  
  coutP(InputArguments) << "ExtendedModel::ExtendedModel(" << fName
                        << ") created" << std::endl;
}

// _____________________________________________________________________________

RooFitUtils::ExtendedModel::~ExtendedModel() {
  // Destructor
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::initialise() {
  // Load all model information from specified file
  coutP(InputArguments) << "Opening file " << fFileName << std::endl;
  fFile = TFile::Open(fFileName.c_str());
  if (!fFile || !fFile->IsOpen()) {
    if (fFile) {
      delete fFile;
      fFile = NULL;
    }
    throw std::runtime_error(
        TString::Format("unable to open file '%s'", fFileName.c_str()).Data());
  }
  
  if(!fPenalty) {
    fPenalty = new RooArgSet();
  }

  TObject *ws = fFile->Get(fWsName.c_str());
  if (!ws) {
    std::stringstream ss;
    ss << "unable to load object '" << fWsName << "', available keys are \n";
    TList *keys = fFile->GetListOfKeys();
    TKey* thekey = 0;
    bool onlyOne = true;
    for (int i = 0; i < keys->GetEntries(); ++i) {
      TKey *k = dynamic_cast<TKey *>(keys->At(i));
      if (!k || strcmp(k->GetClassName(),"RooWorkspace")>0)
        continue;
      if(thekey) onlyOne=false;
      thekey = k;
      ss << "    '" << k->GetName() << "' (" << k->GetClassName() << ")\n";
    }
    if(thekey && onlyOne){
      thekey->Print();
      ws = thekey->ReadObj();
      coutW(InputArguments) << "unable to load object '" << fWsName << "', but found only one RooWorkspace in the file - using '" << ws->GetName() << "' instead!" << std::endl;
    }
    if(!thekey) ss << "    (none)";
    if(!ws){
      throw std::runtime_error(ss.str());
    }
  }
  fWorkspace = dynamic_cast<RooWorkspace *>(ws);
  if (!fWorkspace) {
    throw std::runtime_error(
        TString::Format("object '%s' is not a workspace", fWsName.c_str())
            .Data());
  }

  // Fixes for known features
  if (fBinnedLikelihood) {
    coutP(InputArguments) << "Activating binned likelihood evaluation"
                          << std::endl;
    RooFIter iter = fWorkspace->components().fwdIterator();
    RooAbsArg *arg;
    while ((arg = iter.next())) {
      if (arg->IsA() == RooRealSumPdf::Class()) {
        arg->setAttribute("BinnedLikelihood");
        coutI(InputArguments) << "Activating binned likelihood attribute for "
                              << arg->GetName() << std::endl;
      }
    }
  }

  if (fTagAsMeasurement != "") {
    coutP(InputArguments)
        << "Tagging CMS main measurements to reduce memory consumption"
        << std::endl;
    RooFIter iter = fWorkspace->components().fwdIterator();
    RooAbsArg *arg;
    while ((arg = iter.next())) {
      if (arg->IsA() == RooAddPdf::Class() &&
          TString(arg->GetName()).BeginsWith(fTagAsMeasurement.c_str())) {
        arg->setAttribute("MAIN_MEASUREMENT");
        coutI(InputArguments)
            << "Component " << arg->GetName() << " is a CMS main measurement";
      }
    }
  }

  // Continue loading the model
  coutP(InputArguments) << "Loading ModelConfig " << fModelConfigName
                        << std::endl;
  fModelConfig =
      (RooStats::ModelConfig *)(fWorkspace->obj(fModelConfigName.c_str()));
  if (!fModelConfig) {
    coutE(InputArguments)
        << "Something went wrong when loading the ModelConfig "
        << fModelConfigName << std::endl;
    throw std::runtime_error(TString::Format("unable to load ModelConfig '%s'",
                                             fModelConfigName.c_str())
                                 .Data());
  }

  coutP(InputArguments) << "Grabbing the pdf from the ModelConfig" << std::endl;
  fPdf = (RooAbsPdf *)fModelConfig->GetPdf();
  if (!fPdf) {
    coutE(InputArguments) << "Something went wrong when loading the pdf"
                          << std::endl;
    throw std::runtime_error("unable to obtain pdf");
  }

  fData = (RooAbsData *)(fWorkspace->data(fDataName.c_str()));
  if (!fData) {
    std::stringstream ss;
    ss << "unable to load dataset '" << fDataName
       << "', available datasets are \n";
    for (const auto &d : fWorkspace->allData()) {
      ss << "    '" << d->GetName() << "'\n";
    }
    throw std::runtime_error(ss.str());
  }

  coutP(InputArguments) << "Loading the nuisance parameters" << std::endl;
  fNuis = (RooArgSet *)fModelConfig->GetNuisanceParameters();
  if (!fNuis) {
    coutE(InputArguments)
        << "Something went wrong when loading the nuisance parameters"
        << std::endl;
  } else {
    fAllParams.add(*fNuis);
  }

  coutP(InputArguments) << "Loading the global observables" << std::endl;
  fGlobs = (RooArgSet *)fModelConfig->GetGlobalObservables();
  if (!fGlobs) {
    coutE(InputArguments)
        << "Something went wrong when loading the global observables"
        << std::endl;
  } else {
    fAllParams.add(*fGlobs);
  }

  coutP(InputArguments) << "Loading the parameters of interest" << std::endl;
  fPOIs = (RooArgSet *)fModelConfig->GetParametersOfInterest();
  if (!fPOIs) {
    coutE(InputArguments)
        << "Something went wrong when loading the parameters of interest"
        << std::endl;
    throw std::runtime_error("unable to obtain list of parameters of interest");
  }
	fAllParams.add(*fPOIs);

  coutP(InputArguments) << "Loading the observables" << std::endl;
  fObs = (RooArgSet *)fModelConfig->GetObservables();
  if (!fObs) {
    coutE(InputArguments) << "Something went wrong when loading the observables"
                          << std::endl;
    throw std::runtime_error("unable to obtain list of observables");
  }

  if (fSnapshotName != "") {
    coutP(InputArguments) << "Loading snapshots" << std::endl;
    std::vector<std::string> parsedSnapshots = parseString(fSnapshotName, ",");
    for (size_t i_snapshot = 0; i_snapshot < parsedSnapshots.size();
         ++i_snapshot) {
      std::string thisSnapshot = parsedSnapshots[i_snapshot];
      coutI(InputArguments) << "Loading snapshot " << thisSnapshot << std::endl;
      if(!fWorkspace->loadSnapshot(thisSnapshot.c_str())){
	std::stringstream ss;
	ss << "unable to load snapshot '" << thisSnapshot << "'";
#if ROOT_VERSION_CODE > ROOT_VERSION(6,25,0)	
	ss << ", available snapshots are are \n";
	for(const auto& snsh : fWorkspace->getSnapshots()){
	  ss << "  " << snsh->GetName() << "\n";
	}
#endif
	throw std::runtime_error(ss.str().c_str());
      }
    }
  }
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::floatParametersOfInterest(const std::string &floatName) {
  // Float a subset of the POIs at the specified values
  floatParametersOfInterest(parseString(floatName, ","));
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::floatParametersOfInterest(const std::vector<std::string> &parsed) {
  // Float a subset of the POIs at the specified values
  floatParameters(parsed,fModelConfig->GetParametersOfInterest());
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::floatNuisanceParameters(const std::string &floatName) {
  // Float a subset of the nuisance parameters at the specified values
  floatNuisanceParameters(parseString(floatName, ","));
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::floatNuisanceParameters(const std::vector<std::string> &parsed) {
  // Float a subset of the nuisance parameters at the specified values
  floatParameters(parsed,fModelConfig->GetNuisanceParameters());
}


// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::fixParametersOfInterest(const std::string &fixName) {
  // Fix a subset of the POIs at the specified values
  fixParametersOfInterest(parseString(fixName, ","));
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::fixParametersOfInterest(const std::vector<std::string> &parsed) {
  // Fix a subset of the POIs at the specified values
  fixParameters(parsed,fModelConfig->GetParametersOfInterest());
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::fixNuisanceParameters(const std::string &fixName) {
  // Fix a subset of the nuisance parameters at the specified values
  fixNuisanceParameters(parseString(fixName, ","));
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::fixNuisanceParameters(const std::vector<std::string> &parsed) {
  // Fix a subset of the nuisance parameters at the specified values
  floatParameters(parsed,fModelConfig->GetNuisanceParameters());
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::randomizeParametersOfInterest(const std::string &randomizeName) {
  // Randomize a subset of the POIs at the specified values
  randomizeParametersOfInterest(parseString(randomizeName, ","));
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::randomizeParametersOfInterest(const std::vector<std::string> &parsed) {
  // Randomize a subset of the POIs at the specified values
  randomizeParameters(parsed,fModelConfig->GetParametersOfInterest());
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::randomizeNuisanceParameters(const std::string &randomizeName) {
  // Randomize a subset of the nuisance parameters at the specified values
  randomizeNuisanceParameters(parseString(randomizeName, ","));
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::randomizeNuisanceParameters(const std::vector<std::string> &parsed) {
  // Randomize a subset of the nuisance parameters at the specified values
  randomizeParameters(parsed,fModelConfig->GetNuisanceParameters());
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::fixParameters(const std::string &fixName){
  // Fix a subset of the parameters at the specified values
  fixParameters(parseString(fixName,","));
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::fixParameters(const std::vector<std::string> &parsed) {
  // Fix a subset of the parameters at the specified values
  RooArgSet allVars(fWorkspace->allVars());
  fixParameters(parsed,&allVars);
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::floatParameters(const std::string &floatName){
  // Float a subset of the parameters
  floatParameters(parseString(floatName,","));
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::floatParameters(const std::vector<std::string> &parsed) {
  // Float a subset of the parameters at the specified values
  RooArgSet allVars(fWorkspace->allVars());
  floatParameters(parsed,&allVars);
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::randomizeParameters(const std::string &randomizeName){
  // Randomize a subset of the parameters
  randomizeParameters(parseString(randomizeName,","));
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::randomizeParameters(const std::vector<std::string> &parsed) {
  // Randomize a subset of the parameters at the specified values
  RooArgSet allVars(fWorkspace->allVars());
  randomizeParameters(parsed,&allVars);
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::fixParameters(const std::vector<std::string> &parsed, const RooArgSet* params) {
  // Fix a subset of the parameters at the specified values
  if(!params || params->getSize() == 0) return;
  for(const auto&s:parsed){
    auto nameEnd = s.find_first_of("=[");
    std::string pat(s.substr(0,nameEnd));
    int found = 0;
    double value = 0;
    bool valSet = false;
    if(s.size() != pat.size()){
      auto valEnd = pat.find_first_not_of("0123456789.",nameEnd);
      std::string val(s.substr(nameEnd+1,valEnd-nameEnd));
      valSet = true;
      value = atof(val.c_str());
    }
    for(auto obj:*params){
      if(RooFitUtils::matches(obj->GetName(),pat)){
        RooRealVar *par = dynamic_cast<RooRealVar *>(obj);
        if(par){
          found ++;
          if(valSet){
            par->setVal(value);
          }
          coutI(ObjectHandling) << "Fixing parameter " << par->GetName() << " at value " << par->getVal() << std::endl;
          par->setConstant(1);					
        }
      }
    }
    if (found == 0) {
      throw std::runtime_error("Parameter '" + pat + "' does not exist.");
    }
  }
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::floatParameters(const std::vector<std::string> &parsed, const RooArgSet* params) {
  // Fix a subset of the parameters at the specified values
  if(!params || params->getSize() == 0) return;  
  for(const auto&pat:parsed){
    int found = 0;
    for(auto obj:*params){
      if(RooFitUtils::matches(obj->GetName(),pat)){
        RooRealVar *par = dynamic_cast<RooRealVar *>(obj);
        if(par){
          found ++;
          coutI(ObjectHandling) << "Floating parameter " << par->GetName() << std::endl;
          par->setConstant(0);
        }
      }
    }
    if (found == 0) {
      throw std::runtime_error("Parameter '" + pat + "' does not exist.");
    }
  }
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::randomizeParameters(const std::vector<std::string> &parsed, const RooArgSet* params) {
  // Randomize a subset of the parameters at the specified values
  TRandom rnd;
  TObject* obj;
  for(const auto&pat:parsed){
    RooFIter itr(params->fwdIterator());
    int found = 0;
    while((obj = itr.next())){
      if(RooFitUtils::matches(obj->GetName(),pat)){
        RooRealVar *par = dynamic_cast<RooRealVar *>(obj);
        if(par){
          found ++;
          double err = par->getError();
          if(err < 1e-9) err = par->getVal();
          if(err < 1e-9) err = 1.;
          par->setVal(rnd.Gaus(par->getVal(),err));
          coutI(ObjectHandling) << "Randomizing parameter " << par->GetName() << ", new value is " << par->getVal() << std::endl;
        }
      }
    }
    if (found == 0) {
      coutE(ObjectHandling) << "Parameter " << pat
                            << " does not exist." << std::endl;
      exit(-1);
    }
  }
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::profileParameters(const std::string &profileName) {
  // Fix a subset of the nuisance parameters at the specified values
  this->profileParameters(parseString(profileName, ","));
}

// _____________________________________________________________________________

RooRealVar * RooFitUtils::ExtendedModel::parseParameter(const std::string &pname) {
  // parse a string parameter definition and return the parameter
  int sign = 0;bool useRange = false;double lo = 0;double hi = 0;bool useBoundary = false;double boundary = 0;
  return parseParameter(pname,sign,useRange,lo,hi,useBoundary,boundary);
}


// _____________________________________________________________________________

RooRealVar * RooFitUtils::ExtendedModel::parseParameter(const std::string &pname,int& sign,bool& useRange,double& lo,double& hi,bool& useBoundary,double& boundary) {
  // parse a string parameter definition and return the parameter
  TString thisName(pname);
  
  // Get ranges
  TString range;
  if (thisName.Contains("[")) {
    assert(thisName.Contains("]"));
    TObjArray *thisNameArray = thisName.Tokenize("[");
    thisName = ((TObjString *)thisNameArray->At(0))->GetString();
    range = ((TObjString *)thisNameArray->At(1))->GetString();
    range.ReplaceAll("]", "");
    assert(range.Contains(":"));
    TObjArray *rangeArray = range.Tokenize(":");
    rangeArray->SetOwner(true);
    TString s_lo = ((TObjString *)rangeArray->At(0))->GetString();
    TString s_hi = ((TObjString *)rangeArray->At(1))->GetString();
    delete rangeArray;
    lo = atof(s_lo.Data());
    hi = atof(s_hi.Data());
    useRange = kTRUE;
  }

  // Get sign
  if (thisName.Contains("+")) {
    thisName.ReplaceAll("+", ">0");
  } else if (thisName.Contains("-")) {
    thisName.ReplaceAll("-", "<0");
  }

  // Get boundaries
  if (thisName.Contains(">")) {
    TObjArray *thisNameArray = thisName.Tokenize(">");
    thisName = ((TObjString *)thisNameArray->At(0))->GetString();
    TString boundary_str = ((TObjString *)thisNameArray->At(1))->GetString();
    boundary = atof(boundary_str);
    sign = +1;
    useBoundary = true;
  } else if (thisName.Contains("<")) {
    TObjArray *thisNameArray = thisName.Tokenize("<");
    thisName = ((TObjString *)thisNameArray->At(0))->GetString();
    TString boundary_str = ((TObjString *)thisNameArray->At(1))->GetString();
    boundary = atof(boundary_str);
    sign = -1;
    useBoundary = true;
  } 
  
  RooRealVar *thisPoi = dynamic_cast<RooRealVar *>(this->fWorkspace->var(thisName));
  return thisPoi;
}


// _____________________________________________________________________________

RooRealVar * RooFitUtils::ExtendedModel::configureParameter(const std::string &pname) {
  // Fix a parameter at the specified value and/or apply ranges and boundaries
  TString thisName(pname.c_str());

  int sign = 0;bool useRange = false;double lo = 0;double hi = 0;bool useBoundary = false;double boundary = 0;
  RooRealVar* thisPoi = this->parseParameter(pname,sign,useRange,lo,hi,useBoundary,boundary);
  if(!thisPoi){
    coutE(ObjectHandling) << "Parameter " << thisName << " doesn't exist!"
                          << std::endl;
    return NULL;
  }
  double origVal = thisPoi->getVal();

  
  if(useRange){
    thisPoi->setRange(lo, hi);
    if ((origVal < lo) || (origVal > hi)) {
      double newVal = (hi - lo) / 2;
      thisPoi->setVal(newVal);
      coutI(ObjectHandling) << "Setting value to " << newVal << std::endl;
    }
  }
  
  if (useBoundary) {
    double forigVal = fabs(thisPoi->getVal());
    bool boundaryIsZero = AlmostEqualUlpsAndAbs(boundary, 0.0, 0.0001, 4);

    if (sign > 0) {
      thisPoi->setMin(boundary);
      if (origVal < boundary) {
        thisPoi->setVal(boundary);
      }
      if (boundaryIsZero && origVal < 0) {
        thisPoi->setVal(forigVal);
      }
    } else if (sign < 0) {
      thisPoi->setMax(boundary);
      if (origVal > boundary) {
        thisPoi->setVal(boundary);
      }
      if (boundaryIsZero && origVal > 0) {
        thisPoi->setVal(-forigVal);
      }
    }
  }

  thisPoi->setConstant(0);
//  coutI(ObjectHandling) << thisName.Data() << " = " << thisPoi->getVal()
//                        << " in [" << thisPoi->getMin() << ","
//                        << thisPoi->getMax() << "]" << std::endl;
  return thisPoi;
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::profileParameters(const std::vector<std::string> &parsed) {
  // Fix a subset of the nuisance parameters at the specified values
  for (const auto& parname:parsed){
    RooRealVar *thisPoi = this->configureParameter(parname);
    if (thisPoi) {
      coutI(ObjectHandling) << "Profiling parameter " << thisPoi->GetName() << std::endl;
    }
  }
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::setInitialErrors() {
  // Set initial errors of model parameters depending on constraint terms
  RooArgSet *AllConstraints = new RooArgSet();

  if (fWorkspace->set(Form("CACHE_CONSTR_OF_PDF_%s_FOR_OBS_%s", fPdf->GetName(),
                           RooNameSet(*fData->get()).content()))) {
    // Retrieve constraints from cache
    const RooArgSet *constr = fWorkspace->set(
        Form("CACHE_CONSTR_OF_PDF_%s_FOR_OBS_%s", fPdf->GetName(),
             RooNameSet(*fData->get()).content()));
    AllConstraints->add(*constr);
    delete constr;
  } else {
    // Load information needed to determine attributes from ModelConfig
    RooAbsPdf *tmpPdf = (RooAbsPdf *)fModelConfig->GetPdf();
    RooArgSet *tmpAllNuisanceParameters =
        (RooArgSet *)fModelConfig->GetNuisanceParameters();
    RooArgSet *tmpAllObservables = (RooArgSet *)fModelConfig->GetObservables();

    // Copies, to keep original sets in place after getAllconstraints call
    RooArgSet tmpAllNuisanceParameters2 = *tmpAllNuisanceParameters;
    RooArgSet tmpAllObservables2 = *tmpAllObservables;
    AllConstraints = tmpPdf->getAllConstraints(
        tmpAllObservables2, tmpAllNuisanceParameters2, kFALSE);
  }

  // Take care of the case where we have a product of constraint terms
  TIterator *ConstraintItrAll = AllConstraints->createIterator();
  RooAbsArg *nextConstraint;
  RooArgSet *tmpAllConstraints = new RooArgSet(AllConstraints->GetName());
  while ((nextConstraint = (RooAbsArg *)ConstraintItrAll->Next())) {
    if (nextConstraint->IsA() == RooProdPdf::Class()) {
      RooArgSet thisComponents;
      FindUniqueProdComponents((RooProdPdf *)nextConstraint, thisComponents);
      tmpAllConstraints->add(thisComponents);
    } else {
      coutI(ObjectHandling)
          << "Adding constraint " << nextConstraint->GetName() << std::endl;
      tmpAllConstraints->add(*nextConstraint);
    }
  }

  for (RooLinkedListIter it = fNuis->iterator();
       RooRealVar *nuip = dynamic_cast<RooRealVar *>(it.Next());) {
    coutI(ObjectHandling) << "On nuisance parameter " << nuip->GetName();
    double prefitvariation = 1.0;

    TIterator *ConstraintItr = tmpAllConstraints->createIterator();
    bool foundConstraint = kFALSE;
    bool foundGaussianConstraint = kFALSE;
    while ((nextConstraint = (RooAbsArg *)ConstraintItr->Next()) &&
           !foundConstraint) {
      if (nextConstraint->dependsOn(*nuip)) {
        foundConstraint = kTRUE;

        // Loop over global observables to match nuisance parameter and
        // global observable in case of a constrained nuisance parameter
        TIterator *GlobsItr = fGlobs->createIterator();
        RooRealVar *nextGlobalObservable;
        bool foundGlobalObservable = kFALSE;
        while ((nextGlobalObservable = (RooRealVar *)GlobsItr->Next()) &&
               !foundGlobalObservable) {
          if (nextConstraint->dependsOn(*nextGlobalObservable)) {
            foundGlobalObservable = kTRUE;

            // find constraint width in case of a Gaussian
            if (nextConstraint->IsA() == RooGaussian::Class()) {
              foundGaussianConstraint = kTRUE;
              double oldSigmaVal = 1.0;
              TIterator *ServerItr = nextConstraint->serverIterator();
              RooRealVar *nextServer;
              bool foundSigma = kFALSE;
              while ((nextServer = (RooRealVar *)ServerItr->Next()) &&
                     !foundSigma) {
                if (nextServer != nextGlobalObservable && nextServer != nuip) {
                  oldSigmaVal = nextServer->getVal();
                  foundSigma = kTRUE;
                }
              }

              if (AlmostEqualUlpsAndAbs(oldSigmaVal, 1.0, 0.001, 4)) {
                oldSigmaVal = 1.0;
              }

              if (!foundSigma) {
                coutI(ObjectHandling)
                    << "Sigma for pdf " << nextConstraint->GetName()
                    << " not found. Using 1.0." << std::endl;
              } else {
                coutI(ObjectHandling)
                    << "Using " << oldSigmaVal << " for sigma of pdf "
                    << nextConstraint->GetName() << std::endl;
              }

              prefitvariation = oldSigmaVal;
            }
          }
        }
        delete GlobsItr;
      }
    }
    delete ConstraintItr;

    if (foundGaussianConstraint) {
      coutP(ObjectHandling)
          << "Changing error of " << nuip->GetName() << " from "
          << nuip->getError() << " to " << prefitvariation << std::endl;
      nuip->setError(prefitvariation);
      nuip->removeRange();
    }
  }
}
