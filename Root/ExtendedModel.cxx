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
#elif ROOT_VERSION_CODE < ROOT_VERSION(6,27,0)
#include "RooFitLegacy/RooNameSet.h"
#else
#define HAS_NO_RooNameSet
#endif

#if ROOT_VERSION_CODE > ROOT_VERSION(6,26,0)
#include "RooFitHS3/RooJSONFactoryWSTool.h"
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

namespace {
  bool hasEnding (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
      return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
      return false;
    }
  }
}

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,30,0)
namespace RooHelpers {
  std::string getColonSeparatedNameString(RooArgSet const &argSet, char delim = ':')
  {
    
    RooArgList tmp(argSet);
    tmp.sort();
    
    std::string content;
    for (auto const &arg : tmp) {
      content += arg->GetName();
      content += delim;
    }
    if (!content.empty()) {
      content.pop_back();
    }
    return content;
  }
}
#endif
  
// _____________________________________________________________________________

RooFitUtils::ExtendedModel::ExtendedModel(
					  const std::string &ModelName, const std::string &FileName,
					  const std::string &WsName, const std::string &ModelConfigName,
					  const std::string &DataName, const std::string &SnapshotName, const std::string& PdfName,
					  bool binnedLikelihood, RooArgSet *penalty, const std::string &TagAsMeasurement)
: TNamed(ModelName.c_str(), ModelName.c_str()), fFileName(FileName),
  fWsName(WsName), fModelConfigName(ModelConfigName), fPdfName(PdfName), fDataName(DataName),
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

  if(::hasEnding(fFileName,".root")){
    fFile = TFile::Open(fFileName.c_str());
    if (!fFile || !fFile->IsOpen()) {
      if (fFile) {
	delete fFile;
	fFile = NULL;
      }
      throw std::runtime_error(
			       TString::Format("unable to open file '%s'", fFileName.c_str()).Data());
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
	ws = thekey->ReadObj();
	coutW(InputArguments) << "unable to load object '" << fWsName << "', but found only one RooWorkspace in the file - using '" << ws->GetName() << "' instead!" << std::endl;
      }
      if(!thekey) ss << "    (none)";
      if(!ws){
	throw std::runtime_error(ss.str());
      }
    }
    fWorkspace = dynamic_cast<RooWorkspace *>(ws);
    
  } else if(::hasEnding(fFileName,".json")){
#if ROOT_VERSION_CODE > ROOT_VERSION(6,26,0)
    fWorkspace = new RooWorkspace("workspace");
    RooJSONFactoryWSTool tool(*fWorkspace);
    tool.importJSON(fFileName);
#else
    throw std::runtime_error("unable to read JSON files in this ROOT version!");
#endif
  }
  
  if(!fPenalty) {
    fPenalty = new RooArgSet();
  }

  if (!fWorkspace) {
    throw std::runtime_error(
        TString::Format("object '%s' is not a workspace", fWsName.c_str())
            .Data());
  }

  // Fixes for known features
  coutP(InputArguments) << (fBinnedLikelihood ? "Activating" : "Deactivating") << " binned likelihood evaluation"  << std::endl;
  {
    for(auto* arg : fWorkspace->components()){
      if (arg->IsA() == RooRealSumPdf::Class()) {
	arg->setAttribute("BinnedLikelihood",fBinnedLikelihood);
	coutP(InputArguments) << (fBinnedLikelihood ? "Activating" : "Deactivating") << " binned likelihood attribute for " << arg->GetName() << std::endl;
      }
    }
  }

  if (fTagAsMeasurement != "") {
    coutP(InputArguments)
        << "Tagging CMS main measurements to reduce memory consumption"
        << std::endl;
    for(auto* arg : fWorkspace->components()){
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
  }

  if(fModelConfig){
    coutP(InputArguments) << "Grabbing the pdf from the ModelConfig" << std::endl;
    fPdf = (RooAbsPdf *)fModelConfig->GetPdf();
  } else {
    fPdf = fWorkspace->pdf(fPdfName.c_str());
  }
  if (!fPdf) {
    coutE(InputArguments) << "Something went wrong when loading the pdf - none given in Model Config, and none named '" << fPdfName << "' available"
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

  if(fModelConfig){
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
  if(fModelConfig) floatParameters(parsed,fModelConfig->GetParametersOfInterest());
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::floatNuisanceParameters(const std::string &floatName) {
  // Float a subset of the nuisance parameters at the specified values
  if(fModelConfig) floatNuisanceParameters(parseString(floatName, ","));
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::floatNuisanceParameters(const std::vector<std::string> &parsed) {
  // Float a subset of the nuisance parameters at the specified values
  if(fModelConfig) floatParameters(parsed,fModelConfig->GetNuisanceParameters());
}


// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::fixParametersOfInterest(const std::string &fixName) {
  // Fix a subset of the POIs at the specified values
  fixParametersOfInterest(parseString(fixName, ","));
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::fixParametersOfInterest(const std::vector<std::string> &parsed) {
  // Fix a subset of the POIs at the specified values
  if(fModelConfig) fixParameters(parsed,fModelConfig->GetParametersOfInterest());
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::fixNuisanceParameters(const std::string &fixName) {
  // Fix a subset of the nuisance parameters at the specified values
  fixNuisanceParameters(parseString(fixName, ","));
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::fixNuisanceParameters(const std::vector<std::string> &parsed) {
  // Fix a subset of the nuisance parameters at the specified values
  if(fModelConfig){
    fixParameters(parsed,fModelConfig->GetNuisanceParameters());
  } else {
    RooArgSet allVars(fWorkspace->allVars());
    fixParameters(parsed,&allVars);

  }
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::randomizeParametersOfInterest(const std::string &randomizeName) {
  // Randomize a subset of the POIs at the specified values
  randomizeParametersOfInterest(parseString(randomizeName, ","));
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::randomizeParametersOfInterest(const std::vector<std::string> &parsed) {
  // Randomize a subset of the POIs at the specified values
  if(fModelConfig) randomizeParameters(parsed,fModelConfig->GetParametersOfInterest());
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::randomizeNuisanceParameters(const std::string &randomizeName) {
  // Randomize a subset of the nuisance parameters at the specified values
  randomizeNuisanceParameters(parseString(randomizeName, ","));
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::randomizeNuisanceParameters(const std::vector<std::string> &parsed) {
  // Randomize a subset of the nuisance parameters at the specified values
  if(fModelConfig) randomizeParameters(parsed,fModelConfig->GetNuisanceParameters());
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
  for(const auto&pat:parsed){
    int found = 0;
    for(auto* obj: *params){
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
  int sign = 0;bool useRange = false;double lo = 0;double hi = 0;bool useBoundary = false;double boundary = 0;double val = 0;
  return parseParameter(pname,val,sign,useRange,lo,hi,useBoundary,boundary);
}


// _____________________________________________________________________________

RooRealVar * RooFitUtils::ExtendedModel::parseParameter(const std::string &pname,double& val,int& sign,bool& useRange,double& lo,double& hi,bool& useBoundary,double& boundary) {
  // parse a string parameter definition and return the parameter
  TString thisName(pname);

  // Prepare
  thisName.ReplaceAll("+", ">0");
  thisName.ReplaceAll("-", "<0");

  // Get starting values
  bool has_val = false;
  if (thisName.Contains("=")){
    int pos = thisName.First("=");
    TString val_str = thisName(pos+1,thisName.Length()-pos);
    val = atof(val_str);
    thisName.Remove(pos);
    has_val = true;
  }
      
  // Get ranges
  TString range;
  if (thisName.Contains("[")) {
    assert(thisName.Contains("]"));
    int begin = thisName.First("[");
    int mid = thisName.First(":");    
    int end = thisName.First("]");
    TString s_lo = thisName(begin,mid-begin);
    TString s_hi = thisName(begin,end-mid);    
    lo = atof(s_lo.Data());
    hi = atof(s_hi.Data());
    useRange = kTRUE;
    thisName.Remove(begin,end-begin);
  }

  // Get boundaries
  if (thisName.Contains(">")) {
    int pos = thisName.First(">");
    TString boundary_str = thisName(pos+1,thisName.Length()-pos);
    boundary = atof(boundary_str);
    sign = +1;
    useBoundary = true;
    thisName.Remove(pos);
  } else if (thisName.Contains("<")) {
    int pos = thisName.First("<");
    TString boundary_str = thisName(pos+1,thisName.Length()-pos);
    boundary = atof(boundary_str);
    sign = -1;
    useBoundary = true;
    thisName.Remove(pos);
  } 
  
  RooRealVar *thisPoi = dynamic_cast<RooRealVar *>(this->fWorkspace->var(thisName));
  if(!has_val) val = thisPoi->getVal();
  return thisPoi;
}


// _____________________________________________________________________________

RooArgSet RooFitUtils::ExtendedModel::configureParameter(const std::string &pname) {
  // Fix a parameter at the specified value and/or apply ranges and boundaries
  TString thisName(pname.c_str());

  RooArgSet retval;
  
  int sign = 0;double origVal=0;bool useRange = false;double lo = 0;double hi = 0;bool useBoundary = false;double boundary = 0;
  if(thisName.Contains("*")){
    for(auto* var:this->fWorkspace->allVars()){
      if(RooFitUtils::matches(var->GetName(),pname)){
	static_cast<RooRealVar*>(var)->setConstant(0);
	retval.add(*var);
      }
    }
    return retval;
  }
  
  RooRealVar* thisPoi = this->parseParameter(pname,origVal,sign,useRange,lo,hi,useBoundary,boundary);
  if(!thisPoi){
    coutE(ObjectHandling) << "Parameter " << thisName << " doesn't exist!"
                          << std::endl;
    return retval;
  }
  retval.add(*thisPoi);
  
  if(useRange){
    thisPoi->setRange(lo, hi);
    if ((origVal < lo) || (origVal > hi)) {
      double newVal = (hi - lo) / 2;
      thisPoi->setVal(newVal);
      coutI(ObjectHandling) << "Setting value to " << newVal << std::endl;
    }
  }
  
  thisPoi->setVal(origVal);  
  
  if (useBoundary) {
    bool boundaryIsZero = AlmostEqualUlpsAndAbs(boundary, 0.0, 0.0001, 4);

    if (sign > 0) {
      thisPoi->setMin(boundary);
      if (origVal < boundary) {
        thisPoi->setVal(boundary);
      }
      if (boundaryIsZero && origVal < 0) {
        thisPoi->setVal(origVal);
      }
    } else if (sign < 0) {
      thisPoi->setMax(boundary);
      if (origVal > boundary) {
        thisPoi->setVal(boundary);
      }
      if (boundaryIsZero && origVal > 0) {
        thisPoi->setVal(-origVal);
      }
    }
  }

  thisPoi->setConstant(0);
//  coutI(ObjectHandling) << thisName.Data() << " = " << thisPoi->getVal()
//                        << " in [" << thisPoi->getMin() << ","
//                        << thisPoi->getMax() << "]" << std::endl;

  return retval;
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::profileParameters(const std::vector<std::string> &parsed) {
  // Fix a subset of the nuisance parameters at the specified values
  for (const auto& parname:parsed){
    auto poi = this->configureParameter(parname);
    if (poi.getSize()>0) {
      coutI(ObjectHandling) << "Profiling parameter " << poi.first()->GetName() << std::endl;
    }
  }
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::setInitialErrors() {
  // Set initial errors of model parameters depending on constraint terms
  if(fPOIs){
    for(auto* it : *fPOIs){
      RooRealVar *poi = dynamic_cast<RooRealVar *>(it);
      poi->setError(0.* (poi->getMax() - poi->getMin()));
    }
  }

  if(fNuis){
    RooArgSet *AllConstraints = new RooArgSet();

#ifndef HAS_NO_RooNameSet
    if (fWorkspace->set(Form("CACHE_CONSTR_OF_PDF_%s_FOR_OBS_%s", fPdf->GetName(),
			     RooNameSet(*fData->get()).content()))) {
      // Retrieve constraints from cache
      const RooArgSet *constr = fWorkspace->set(
						Form("CACHE_CONSTR_OF_PDF_%s_FOR_OBS_%s", fPdf->GetName(),
						     RooNameSet(*fData->get()).content()));
#else
      if (fWorkspace->set(Form("CACHE_CONSTR_OF_PDF_%s_FOR_OBS_%s", fPdf->GetName(),
			       RooHelpers::getColonSeparatedNameString(*fData->get()).c_str()))) {
	// Retrieve constraints from cache
	const RooArgSet *constr = fWorkspace->set(
						  Form("CACHE_CONSTR_OF_PDF_%s_FOR_OBS_%s", fPdf->GetName(),
						       RooHelpers::getColonSeparatedNameString(*fData->get()).c_str()));
#endif    
	AllConstraints->add(*constr);
	delete constr;
      } else {
	// Load information needed to determine attributes from ModelConfig
	if(fModelConfig){
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
      }

      // Take care of the case where we have a product of constraint terms
      RooArgSet *tmpAllConstraints = new RooArgSet(AllConstraints->GetName());
      for(auto* nextConstraint : *AllConstraints){
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

      for(auto* it : *fNuis){
	RooRealVar *nuip = dynamic_cast<RooRealVar *>(it);
	coutI(ObjectHandling) << "On nuisance parameter " << nuip->GetName();
	double prefitvariation = 1.0;

	bool foundConstraint = kFALSE;
	bool foundGaussianConstraint = kFALSE;
	for(auto* next : *tmpAllConstraints){
	  RooAbsReal* nextConstraint = static_cast<RooAbsReal*>(next);
	  if(foundConstraint) break;
	  if (nextConstraint->dependsOn(*nuip)) {
	    foundConstraint = kTRUE;

	    // Loop over global observables to match nuisance parameter and
	    // global observable in case of a constrained nuisance parameter
	    for(auto* next : *fGlobs){
	      RooRealVar *nextGlobalObservable = static_cast<RooRealVar*>(next);
	      if (nextConstraint->dependsOn(*nextGlobalObservable)) {

		// find constraint width in case of a Gaussian
		if (nextConstraint->IsA() == RooGaussian::Class()) {
		  foundGaussianConstraint = kTRUE;
		  double oldSigmaVal = 1.0;
		  bool foundSigma = kFALSE;
		  for(auto* next : servers(nextConstraint)){
		    RooRealVar* nextServer = dynamic_cast<RooRealVar*>(next);
		    if(foundSigma) break;
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
	  }
	}

	if (foundGaussianConstraint) {
	  coutP(ObjectHandling)
	    << "Changing error of " << nuip->GetName() << " from "
	    << nuip->getError() << " to " << prefitvariation << std::endl;
	  nuip->setError(prefitvariation);
	  nuip->removeRange();
	}
      }
    }
  }
