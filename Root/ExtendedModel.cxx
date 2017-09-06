// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2016-03-17
// Description : Load models from ROOT file and prepare them for fits

#include "RooRealSumPdf.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooRealVar.h"

#include "RooFitUtils/ExtendedModel.h"
#include "RooFitUtils/Log.h"
#include "RooFitUtils/Utils.h"

#include <chrono>

#include "TROOT.h"
#include "TTime.h"
#include "TSystem.h"
#include "TMath.h"

#include "RooProdPdf.h"
#include "RooArgSet.h"

#include "TKey.h"

#include <iostream>
#include <sstream>

using namespace std;
using namespace RooFit;
using namespace RooStats;

ClassImp(ExtendedModel)

// _____________________________________________________________________________

RooFitUtils::ExtendedModel::ExtendedModel( const std::string& ModelName, const std::string& FileName, const std::string& WsName, const std::string& ModelConfigName, const std::string& DataName, const std::string& SnapshotName, bool binnedLikelihood, const std::string& TagAsMeasurement, bool FixCache, bool FixMulti )
  :
  TNamed( ModelName.c_str(), ModelName.c_str() ),
  fFileName( FileName ),
  fWsName( WsName ),
  fModelConfigName( ModelConfigName ),
  fDataName( DataName ),
  fSnapshotName( SnapshotName ),
  fBinnedLikelihood( binnedLikelihood ),
  fTagAsMeasurement( TagAsMeasurement )
{
  // Constructor
  initialise(FixCache,FixMulti);

  coutP(InputArguments) << "ExtendedModel::ExtendedModel(" << fName <<") created" << endl;
}

// _____________________________________________________________________________

RooFitUtils::ExtendedModel::~ExtendedModel()
{
  // Destructor
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::initialise(bool fixCache, bool fixMulti)
{
  // Load all model information from specified file
  coutP(InputArguments) << "Opening file " << fFileName << endl;
  fFile = TFile::Open(fFileName.c_str());
  if(!fFile || !fFile->IsOpen()){
    if(fFile){
      delete fFile;
      fFile=NULL;
    }
    throw std::runtime_error(TString::Format("unable to open file '%s'",fFileName.c_str()).Data());
  }

  TObject* ws = fFile->Get(fWsName.c_str());
  if(!ws){
    std::stringstream ss;
    ss << "unable to load object '" << fWsName << "', available keys are \n";
    TList* keys = fFile->GetListOfKeys();
    for(int i=0; i<keys->GetEntries(); ++i){
      TKey* k = dynamic_cast<TKey*>( keys->At(i) );
      if(!k) continue;
      ss << "    '" <<  k->GetName() << "' (" << k->GetClassName() << ")\n";
    }
    throw std::runtime_error(ss.str());
  }
  fWorkspace = dynamic_cast<RooWorkspace*>(ws);
  if (!fWorkspace) {
    throw std::runtime_error(TString::Format("object '%s' is not a workspace",fWsName.c_str()).Data());
  }

  // Fixes for known features
  if (fBinnedLikelihood) {
    coutP(InputArguments) << "Activating binned likelihood evaluation" << endl;
    RooFIter iter = fWorkspace->components().fwdIterator();
    RooAbsArg* arg;
    while ((arg = iter.next())) {
      if (arg->IsA() == RooRealSumPdf::Class()) {
        arg->setAttribute("BinnedLikelihood");
        coutI(InputArguments) << "Activating binned likelihood attribute for " << arg->GetName() << endl;
      }
    }
  }

  if (fTagAsMeasurement != "") {
    coutP(InputArguments) << "Tagging CMS main measurements to reduce memory consumption" << endl;
    RooFIter iter = fWorkspace->components().fwdIterator() ;
    RooAbsArg* arg ;
    while ((arg = iter.next())) {
      if (arg->IsA()==RooAddPdf::Class() && TString(arg->GetName()).BeginsWith(fTagAsMeasurement.c_str())) {
      arg->setAttribute("MAIN_MEASUREMENT") ;
      coutI(InputArguments) << "Component " << arg->GetName() << " is a CMS main measurement";
      }
    }
  }

//   if (fixMulti) {
//     coutP(InputArguments) << "De-activating level 2 constant term optimization for RooMultiPdf" << endl;
//     RooFIter iter = fWorkspace->components().fwdIterator();
//     RooAbsArg* arg;
//     while ((arg = iter.next())) {
//       if (arg->IsA() == RooMultiPdf::Class()) {
//         arg->setAttribute("NOCacheAndTrack");
//         coutI(InputArguments) << "De-activation of level 2 constant term optimization for " << arg->GetName() << endl;
//       }
//     }
//   }

  if (kTRUE) {
    coutP(InputArguments) << "De-activating level 2 constant term optimization for specified pdfs" << endl;
    RooFIter iter = fWorkspace->components().fwdIterator();
    RooAbsArg* arg;
    int n = 0;
    while ((arg = iter.next())) {
      TString aname(arg->GetName());
      if (arg->InheritsFrom(RooAbsPdf::Class()) && (aname.EndsWith("_mm") || aname.Contains("mumu_atlas"))) {
	n++;
        arg->setAttribute("NOCacheAndTrack");
        coutI(InputArguments) << "De-activation of level 2 constant term optimization for " << arg->GetName();
      }
    }
  }

  // Continue loading the model
  coutP(InputArguments) << "Loading ModelConfig " << fModelConfigName << endl;
  fModelConfig = (ModelConfig*)(fWorkspace->obj(fModelConfigName.c_str()));
  if (!fModelConfig) {
    coutE(InputArguments) << "Something went wrong when loading the ModelConfig " << fModelConfigName << endl;
    throw std::runtime_error(TString::Format("unable to load ModelConfig '%s'",fModelConfigName.c_str()).Data());
  }

  coutP(InputArguments) << "Grabbing the pdf from the ModelConfig" << endl;
  fPdf = (RooAbsPdf*)fModelConfig->GetPdf();
  if (!fPdf) {
    coutE(InputArguments) << "Something went wrong when loading the pdf" << endl;
    throw std::runtime_error("unable to obtain pdf");
  }

  fData = (RooAbsData*)(fWorkspace->data(fDataName.c_str()));
  if (!fData) {
    std::stringstream ss;
    ss << "unable to load dataset '" << fDataName << "', available datasets are \n";
    for(const auto& d:fWorkspace->allData()){
      ss << "    '" <<  d->GetName() << "'\n";
    }
    throw std::runtime_error(ss.str());
  }

  coutP(InputArguments) << "Loading the nuisance parameters" << endl;
  fNuis = (RooArgSet*)fModelConfig->GetNuisanceParameters();
  if (!fNuis) {
    coutE(InputArguments) << "Something went wrong when loading the nuisance parameters" << endl;
    throw std::runtime_error("unable to obtain list of nuisance parameters");
  }

  coutP(InputArguments) << "Loading the global observables" << endl;
  fGlobs = (RooArgSet*)fModelConfig->GetGlobalObservables();
  if (!fGlobs) {
    coutE(InputArguments) << "Something went wrong when loading the global observables" << endl;
    throw std::runtime_error("unable to obtain list of global observables");
  }

  coutP(InputArguments) << "Loading the parameters of interest" << endl;
  fPOIs = (RooArgSet*)fModelConfig->GetParametersOfInterest();
  if (!fPOIs) {
    coutE(InputArguments) << "Something went wrong when loading the parameters of interest" << endl;
    throw std::runtime_error("unable to obtain list of parameters of interest");
  }

  coutP(InputArguments) << "Loading the observables" << endl;
  fObs = (RooArgSet*)fModelConfig->GetObservables();
  if (!fObs) {
    coutE(InputArguments) << "Something went wrong when loading the observables" << endl;
    throw std::runtime_error("unable to obtain list of observables");
  }

  if (fSnapshotName != "") {
    coutP(InputArguments) << "Loading snapshots" << endl;
    vector<string> parsedSnapshots = parseString(fSnapshotName, ",");
    for (size_t i_snapshot = 0; i_snapshot < parsedSnapshots.size(); ++i_snapshot) {
      string thisSnapshot = parsedSnapshots[i_snapshot];
      coutI(InputArguments) << "Loading snapshot " << thisSnapshot << endl;
      fWorkspace->loadSnapshot(thisSnapshot.c_str());
    }
  }
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::fixNuisanceParameters()
{
  // Fix all nuisance parameters
  for (RooLinkedListIter it = fNuis->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
    Double_t value = v->getVal();
    string name = v->GetName();
    coutI(ObjectHandling) << "Fixing nuisance parameter " << name << " at value " << value << endl;
    v->setConstant(1);
  }
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::fixParametersOfInterest()
{
  // Fix all parameters of interest
  for (RooLinkedListIter it = fPOIs->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
    Double_t value = v->getVal();
    string name = v->GetName();
    coutI(ObjectHandling) << "Fixing parameter of interest " << name << " at value " << value << endl;
    v->setConstant(1);
  }
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::fixNuisanceParameters( const std::string& fixName )
{
  // Fix a subset of the nuisance parameters at the specified values
  fixNuisanceParameters(parseString(fixName, ","));
}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::fixNuisanceParameters( const std::vector<std::string>& parsed )
{
  // Fix a subset of the nuisance parameters at the specified values  
  for (size_t i = 0; i < parsed.size(); i++) {
     TString thisName = parsed[i].c_str();
     TString thisVal;
     if (thisName.Contains("[")) {
       assert(thisName.Contains("]"));
       TObjArray* thisNameArray = thisName.Tokenize("[");
       thisName = ((TObjString*)thisNameArray->At(0))->GetString();
       thisVal = ((TObjString*)thisNameArray->At(1))->GetString();
       thisVal.ReplaceAll("]","");
     }

     RooRealVar* par = (RooRealVar*)fWorkspace->var(thisName.Data());
     if (!par) {
       coutE(ObjectHandling) << "Nuisance parameter " << thisName.Data() << " does not exist." << endl;
       exit(-1);
     }

     double value = par->getVal();
     if (thisVal.IsFloat()) {
       value = thisVal.Atof();
       par->setVal(value);
     }

     coutI(ObjectHandling) << "Fixing nuisance parameter " << thisName.Data() << " at value " << value << endl;
     par->setConstant(1);
   }

}

// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::profileParameters( const std::string& profileName )
{
  // Fix a subset of the nuisance parameters at the specified values
  profileParameters(parseString(profileName, ","));
}

// _____________________________________________________________________________

RooRealVar* RooFitUtils::ExtendedModel::configureParameter(const std::string& pname){
  // Fix a subset of the nuisance parameters at the specified values
  TString thisName(pname.c_str());
    TString range;
    TString boundary;
    int sign = 0;

    bool useRange    = kFALSE;
    bool useBoundary = kFALSE;

    // Get ranges
    if (thisName.Contains("[")) {
      assert(thisName.Contains("]"));
      TObjArray* thisNameArray = thisName.Tokenize("[");
      thisName = ((TObjString*)thisNameArray->At(0))->GetString();
      range = ((TObjString*)thisNameArray->At(1))->GetString();
      range.ReplaceAll("]","");
      assert(range.Contains(":"));
      useRange = kTRUE;
    }

    // Get sign
    if (thisName.Contains("+")) {
      thisName.ReplaceAll("+",">0");
    } else if (thisName.Contains("-")) {
      thisName.ReplaceAll("-","<0");
    }

    // Get boundaries
    if (thisName.Contains(">")) {
      TObjArray* thisNameArray = thisName.Tokenize(">");
      thisName = ((TObjString*)thisNameArray->At(0))->GetString();
      boundary = ((TObjString*)thisNameArray->At(1))->GetString();
      sign = +1;
      useBoundary = kTRUE;
    } else if (thisName.Contains("<")) {
      TObjArray* thisNameArray = thisName.Tokenize("<");
      thisName = ((TObjString*)thisNameArray->At(0))->GetString();
      boundary = ((TObjString*)thisNameArray->At(1))->GetString();
      sign = -1;
      useBoundary = kTRUE;
    }

    RooRealVar* thisPoi = (RooRealVar*)fWorkspace->var(thisName);
    if (!thisPoi) {
      coutE(ObjectHandling) << "Parameter " << thisName << " doesn't exist!" << endl;
      return NULL;
    }

    if (useRange) {
      double origVal = thisPoi->getVal();
      TObjArray* rangeArray = range.Tokenize(":");
      TString s_lo = ((TObjString*)rangeArray->At(0))->GetString();
      TString s_hi = ((TObjString*)rangeArray->At(1))->GetString();
      double lo = atof(s_lo.Data());
      double hi = atof(s_hi.Data());
      thisPoi->setRange(lo, hi);
      if ((origVal < lo) || (origVal > hi)) {
        double newVal = (hi - lo) / 2;
        thisPoi->setVal(newVal);
        coutI(ObjectHandling) << "Setting value to " << newVal << endl;
      }
    }

    if (useBoundary) {
      double tmpBoundary = atof(boundary.Data());
      double origVal = thisPoi->getVal();
      double forigVal = fabs(thisPoi->getVal());
      bool boundaryIsZero = AlmostEqualUlpsAndAbs(tmpBoundary, 0.0, 0.0001, 4);

      if (sign > 0) {
        thisPoi->setMin(tmpBoundary);
        if (origVal < tmpBoundary) {
          thisPoi->setVal(tmpBoundary);
        }
        if (boundaryIsZero && origVal < 0) {
          thisPoi->setVal(forigVal);
        }
      } else if (sign < 0) {
        thisPoi->setMax(tmpBoundary);
        if (origVal > tmpBoundary) {
          thisPoi->setVal(tmpBoundary);
        }
        if (boundaryIsZero && origVal > 0) {
          thisPoi->setVal(-forigVal);
        }
      }
    }

    thisPoi->setConstant(0);
    coutI(ObjectHandling) << thisName.Data() << " = " << thisPoi->getVal() << " in [" << thisPoi->getMin() << "," << thisPoi->getMax() << "]" << endl;
    return thisPoi;

}


// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::profileParameters( const std::vector<std::string>& parsed )
{
  // Fix a subset of the nuisance parameters at the specified values
  for (size_t i = 0; i < parsed.size(); i++) {
    RooRealVar* thisPoi = this->configureParameter(parsed[i]);
    if(thisPoi){
      coutI(ObjectHandling) << "Profiling parameter " << thisPoi->GetName() << endl;
    }
  }
}


// _____________________________________________________________________________

void RooFitUtils::ExtendedModel::setInitialErrors()
{
  // Set initial errors of model parameters depending on constraint terms
  RooArgSet* AllConstraints = new RooArgSet();

  if (fWorkspace->set(Form("CACHE_CONSTR_OF_PDF_%s_FOR_OBS_%s", fPdf->GetName(), RooNameSet(*fData->get()).content()))) {
    // Retrieve constraints from cache
    const RooArgSet* constr = fWorkspace->set(Form("CACHE_CONSTR_OF_PDF_%s_FOR_OBS_%s", fPdf->GetName(), RooNameSet(*fData->get()).content()));
    AllConstraints->add(*constr);
    delete constr;
  } else {
    // Load information needed to determine attributes from ModelConfig
    RooAbsPdf* tmpPdf = (RooAbsPdf*)fModelConfig->GetPdf();
    RooArgSet* tmpAllNuisanceParameters = (RooArgSet*)fModelConfig->GetNuisanceParameters();
    RooArgSet* tmpAllObservables = (RooArgSet*)fModelConfig->GetObservables();

    // Copies, to keep original sets in place after getAllconstraints call
    RooArgSet tmpAllNuisanceParameters2 = *tmpAllNuisanceParameters;
    RooArgSet tmpAllObservables2 = *tmpAllObservables;
    AllConstraints = tmpPdf->getAllConstraints(tmpAllObservables2, tmpAllNuisanceParameters2, kFALSE);
  }

  // Take care of the case where we have a product of constraint terms
  TIterator* ConstraintItrAll = AllConstraints->createIterator();
  RooAbsArg* nextConstraint;
  RooArgSet* tmpAllConstraints = new RooArgSet(AllConstraints->GetName());
  while ((nextConstraint = (RooAbsArg*)ConstraintItrAll->Next())) {
    if (nextConstraint->IsA() == RooProdPdf::Class()) {
      RooArgSet thisComponents;
      FindUniqueProdComponents((RooProdPdf*)nextConstraint, thisComponents);
      tmpAllConstraints->add(thisComponents);
    } else {
      coutI(ObjectHandling) << "Adding constraint " << nextConstraint->GetName() << endl;
      tmpAllConstraints->add(*nextConstraint);
    }
  }

  for (RooLinkedListIter it = fNuis->iterator(); RooRealVar* nuip = dynamic_cast<RooRealVar*>(it.Next());) {
    coutI(ObjectHandling) << "On nuisance parameter " << nuip->GetName();
    double prefitvariation = 1.0;

    TIterator* ConstraintItr = tmpAllConstraints->createIterator();
    bool foundConstraint = kFALSE;
    bool foundGaussianConstraint = kFALSE;
    while ((nextConstraint = (RooAbsArg*)ConstraintItr->Next()) && !foundConstraint) {
      if (nextConstraint->dependsOn(*nuip)) {
        foundConstraint = kTRUE;

        // Loop over global observables to match nuisance parameter and
        // global observable in case of a constrained nuisance parameter
        TIterator* GlobsItr = fGlobs->createIterator();
        RooRealVar* nextGlobalObservable;
        bool foundGlobalObservable = kFALSE;
        while ((nextGlobalObservable = (RooRealVar*)GlobsItr->Next()) && !foundGlobalObservable) {
          if (nextConstraint->dependsOn(*nextGlobalObservable)) {
            foundGlobalObservable = kTRUE;

            // find constraint width in case of a Gaussian
            if (nextConstraint->IsA() == RooGaussian::Class()) {
              foundGaussianConstraint = kTRUE;
              double oldSigmaVal = 1.0;
              TIterator* ServerItr = nextConstraint->serverIterator();
              RooRealVar* nextServer;
              bool foundSigma = kFALSE;
              while ((nextServer = (RooRealVar*)ServerItr->Next()) && !foundSigma) {
                if (nextServer != nextGlobalObservable && nextServer != nuip) {
                  oldSigmaVal = nextServer->getVal();
                  foundSigma = kTRUE;
                }
              }

              if (AlmostEqualUlpsAndAbs(oldSigmaVal, 1.0, 0.001, 4)) {
                oldSigmaVal = 1.0;
              }

              if (!foundSigma) {
                coutI(ObjectHandling) << "Sigma for pdf " << nextConstraint->GetName() << " not found. Using 1.0." << endl;
              } else {
                coutI(ObjectHandling)  << "Using " << oldSigmaVal << " for sigma of pdf " << nextConstraint->GetName() << endl;
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
      coutP(ObjectHandling) << "Changing error of " << nuip->GetName() << " from " << nuip->getError() << " to " << prefitvariation << endl;
      nuip->setError(prefitvariation);
      nuip->removeRange();
    }
  }
}
