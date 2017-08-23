// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2016-03-17
// Description : Load models from ROOT file and prepare them for fits

#ifndef EXTENDEDMODEL
#define EXTENDEDMODEL

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "TNamed.h"
#include "TFile.h"
#include "Math/MinimizerOptions.h"

#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooArgSet.h"
#include "RooWorkspace.h"
#include "RooMsgService.h"

#include "RooStats/ModelConfig.h"

#ifdef __MAKECINT__
#pragma link C++ class ExtendedModel;
#endif

class ExtendedModel : public TNamed {
  
// _____________________________________________________________________________
public:

  // Constructor and destructor
  ExtendedModel( const std::string& ModelName, const std::string& FileName, const std::string& WsName,
                 const std::string& ModelConfigName, const std::string& DataName, const std::string& SnapshotName,
                 bool binnedLikelihood = true, const std::string& TagAsMeasurement = "pdf_",
                 bool FixCache = true, bool FixMulti = true );
  virtual ~ExtendedModel();
  RooWorkspace* GetWorkspace() { return fWorkspace; }
  RooStats::ModelConfig* GetModelConfig() { return fModelConfig; }
  RooAbsPdf* GetPdf() { return fPdf; }
  RooAbsData* GetData() { return fData; }
  RooArgSet* GetNuisanceParameters() { return fNuis; }
  RooArgSet* GetGlobalObservables() { return fGlobs; }
  RooArgSet* GetParametersOfInterest() { return fPOIs; }
  RooArgSet* GetObservables() { return fObs; }

  void fixNuisanceParameters();
  void fixNuisanceParameters( const std::string& fixName );
  void fixNuisanceParameters( const std::vector<std::string>& parsed );
  void fixParametersOfInterest();
  void profileParameters( const std::string& profileName );
  void profileParameters( const std::vector<std::string>& profileName );
  RooRealVar* configureParameter(const std::string& pname);
  void setInitialErrors();

// _____________________________________________________________________________
protected:

  void initialise(bool fixCache, bool fixMulti);

// _____________________________________________________________________________
private:

  std::string fFileName;
  std::string fWsName;
  std::string fModelConfigName;
  std::string fDataName;
  std::string fSnapshotName;
  bool fBinnedLikelihood;
  std::string fTagAsMeasurement;

  TFile* fFile;
  RooWorkspace* fWorkspace;
  RooStats::ModelConfig* fModelConfig;
  RooAbsPdf* fPdf;
  RooAbsData* fData;
  RooArgSet* fNuis;
  RooArgSet* fGlobs;
  RooArgSet* fPOIs;
  RooArgSet* fObs;

// _____________________________________________________________________________
protected:

  ClassDef(ExtendedModel, 0)

};

#endif
