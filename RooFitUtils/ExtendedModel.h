// this file looks like plain C, but it's actually -*- c++ -*-
#ifndef EXTENDEDMODEL
#define EXTENDEDMODEL

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Math/MinimizerOptions.h"
#include "TFile.h"
#include "TNamed.h"

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooMsgService.h"
#include "RooWorkspace.h"

#include "RooStats/ModelConfig.h"

namespace RooFitUtils {

class ExtendedModel : public TNamed {

  // _____________________________________________________________________________
public:
  // Constructor and destructor
  ExtendedModel(const std::string &ModelName, const std::string &FileName,
                const std::string &WsName, const std::string &ModelConfigName,
                const std::string &DataName, const std::string &SnapshotName,
                bool binnedLikelihood = true,
                const std::string &TagAsMeasurement = "pdf_",
                bool FixCache = true, bool FixMulti = true);
  virtual ~ExtendedModel();
  RooWorkspace *GetWorkspace() { return fWorkspace; }
  RooStats::ModelConfig *GetModelConfig() { return fModelConfig; }
  RooAbsPdf *GetPdf() { return fPdf; }
  RooAbsData *GetData() { return fData; }
  RooArgSet *GetNuisanceParameters() { return fNuis; }
  RooArgSet *GetGlobalObservables() { return fGlobs; }
  RooArgSet *GetParametersOfInterest() { return fPOIs; }
  RooArgSet *GetParameterSet() { return &fAllParams; }
  RooArgSet *GetObservables() { return fObs; }

  void fixNuisanceParameters();
  void fixParametersOfInterest();
  void fixNuisanceParameters(const std::string &fixName);
  void fixNuisanceParameters(const std::vector<std::string> &parsed);
  void fixParametersOfInterest(const std::string &fixName);
  void fixParametersOfInterest(const std::vector<std::string> &parsed);
  void fixParameters(const std::vector<std::string> &parsed);  
  void profileParameters(const std::string &profileName);
  void profileParameters(const std::vector<std::string> &profileName);
  RooRealVar * parseParameter(const std::string &pname);
  RooRealVar *configureParameter(const std::string &pname);
  void setInitialErrors();

  // _____________________________________________________________________________
protected:
  void initialise(bool fixCache, bool fixMulti);
  RooRealVar * parseParameter(const std::string &pname,int& sign,bool& useRange,double& lo,double& hi,bool& useBoundary,double& boundary);
  
  // _____________________________________________________________________________
private:
  std::string fFileName;
  std::string fWsName;
  std::string fModelConfigName;
  std::string fDataName;
  std::string fSnapshotName;
  bool fBinnedLikelihood;
  std::string fTagAsMeasurement;

  TFile *fFile;
  RooWorkspace *fWorkspace;
  RooStats::ModelConfig *fModelConfig;
  RooAbsPdf *fPdf;
  RooAbsData *fData;
  RooArgSet *fNuis;
  RooArgSet *fGlobs;
  RooArgSet *fPOIs;
  RooArgSet *fObs;
	RooArgSet fAllParams;

  // _____________________________________________________________________________
protected:
  ClassDef(ExtendedModel, 0)
};
}

#endif
