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
#include "RooFormulaVar.h"
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
                const std::string &DataName, const std::string &SnapshotName, const std::string &PdfName,
                bool binnedLikelihood = true, RooArgSet* penalty = NULL,
                const std::string &TagAsMeasurement = "pdf_");
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
  RooArgSet *GetPenalty() { return fPenalty && fPenalty->getSize() > 0 ? fPenalty : 0; }
  void addPenalty(RooAbsReal *x) {   
	fPenalty->add(*x);
  }

  
  void fixNuisanceParameters(const std::string &fixName = "*");
  void fixNuisanceParameters(const std::vector<std::string> &parsed);
  void fixParametersOfInterest(const std::string &fixName = "*");
  void fixParametersOfInterest(const std::vector<std::string> &parsed);
  void fixParameters(const std::string &fixName);	
  void fixParameters(const std::vector<std::string> &parsed);
	
  void floatNuisanceParameters(const std::string &floatName = "*");
  void floatNuisanceParameters(const std::vector<std::string> &parsed);
  void floatParametersOfInterest(const std::string &floatName = "*");
  void floatParametersOfInterest(const std::vector<std::string> &parsed);
  void floatParameters(const std::string &floatName);
  void floatParameters(const std::vector<std::string> &parsed);

  void randomizeNuisanceParameters(const std::string &floatName = "*");
  void randomizeNuisanceParameters(const std::vector<std::string> &parsed);
  void randomizeParametersOfInterest(const std::string &floatName = "*");
  void randomizeParametersOfInterest(const std::vector<std::string> &parsed);
  void randomizeParameters(const std::string &floatName);
  void randomizeParameters(const std::vector<std::string> &parsed);
	
  void profileParameters(const std::string &profileName);
  void profileParameters(const std::vector<std::string> &profileName);
  RooRealVar * parseParameter(const std::string &pname);
  RooArgSet configureParameter(const std::string &pname);
  void setInitialErrors();

  // _____________________________________________________________________________
protected:
  void initialise();
  RooRealVar * parseParameter(const std::string &pname,double& val,int& sign,bool& useRange,double& lo,double& hi,bool& useBoundary,double& boundary);
  void fixParameters(const std::vector<std::string> &parsed, const RooArgSet* params);
  void floatParameters(const std::vector<std::string> &parsed, const RooArgSet* params);
  void randomizeParameters(const std::vector<std::string> &parsed, const RooArgSet* params);	  
	
  // _____________________________________________________________________________
private:
  std::string fFileName;
  std::string fWsName;
  std::string fModelConfigName;
  std::string fPdfName;  
  std::string fDataName;
  std::string fSnapshotName;
  bool fBinnedLikelihood;
  std::string fTagAsMeasurement;

  TFile *fFile = 0;
  RooWorkspace *fWorkspace = 0;
  RooStats::ModelConfig *fModelConfig = 0;
  RooAbsPdf *fPdf = 0;
  RooAbsData *fData = 0;
  RooArgSet *fNuis = 0;
  RooArgSet *fGlobs = 0;
  RooArgSet *fPOIs = 0;
  RooArgSet *fObs = 0;
  RooArgSet *fPenalty = 0;
  RooArgSet fAllParams;
 

  // _____________________________________________________________________________
protected:
  ClassDef(ExtendedModel, 0)
};
}

#endif
