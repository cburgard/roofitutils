//this file looks like plain C, but it's actually -*- c++ -*-
#ifndef CHANNEL
#define CHANNEL

#include <string>
#include <map>
#include <iostream>
#include <sstream>

#include "TNamed.h"
#include "TString.h"
#include "TMatrixDSym.h"
#include "TIterator.h"

#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsArg.h"
#include "RooMsgService.h"

#include "RooFitUtils/RenamingMap.h"

namespace RooFitUtils {
class Channel : public TNamed {

public:

  // Constructor and destructor
  Channel( std::string ChannelName, RooAbsPdf* Pdf, RooAbsData* Data, std::string ParentName );
  ~Channel();

  // Accessors
  void SetPdf( RooAbsPdf* Pdf ) { fPdf = Pdf; }
  RooAbsPdf* GetPdf() { return fPdf; }

  void SetData( RooAbsData* Data ) { fData = Data; }
  RooAbsData* GetData() { return fData; }

  void SetParentName( const std::string& ParentName ) { fParentMeasurement = ParentName; }
  std::string GetParentName() { return fParentMeasurement; }

  void SetRenamingMap( RenamingMap& Map ) { fRenamingMap = Map; }
  RenamingMap GetRenamingMap() { return fRenamingMap; }

  void SetCorrelationFactors( std::map< std::string, std::pair< TString, TMatrixDSym > >& Factors ) { fCorrelationFactors = Factors; }
  std::map< std::string, std::pair< TString, TMatrixDSym > > GetCorrelationFactors() { return fCorrelationFactors; }

  void SetNuisanceParameters( RooArgSet* NuisanceParameters ) { fNuisanceParameters = NuisanceParameters; }
  RooArgSet* GetNuisanceParameters() { return fNuisanceParameters; }

  void SetObservables( RooArgSet* Observables ) { fObservables = Observables; }
  RooArgSet* GetObservables() { return fObservables; }

  void SetGlobalObservables( RooArgSet* GlobalObservables ) { fGlobalObservables = GlobalObservables; }
  RooArgSet* GetGlobalObservables() { return fGlobalObservables; }

  // Steering
  void RegulariseChannel();

private:

  RooWorkspace* fWorkSpace;
  RooAbsPdf* fPdf;
  RooAbsData* fData;
  std::string fParentMeasurement;
  RenamingMap fRenamingMap;
  std::map< std::string, std::pair< TString, TMatrixDSym > > fCorrelationFactors;
  RooArgSet* fNuisanceParameters;
  RooArgSet* fObservables;
  RooArgSet* fGlobalObservables;

protected:

  ClassDefOverride(Channel, 1)

};

}

#endif
