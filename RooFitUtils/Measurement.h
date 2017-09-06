//this file looks like plain C, but it's actually -*- c++ -*-
#ifndef MEASUREMENT
#define MEASUREMENT

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "TNamed.h"
#include "TFile.h"
#include "TList.h"
#include "TIterator.h"

#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooArgSet.h"
#include "RooWorkspace.h"
#include "RooAbsArg.h"
#include "RooLinkedListIter.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooProdPdf.h"
#include "RooArgList.h"
#include "RooMsgService.h"
#include "RooRealSumPdf.h"

#include "RooStats/ModelConfig.h"

#include "RooFitUtils/AbsMeasurement.h"
#include "RooFitUtils/Channel.h"
#include "RooFitUtils/RenamingMap.h"

namespace RooFitUtils {

class Measurement : public AbsMeasurement {

// ____________________________________________________________________________|__________
public:

  // Constructor and destructor
  Measurement( const std::string& MeasurementName = "", const std::string& FileName = "", const std::string& WorkspaceName = "", const std::string& ModelConfigName = "", const std::string& DataName = "", const std::string& SnapshotName = "nominalNuis", bool BinnedLikelihood = kTRUE );
  Measurement( const std::string& MeasurementName, RooWorkspace* ws, const std::string& ModelConfigName = "", const std::string& DataName = "", const std::string& SnapshotName = "nominalNuis", bool BinnedLikelihood = kTRUE);
  ~Measurement();

  // Accessors
  void SetRenamingMap( const RenamingMap& Map ) { fRenamingMap = Map; }
  RenamingMap GetRenamingMap() { return fRenamingMap; }

  void SetCorrelationFactors( const std::map< std::string, std::pair< TString, TMatrixDSym > >& Factors ) { fCorrelationFactors = Factors; }
  std::map< std::string, std::pair< TString, TMatrixDSym > > GetCorrelationFactors() { return fCorrelationFactors; }

  void SetChannels( std::vector< Channel* > Channels ) { fChannels = Channels; }
  std::vector< Channel* > GetChannels() { return fChannels; }

  void SetBinnedLikelihood( bool BinnedLikelihood ) { fBinnedLikelihood = BinnedLikelihood; }
  bool GetBinnedLikelihood() { return fBinnedLikelihood; }

  void SetSnapshotName( const std::string& SnapshotName ) { fSnapshotName = SnapshotName; }
  std::string GetSnapshotName() { return fSnapshotName; }

  void SetChannelFilter( const std::string& Filter ) { fChannelFilter = Filter; }
  std::string GetChannelFilter() { return fChannelFilter; }

  // Steering
  virtual void initialise() override;
  void CollectChannels();
  RenamingMap::ConstraintType DetermineConstraintType( RooRealVar* Parameter );
  using TNamed::Print;
  void Print();

// ____________________________________________________________________________|__________
protected:

  void FindUniqueProdComponents( RooProdPdf* Pdf, RooArgSet& Components );
  void TagPrunableParameters();
  void hist2dataset( RooSimultaneous* thisSim, RooCategory* cat );

// ____________________________________________________________________________|__________
private:

  std::string fSnapshotName;

  TList* dataList = NULL;
  RooArgSet* fConstraints = NULL;
  bool fBinnedLikelihood = true;

  std::vector< Channel* > fChannels;
  std::vector< std::string > fExcludedChannels;
  std::string fChannelFilter = ".";
  RenamingMap fRenamingMap;
  std::map< std::string, std::pair< TString, TMatrixDSym > > fCorrelationFactors;

// ____________________________________________________________________________|__________
protected:

  ClassDefOverride(Measurement, 1)

};

}

#endif
