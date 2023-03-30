// this file looks like plain C, but it's actually -*- c++ -*-
#ifndef COMBINEDMEASUREMENT
#define COMBINEDMEASUREMENT

#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "Math/MinimizerOptions.h"
#include "TFile.h"
#include "TIterator.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TString.h"

#include "RooAbsCollection.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooLinkedListIter.h"
#include "RooMinimizer.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"

#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/ModelConfig.h"

#include "RooFitUtils/AbsMeasurement.h"
#include "RooFitUtils/CorrelationScheme.h"
#include "RooFitUtils/Measurement.h"
#include "RooFitUtils/ParametrisationSequence.h"

namespace RooFitUtils {

class CombinedMeasurement : public AbsMeasurement {

  // ____________________________________________________________________________|__________
public:
  // Constructor and destructor
  CombinedMeasurement(const std::string &CombinedMeasurementName = "master",
                      const std::string &WorkspaceName = "combined",
                      const std::string &ModelConfigName = "ModelConfig",
                      const std::string &DataName = "combData");
  CombinedMeasurement(const std::string &CombinedMeasurementName,
                      const std::string &FileName,
                      const std::string &WorkspaceName,
                      const std::string &ModelConfigName,
                      const std::string &DataName);
  ~CombinedMeasurement();
  enum SnapshotName { background, nominal, ucmles };

  // Accessors
  void SetCorrelationScheme(CorrelationScheme *Scheme) {
    fCorrelationScheme = Scheme;
    fParametersOfInterestString = Scheme->GetParametersOfInterest();
  }
  CorrelationScheme *GetCorrelationScheme() { return fCorrelationScheme; }

  void SetParametrisationSequence(ParametrisationSequence *Sequence) {
    fParametrisationSequence = Sequence;
    fParametersOfInterestString = Sequence->GetParametersOfInterest();
  }
  ParametrisationSequence *GetParametrisationSequence() {
    return fParametrisationSequence;
  }

  void SetParametersOfInterest(std::string ParametersOfInterest) {
    fParametersOfInterestString = ParametersOfInterest;
  }
  std::string GetParametersOfInterest() { return fParametersOfInterestString; }

  // Steering
  virtual void initialise() override;
  void AddMeasurement(Measurement &aMeasurement) {
    fMeasurements[aMeasurement.GetName()] = (new Measurement(aMeasurement));
  }
  void AddMeasurement(Measurement *aMeasurement) {
    fMeasurements[aMeasurement->GetName()] = aMeasurement;
  }
  void CollectMeasurements();
  void CombineMeasurements();
  void ParametriseMeasurements();
  void DefineParametersOfInterest(const std::string &ParametersOfInterest,
                                  RooStats::ModelConfig *tmpModelConfig);
  void MakeCleanWorkspace();
  void MakeAsimovData(bool Conditional,
                      CombinedMeasurement::SnapshotName profileGenerateAt);
  void MakeAsimovData(bool Conditional,
                      CombinedMeasurement::SnapshotName profileAt,
                      CombinedMeasurement::SnapshotName generateAt);
  void MakeSnapshots(std::set<CombinedMeasurement::SnapshotName> Snapshots,
                     bool Conditional);
  void MakeSnapshots(CombinedMeasurement::SnapshotName Snapshot,
                     bool Conditional);
  using TNamed::Print;
  void Print();
  static void PrintCollection(RooAbsCollection *collection);

  // ____________________________________________________________________________|__________
protected:
  void DetermineAutoCorrelations(
				 std::map<std::string, RooArgSet *> &tmpAllNuisanceParameters, bool useConstraints = true);
  void UnfoldConstraints(RooArgSet &initial, RooArgSet &final, RooArgSet &obs,
                         RooArgSet &nuis, int &counter);

  // ____________________________________________________________________________|__________
private:
  std::string fParametersOfInterestString;

  CorrelationScheme *fCorrelationScheme;
  ParametrisationSequence *fParametrisationSequence;
  std::map<std::string, Measurement *> fMeasurements;

  // ____________________________________________________________________________|__________
protected:
  ClassDefOverride(CombinedMeasurement, 1)
};
}

#endif
