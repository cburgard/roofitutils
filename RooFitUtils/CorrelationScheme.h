// this file looks like plain C, but it's actually -*- c++ -*-
#ifndef CORRELATIONSCHEME
#define CORRELATIONSCHEME

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

#include "TMatrixDSym.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TString.h"

#include "RooMsgService.h"

#include "RooFitUtils/RenamingMap.h"

namespace RooFitUtils {
  class CorrelationScheme : public TNamed {

  public:
    // Constructor and destructor
    CorrelationScheme(const std::string &SchemeName = "",
                      const std::string &ParametersOfInterest = "",
                      bool AutoCorrelation = kFALSE);
    ~CorrelationScheme();

    // Accessors
    void
    SetCorrelationMap(const std::map<std::string, RenamingMap> &CorrelationMap) {
      fCorrelationMap = CorrelationMap;
    }
    std::map<std::string, RenamingMap> GetCorrelationMap() {
      return fCorrelationMap;
    }

    void SetRenamingMap(const std::string &MeasurementName,
                        const RenamingMap &Map) {
      fCorrelationMap[MeasurementName] = Map;
    }
    RenamingMap GetRenamingMap(const std::string &MeasurementName) {
      return fCorrelationMap[MeasurementName];
    }

    void SetCorrelationFactors(
                               const std::string &MeasurementName,
                               const std::map<std::string, std::pair<TString, TMatrixDSym>>
                               &CorrelationFactors) {
      fCorrelationFactors[MeasurementName] = CorrelationFactors;
    }
    std::map<std::string, std::pair<TString, TMatrixDSym>>
      GetCorrelationFactors(const std::string &MeasurementName) {
      return fCorrelationFactors[MeasurementName];
    }

    void SetParametersOfInterest(std::string ParametersOfInterest) {
      fParametersOfInterest = ParametersOfInterest;
    }
    std::string GetParametersOfInterest() { return fParametersOfInterest; }

    void SetAutoCorrelation(bool setting);
    bool GetAutoCorrelation() { return fAutoCorrelation; }

    // Steering
    void CorrelateParameter(
                            const char *OldParameterNamePlusMeasurement, const char *NewParameterName,
                            RenamingMap::ConstraintType thisConstraintType = RenamingMap::automatic);
    void RenameParameter(
                         const char *MeasurementName, const char *OldParameterName,
                         const char *NewParameterName,
                         RenamingMap::ConstraintType thisConstraintType = RenamingMap::automatic);

    void CorrelateParameter(const char *OldParameterNamePlusMeasurement,
                            const char *NewParameterName, double rho);
    void CorrelateParameter(const char *OldParameterNamePlusMeasurement,
                            const char *NewParameterName, TMatrixDSym cov);
    void IntroduceCorrelation(const char *MeasurementName,
                              std::vector<TString> NewParameterNames,
                              TMatrixDSym cov);

    void ParseInputs(std::string &InputName, std::string &InputConstraintName,
                     std::string &InputObservableName,
                     std::string &InputObservableRange,
                     std::string &InputGlobalObservableName,
                     std::string &InputGlobalObservableRange,
                     std::string &InputSigmaName, std::string &InputSigmaRange);
    void DecomposeVariable(std::string &InputName, std::string &InputVariableName,
                           std::string &InputVariableRange);

    using TNamed::Print;
    void Print();
    void Print(const std::set<std::string>& thisMeasurements);
    void printToFile(const char* filename, bool standalone=false);
    void printToStream(std::ostream& out);
    void printToStream(std::ostream& out, const std::set<std::string>& thisMeasurements);

    // ____________________________________________________________________________|__________
  private:
    std::map<std::string, RenamingMap> fCorrelationMap;
    std::map<std::string, std::map<std::string, std::pair<TString, TMatrixDSym>>>
      fCorrelationFactors;
    std::string fParametersOfInterest;
    bool fAutoCorrelation;

    // ____________________________________________________________________________|__________
  protected:
    ClassDefOverride(CorrelationScheme, 1)
  };
}

#endif
