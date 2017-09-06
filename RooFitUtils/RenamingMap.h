// this file looks like plain C, but it's actually -*- c++ -*-
#ifndef RENAMINGMAP
#define RENAMINGMAP

#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "TNamed.h"
#include "TObjArray.h"
#include "TString.h"

#include "RooArgSet.h"
#include "RooMsgService.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"

#include "RooStats/ModelConfig.h"

namespace RooFitUtils {
class RenamingMap : public TNamed {

public:
  // Constructor and destructor
  RenamingMap(std::string MapName = "");
  ~RenamingMap();
  enum Attribute {
    Observable,
    ObservableRange,
    GlobalObservable,
    GlobalObservableRange,
    Sigma,
    SigmaRange,
    Constraint,
    Type
  };
  static const std::string AttributeNames[];
  enum MeasurementType { individual, combined };
  enum ConstraintType {
    automatic,
    unconstrained,
    Gaussian,
    Poisson,
    Lognormal,
    unknown
  };
  static const std::string ConstraintTypeNames[];

  // Accessors
  void SetAttribute(const std::string &ParameterName, Attribute thisAttribute,
                    const std::string &AttributeValue,
                    MeasurementType thisMeasurementType);
  std::string GetAttribute(const std::string &ParameterName,
                           Attribute thisAttribute,
                           MeasurementType thisMeasurementType);
  void AddAttributes(RooStats::ModelConfig *tmpModelConfig);

  std::map<std::string, std::string> GetRenamingMap() { return fRenamingMap; }
  void SetRenamingMap(std::map<std::string, std::string> RenamingMap) {
    fRenamingMap = RenamingMap;
  }

  // Steering
  void RenameParameter(const std::string &OldParameterName,
                       const std::string &NewParameterName);
  std::string FactoryExpression(const std::string &ParameterName,
                                MeasurementType);
  void FindUniqueProdComponents(RooProdPdf *Pdf, RooArgSet &Components);
  using TNamed::Print;
  void Print();

private:
  std::map<std::string /* OldParameterName */,
           std::string /* NewParameterName */>
      fRenamingMap;
  std::map<
      std::string /* OldParameterName */,
      std::map<Attribute /* obs, globs, pdf, ... */, std::string /* values */>>
      fOldParameterMap;
  std::map<
      std::string /* NewParameterName */,
      std::map<Attribute /* obs, globs, pdf, ... */, std::string /* values */>>
      fNewParameterMap;

protected:
  ClassDefOverride(RenamingMap, 1)
};
}

#endif
