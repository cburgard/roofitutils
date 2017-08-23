#ifndef RENAMINGMAP
#define RENAMINGMAP

#include <string>
#include <map>
#include <iostream>
#include <sstream>

#include "TNamed.h"
#include "TString.h"
#include "TObjArray.h"

#include "RooMsgService.h"
#include "RooWorkspace.h"
#include "RooArgSet.h"
#include "RooSimultaneous.h"
#include "RooRealVar.h"
#include "RooProdPdf.h"

#include "RooStats/ModelConfig.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

class RenamingMap : public TNamed {

// ____________________________________________________________________________|__________
public:

  // Constructor and destructor
  RenamingMap( std::string MapName = "" );
  ~RenamingMap();
  enum Attribute { Observable, ObservableRange, GlobalObservable, GlobalObservableRange, Sigma, SigmaRange, Constraint, Type };
  static const std::string AttributeNames[];
  enum MeasurementType { individual, combined };
  enum ConstraintType { automatic, unconstrained, Gaussian, Poisson, Lognormal, unknown };
  static const std::string ConstraintTypeNames[];

  // Accessors
  void SetAttribute( std::string ParameterName, Attribute thisAttribute, std::string AttributeValue, MeasurementType thisMeasurementType );
  std::string GetAttribute( std::string ParameterName, Attribute thisAttribute, MeasurementType thisMeasurementType );
  void AddAttributes( ModelConfig* tmpModelConfig );

  std::map< std::string, std::string > GetRenamingMap() { return fRenamingMap; }
  void SetRenamingMap( std::map< std::string, std::string > RenamingMap) { fRenamingMap = RenamingMap; }

  // Steering
  void RenameParameter( std::string OldParameterName, std::string NewParameterName );
  std::string FactoryExpression( std::string ParameterName, MeasurementType  );
  void FindUniqueProdComponents( RooProdPdf* Pdf, RooArgSet& Components );
  using TNamed::Print;
  void Print();

// ____________________________________________________________________________|__________
private:

  std::map< std::string /* OldParameterName */, std::string /* NewParameterName */ > fRenamingMap;
  std::map< std::string /* OldParameterName */, std::map< Attribute /* obs, globs, pdf, ... */, std::string /* values */ > > fOldParameterMap;
  std::map< std::string /* NewParameterName */, std::map< Attribute /* obs, globs, pdf, ... */, std::string /* values */ > > fNewParameterMap;

// ____________________________________________________________________________|__________
protected:

  ClassDefOverride(RenamingMap, 1)

};

#endif
