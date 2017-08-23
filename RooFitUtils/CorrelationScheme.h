#ifndef CORRELATIONSCHEME
#define CORRELATIONSCHEME

#include <map>
#include <string>
#include <set>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "TNamed.h"
#include "TString.h"
#include "TObjString.h"
#include "TMatrixDSym.h"

#include "RooMsgService.h"

#include "RooFitUtils/RenamingMap.h"

using namespace std;
using namespace RooFit;

class CorrelationScheme : public TNamed {

// ____________________________________________________________________________|__________
public:

  // Constructor and destructor
  CorrelationScheme( std::string SchemeName = "", std::string ParametersOfInterest = "", bool AutoCorrelation = kFALSE );
  ~CorrelationScheme();

  // Accessors
  void SetCorrelationMap( const std::map< std::string, RenamingMap >& CorrelationMap ) { fCorrelationMap = CorrelationMap; }
  std::map< std::string, RenamingMap > GetCorrelationMap() { return fCorrelationMap; }

  void SetRenamingMap( const std::string& MeasurementName, const RenamingMap& Map ) { fCorrelationMap[MeasurementName] = Map; }
  RenamingMap GetRenamingMap( const std::string& MeasurementName ) { return fCorrelationMap[MeasurementName]; }

  void SetCorrelationFactors( const std::string& MeasurementName, const std::map< std::string, std::pair< TString, TMatrixDSym > >& CorrelationFactors ) { fCorrelationFactors[MeasurementName] = CorrelationFactors; }
  std::map< std::string, std::pair< TString, TMatrixDSym > > GetCorrelationFactors( const std::string& MeasurementName ) { return fCorrelationFactors[MeasurementName]; }

  void SetParametersOfInterest( std::string ParametersOfInterest ) { fParametersOfInterest = ParametersOfInterest; }
  std::string GetParametersOfInterest() { return fParametersOfInterest; }

  void SetAutoCorrelation( bool setting );
	bool GetAutoCorrelation() { return fAutoCorrelation; }

  // Steering
  void CorrelateParameter( const char* OldParameterNamePlusMeasurement, const char* NewParameterName, RenamingMap::ConstraintType thisConstraintType = RenamingMap::automatic );
  void RenameParameter( const char* MeasurementName, const char* OldParameterName, const char* NewParameterName, RenamingMap::ConstraintType thisConstraintType = RenamingMap::automatic );

  void CorrelateParameter( const char* OldParameterNamePlusMeasurement, const char* NewParameterName, double rho );
  void CorrelateParameter( const char* OldParameterNamePlusMeasurement, const char* NewParameterName, TMatrixDSym cov );
  void IntroduceCorrelation( const char* MeasurementName, std::vector<TString> NewParameterNames, TMatrixDSym cov );

  void ParseInputs( std::string& InputName, std::string& InputConstraintName, std::string& InputObservableName, std::string& InputObservableRange, std::string& InputGlobalObservableName, std::string& InputGlobalObservableRange, std::string& InputSigmaName, std::string& InputSigmaRange );
  void DecomposeVariable( std::string& InputName, std::string& InputVariableName, std::string& InputVariableRange );

  using TNamed::Print;
  void Print();
  void Print( std::set<std::string> thisMeasurements );
  static void PrintTable(string* firstCol, std::string** matrix, std::string** matrixErr, string* header, int nrRows, int nrCols, int nSigFig, ostream& ost, std::string indent = "", std::string delim = " & ", std::string ending = " \\\\");

// ____________________________________________________________________________|__________
private:

  std::map< std::string, RenamingMap > fCorrelationMap;
  std::map< std::string, std::map< std::string, std::pair< TString, TMatrixDSym > > > fCorrelationFactors;
  std::string fParametersOfInterest;
	bool fAutoCorrelation;

// ____________________________________________________________________________|__________
protected:

  ClassDefOverride(CorrelationScheme, 1)

};

#endif
