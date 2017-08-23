#ifndef PARAMETRISATIONSCHEME
#define PARAMETRISATIONSCHEME

#include <string>
#include <list>
#include <iostream>
#include <sstream>

#include "TNamed.h"

#include "RooMsgService.h"

using namespace std;
using namespace RooFit;

class ParametrisationScheme : public TNamed {

// ____________________________________________________________________________|__________
public:

  // Constructor and destructor
  ParametrisationScheme( std::string SchemeName = "" );
  ~ParametrisationScheme();

  // Accessors
  void SetExpressions( const std::list< std::string >& Expressions ) { fExpressions = Expressions; }
  std::list< std::string > GetExpressions() { return fExpressions; }

  void SetNewNuisanceParameters( const std::list< std::string >& NewNuisanceParameters ) { fNewNuisanceParameters = NewNuisanceParameters; }
  std::list< std::string > GetNewNuisanceParameters() { return fNewNuisanceParameters; }

  void SetNewGlobalObservables( const std::list< std::string >& NewGlobalObservables ) { fNewGlobalObservables = NewGlobalObservables; }
  std::list< std::string > GetNewGlobalObservables() { return fNewGlobalObservables; }

  // Steering
  void AddExpression( std::string Expression );
  void AddNewNuisanceParameters( std::string NewNuisanceParameters );
  void AddNewGlobalObservable( std::string NewGlobalObservable );

// ____________________________________________________________________________|__________
private:

  std::list< std::string > fExpressions;
  std::list< std::string > fNewNuisanceParameters;
  std::list< std::string > fNewGlobalObservables;

// ____________________________________________________________________________|__________
protected:

  ClassDefOverride(ParametrisationScheme, 1)

};

#endif
