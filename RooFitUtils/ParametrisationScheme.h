//this file looks like plain C, but it's actually -*- c++ -*-
#ifndef PARAMETRISATIONSCHEME
#define PARAMETRISATIONSCHEME

#include <string>
#include <list>
#include <iostream>
#include <sstream>

#include "TNamed.h"

#include "RooMsgService.h"

namespace RooFitUtils { 
class ParametrisationScheme : public TNamed {

// ____________________________________________________________________________|__________
public:

  // Constructor and destructor
  ParametrisationScheme( const std::string& SchemeName = "" );
  ~ParametrisationScheme();

  // Accessors
  void SetExpressions( const std::list< std::string >& Expressions ) { fExpressions = Expressions; }
  std::list< std::string > GetExpressions() { return fExpressions; }

  void SetNewNuisanceParameters( const std::list< std::string >& NewNuisanceParameters ) { fNewNuisanceParameters = NewNuisanceParameters; }
  std::list< std::string > GetNewNuisanceParameters() { return fNewNuisanceParameters; }

  void SetNewGlobalObservables( const std::list< std::string >& NewGlobalObservables ) { fNewGlobalObservables = NewGlobalObservables; }
  std::list< std::string > GetNewGlobalObservables() { return fNewGlobalObservables; }

  // Steering
  void AddExpression( const std::string& Expression );
  void AddNewNuisanceParameters( const std::string& NewNuisanceParameters );
  void AddNewGlobalObservable( const std::string& NewGlobalObservable );

// ____________________________________________________________________________|__________
private:

  std::list< std::string > fExpressions;
  std::list< std::string > fNewNuisanceParameters;
  std::list< std::string > fNewGlobalObservables;

// ____________________________________________________________________________|__________
protected:

  ClassDefOverride(ParametrisationScheme, 1)

};

}

#endif
