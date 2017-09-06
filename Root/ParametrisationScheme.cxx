#include "RooFitUtils/ParametrisationScheme.h"

// ____________________________________________________________________________|__________

RooFitUtils::ParametrisationScheme::ParametrisationScheme( const std::string& SchemeName )
  :
  TNamed( SchemeName.c_str(), SchemeName.c_str() )
{
  // Constructor
  coutP(InputArguments) << "ParametrisationScheme::ParametrisationScheme(" << fName <<") new ParametrisationScheme" << std::endl;
}

// ____________________________________________________________________________|__________

RooFitUtils::ParametrisationScheme::~ParametrisationScheme()
{
  // Destructor
}

// ____________________________________________________________________________|__________

void RooFitUtils::ParametrisationScheme::AddExpression( const std::string& Expression )
{
  // Interface to add expressions for the re-parametrisation
  coutP(InputArguments) << "ParametrisationScheme::AddExpression(" << fName <<") adding expression " << Expression << std::endl;
  fExpressions.push_back(Expression);
}

// ____________________________________________________________________________|__________

void RooFitUtils::ParametrisationScheme::AddNewNuisanceParameters( const std::string& NewNuisanceParameter )
{
  // Interface to add nuisance parameters for the re-parametrisation
  coutP(InputArguments) << "ParametrisationScheme::AddNewNuisanceParameters(" << fName <<") adding nuisance parameter " << NewNuisanceParameter << std::endl;
  fNewNuisanceParameters.push_back(NewNuisanceParameter);
}

// ____________________________________________________________________________|__________

void RooFitUtils::ParametrisationScheme::AddNewGlobalObservable( const std::string& NewGlobalObservable )
{
  // Interface to add global observables for the re-parametrisation
  coutP(InputArguments) << "ParametrisationScheme::AddNewGlobalObservable(" << fName <<") adding global observable " << NewGlobalObservable << std::endl;
  fNewGlobalObservables.push_back(NewGlobalObservable);
}
