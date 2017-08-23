#include "RooFitUtils/ParametrisationScheme.h"

// ____________________________________________________________________________|__________
// Constructor
ParametrisationScheme::ParametrisationScheme( std::string SchemeName )
  :
  TNamed( SchemeName.c_str(), SchemeName.c_str() )
{
  coutP(InputArguments) << "ParametrisationScheme::ParametrisationScheme(" << fName <<") new ParametrisationScheme" << endl;
}

// ____________________________________________________________________________|__________
// Destructor
ParametrisationScheme::~ParametrisationScheme()
{

}

// ____________________________________________________________________________|__________
// Interface to add expressions for the re-parametrisation
void ParametrisationScheme::AddExpression( std::string Expression )
{
  coutP(InputArguments) << "ParametrisationScheme::AddExpression(" << fName <<") adding expression " << Expression << endl;
  fExpressions.push_back(Expression);
}

// ____________________________________________________________________________|__________
// Interface to add nuisance parameters for the re-parametrisation
void ParametrisationScheme::AddNewNuisanceParameters( std::string NewNuisanceParameter )
{
  coutP(InputArguments) << "ParametrisationScheme::AddNewNuisanceParameters(" << fName <<") adding nuisance parameter " << NewNuisanceParameter << endl;
  fNewNuisanceParameters.push_back(NewNuisanceParameter);
}

// ____________________________________________________________________________|__________
// Interface to add global observables for the re-parametrisation
void ParametrisationScheme::AddNewGlobalObservable( std::string NewGlobalObservable )
{
  coutP(InputArguments) << "ParametrisationScheme::AddNewGlobalObservable(" << fName <<") adding global observable " << NewGlobalObservable << endl;
  fNewGlobalObservables.push_back(NewGlobalObservable);
}
