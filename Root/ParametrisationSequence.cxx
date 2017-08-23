#include "RooFitUtils/ParametrisationSequence.h"

// ____________________________________________________________________________|__________
// Constructor
ParametrisationSequence::ParametrisationSequence( std::string SchemeName, std::string ParametersOfInterest )
  :
  TNamed( SchemeName.c_str(), SchemeName.c_str() ),
  fParametersOfInterest( ParametersOfInterest )
{
  coutP(InputArguments) << "ParametrisationSequence::ParametrisationSequence(" << fName <<") new ParametrisationSequence" << endl;
}

// ____________________________________________________________________________|__________
// Destructor
ParametrisationSequence::~ParametrisationSequence()
{

}

// ____________________________________________________________________________|__________
// Interface to add expressions for the re-parametrisation
void ParametrisationSequence::AddScheme( ParametrisationScheme& Scheme )
{
  coutP(InputArguments) << "ParametrisationSequence::AddScheme(" << fName <<") adding expression " << Scheme.GetName() << endl;
  fSequence.push_back(Scheme);
}
