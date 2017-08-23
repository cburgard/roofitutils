#ifndef PARAMETRISATIONSEQUENCE
#define PARAMETRISATIONSEQUENCE

#include <string>
#include <list>
#include <iostream>
#include <sstream>

#include <TNamed.h>

#include "RooMsgService.h"

#include "ParametrisationScheme.h"

using namespace std;
using namespace RooFit;

class ParametrisationSequence : public TNamed {

// ____________________________________________________________________________|__________
public:

  // Constructor and destructor
  ParametrisationSequence( std::string SequenceName = "", std::string ParametersOfInterest = "" );
  ~ParametrisationSequence();

  // Accessors
  void SetSequence( const std::list< ParametrisationScheme >& Sequence ) { fSequence = Sequence; }
  std::list< ParametrisationScheme > GetSequence() { return fSequence; }

  void SetParametersOfInterest( std::string ParametersOfInterest ) { fParametersOfInterest = ParametersOfInterest; }
  void AppendParametersOfInterest( std::string ParametersOfInterest ) { fParametersOfInterest += "," + ParametersOfInterest; }
  std::string GetParametersOfInterest() { return fParametersOfInterest; }

  // Steering
  void AddScheme( ParametrisationScheme& Scheme );

// ____________________________________________________________________________|__________
private:

  std::list< ParametrisationScheme > fSequence;
  std::string fParametersOfInterest;

// ____________________________________________________________________________|__________
protected:

  ClassDefOverride(ParametrisationSequence, 1)

};

#endif
