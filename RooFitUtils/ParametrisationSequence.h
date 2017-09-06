//this file looks like plain C, but it's actually -*- c++ -*-
#ifndef PARAMETRISATIONSEQUENCE
#define PARAMETRISATIONSEQUENCE

#include <string>
#include <list>
#include <iostream>
#include <sstream>

#include <TNamed.h>

#include "RooMsgService.h"

#include "ParametrisationScheme.h"

namespace RooFitUtils {

class ParametrisationSequence : public TNamed {

public:

  // Constructor and destructor
  ParametrisationSequence( const std::string& SequenceName = "", const std::string& ParametersOfInterest = "" );
  ~ParametrisationSequence();

  // Accessors
  void SetSequence( const std::list< ParametrisationScheme >& Sequence ) { fSequence = Sequence; }
  std::list< ParametrisationScheme > GetSequence() { return fSequence; }

  void SetParametersOfInterest( std::string ParametersOfInterest ) { fParametersOfInterest = ParametersOfInterest; }
  void AppendParametersOfInterest( std::string ParametersOfInterest ) { fParametersOfInterest += "," + ParametersOfInterest; }
  std::string GetParametersOfInterest() { return fParametersOfInterest; }

  // Steering
  void AddScheme( ParametrisationScheme& Scheme );

private:

  std::list< ParametrisationScheme > fSequence;
  std::string fParametersOfInterest;

protected:

  ClassDefOverride(ParametrisationSequence, 1)

};

}

#endif
