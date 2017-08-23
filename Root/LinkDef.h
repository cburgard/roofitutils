//this file is -*- c++ -*- 
#ifndef __RooFitUtilsDICT__
#define __RooFitUtilsDICT__

#include "RooFitUtils/AbsMeasurement.h"  
#include "RooFitUtils/CombinedMeasurement.h"  
#include "RooFitUtils/EditWorkspaces.h"     
#include "RooFitUtils/Measurement.h"            
#include "RooFitUtils/ParametrisationSequence.h"  
#include "RooFitUtils/WildcardList.h"
#include "RooFitUtils/Channel.h"         
#include "RooFitUtils/Log.h"        
#include "RooFitUtils/Utils.h"        
#include "RooFitUtils/ExtendedModel.h"         
#include "RooFitUtils/CorrelationScheme.h"    
#include "RooFitUtils/ExtendedMinimizer.h"
#include "RooFitUtils/RooCustomizerEnhanced.h"  
#include "RooFitUtils/ParametrisationScheme.h"  
#include "RooFitUtils/RenamingMap.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ class TOwnedList+ ;
#pragma link C++ class ExtendedMinimizer+ ;
#pragma link C++ class ExtendedModel+ ;
#pragma link C++ class Log+ ;
#pragma link C++ class ParametrisationScheme+ ;
#pragma link C++ class ParametrisationSequence+ ;
#pragma link C++ class RenamingMap+ ;
#pragma link C++ class CorrelationScheme+ ;
#pragma link C++ class Channel+ ;
#pragma link C++ class AbsMeasurement+ ;
#pragma link C++ class Measurement+ ;
#pragma link C++ class CombinedMeasurement+ ;
#pragma link C++ class RooCustomizerEnhanced+ ;

#pragma link C++ function editws ;
#pragma link C++ function fixRooStarMomentMorph ;

#endif
#endif
