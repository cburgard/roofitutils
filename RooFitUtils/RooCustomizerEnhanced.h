//this file looks like plain C, but it's actually -*- c++ -*-
#ifndef __ROOCUSTOMIZERENHANCED__
#define __ROOCUSTOMIZERENHANCED__

#include <string>
#include <vector>

#include "RooAbsArg.h"
#include "RooArgSet.h"
#include "RooMsgService.h"
#include "RooWorkspace.h"
#include "RooFactoryWSTool.h"
#include "RooRealConstant.h"

//_____________________________________________________________________________

namespace RooFitUtils {

class RooCustomizerEnhanced : public TNamed, public RooPrintable {

public:

  // Constructors, assignment etc
  
  RooCustomizerEnhanced() { _masterLeafListIter=0; _masterBranchListIter=0; } // dummy null constructor so can load from ROOT prompt
  RooCustomizerEnhanced(const RooAbsArg& pdf, const char* name, const char* wild=0) ;
  virtual ~RooCustomizerEnhanced() ;
  
  void setOwning(Bool_t flag) { 
    // If flag is true, make customizer own all created components
    _owning = flag ; 
  }
  
  void replaceArg(const RooAbsArg& orig, const RooAbsArg& subst) ;
  RooAbsArg* build(Bool_t verbose=kFALSE) ;

  const RooArgSet& cloneBranchList() const { 
    // Return list of cloned branch nodes
    return *_cloneBranchList ; 
  }
  const RooArgSet& cloneLeafList() const { 
    // Return list of cloned leaf nodes
    return *_cloneNodeListOwned ; 
  }

  // Printing interface 
  virtual void printName(std::ostream& os) const override;
  virtual void printTitle(std::ostream& os) const override;
  virtual void printClassName(std::ostream& os) const override;
  virtual void printArgs(std::ostream& os) const override;
  virtual void printMultiline(std::ostream& os, Int_t content, Bool_t verbose=kFALSE, TString indent= "") const override;

  inline virtual void Print(Option_t *options= 0) const override{
    // Printing interface
    printStream(defaultPrintStream(),defaultPrintContents(options),defaultPrintStyle(options));
  }

  // Releases ownership of list of cloned branch nodes
  void setCloneBranchSet(RooArgSet& cloneBranchSet) ;

  // Factory interface
  class CustIFace : public RooFactoryWSTool::IFace {
  public:
    virtual ~CustIFace() {} ;
    std::string create(RooFactoryWSTool& ft, const char* typeName, const char* instanceName, std::vector<std::string> args) ;
  } ;

protected:
  
  RooCustomizerEnhanced(const RooCustomizerEnhanced&) ;
  void initialize() ;
  
  RooAbsArg* doBuild(Bool_t verbose) ;

  Bool_t _owning ;  // If true we own all created components
  TString _name ;   // Name of this object

  TList _replaceArgList ; // List of RooAbsArgs to be replaced
  TList _replaceSubList ; // List of replacement RooAbsArgs

  // Master nodes are not owned
  RooAbsArg* _masterPdf ;             // Pointer to input p.d.f

  TIterator* _masterLeafListIter ;    // Iterator over leaf list
  TIterator* _masterBranchListIter ;  // Iterator over branch list

  RooArgSet  _masterBranchList ;      // List of branch nodes
  RooArgSet  _masterLeafList ;        // List of leaf nodes

  RooArgSet  _internalCloneBranchList ; // List of branches of internal clone
  RooArgSet* _cloneBranchList ;         // Pointer to list of cloned branches used

  // Cloned leafs are owned by the user supplied list in the ctor
  RooArgSet* _cloneNodeListAll ;        // List of all cloned nodes
  RooArgSet* _cloneNodeListOwned ;      // List of owned cloned nodes

  std::string _replaceWild ; // Wildcard for replacement

  ClassDefOverride(RooCustomizerEnhanced,0) // Editing tool for RooAbsArg composite object expressions
} ;

}

#endif
