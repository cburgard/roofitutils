/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id: RooCustomizerEnhanced.cxx 239507 2016-06-21
 *10:50:49Z adye $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// RooCustomizerEnhanced is a factory class to produce clones
// of a prototype composite PDF object with the same structure but
// different leaf servers (parameters or dependents)
//
// RooCustomizerEnhanced supports only replacement modifications:
//
// -> replace(leaf_arg,repl_arg)
// replaces each occurence of leaf_arg with repl_arg in the composite pdf.
//
// In addition to replacing leaf nodes, RooCustomizerEnhanced clones all branch
// nodes that depend directly or indirectly on modified leaf nodes, so
// that the input pdf is untouched by each build operation.
//
// The customizer owns all the branch nodes including the returned top
// level node, so the customizer should live as longs as the cloned
// composites are needed.
//
// Can also run from the workspace factory using:
//   SEDIT::name(orig,{wild},origNode=substNode,...)
//   -- Like EDIT::name(orig,origNode=substNode,...), but limits substitutions
//   to components that depend on origNode and match the wildcard
//   (comma-separated names)

//_____________________________________________________________________________

#include "RooFitUtils/RooCustomizerEnhanced.h"
#include "RooFitUtils/Utils.h"
#include <iostream>
#include <stdexcept>

ClassImp(RooFitUtils::RooCustomizerEnhanced)

    namespace {
  static Int_t init() {
    RooFactoryWSTool::IFace *iface =
        new RooFitUtils::RooCustomizerEnhanced::CustIFace;
    RooFactoryWSTool::registerSpecial("SEDIT", iface);
    return 0;
  }
  static Int_t dummy = init();
}

//_____________________________________________________________________________

RooFitUtils::RooCustomizerEnhanced::RooCustomizerEnhanced(const RooAbsArg &pdf,
                                                          const char *name,
                                                          const char *wild)
    : TNamed(pdf.GetName(), pdf.GetTitle()), _owning(kFALSE), _name(name),
      _masterPdf((RooAbsArg *)&pdf), _masterBranchList("masterBranchList"),
      _masterLeafList("masterLeafList"),
      _internalCloneBranchList("cloneBranchList"), _cloneNodeListAll(0),
      _cloneNodeListOwned(0), _replaceWild(wild) {
  // Customizers created by this constructor
  // offer only the replace() method. The supplied 'name' is used as
  // suffix for any cloned branch nodes

  _masterBranchList.setHashTableSize(1000);
  _masterLeafList.setHashTableSize(1000);

  _cloneBranchList = &_internalCloneBranchList;
  _cloneBranchList->setHashTableSize(1000);

  initialize();
}

//_____________________________________________________________________________

void RooFitUtils::RooCustomizerEnhanced::initialize() {
  // Initialize the customizer

  _masterPdf->leafNodeServerList(&_masterLeafList);
  _masterPdf->branchNodeServerList(&_masterBranchList);
}

//_____________________________________________________________________________

RooFitUtils::RooCustomizerEnhanced::~RooCustomizerEnhanced() {
  // Destructor
}

//_____________________________________________________________________________

void RooFitUtils::RooCustomizerEnhanced::replaceArg(const RooAbsArg &orig,
                                                    const RooAbsArg &subst) {
  // Replace any occurence of arg 'orig' with arg 'subst'

  if (_replaceArgList.FindObject(orig.GetName())) {
    coutE(InputArguments) << "RooCustomizerEnhanced(" << GetName()
                          << ") ERROR: multiple replacement rules defined for "
                          << orig.GetName() << " only using first rule"
                          << std::endl;
    return;
  }

  _replaceArgList.Add((RooAbsArg *)&orig);
  _replaceSubList.Add((RooAbsArg *)&subst);
}

//_____________________________________________________________________________

RooAbsArg *RooFitUtils::RooCustomizerEnhanced::build(Bool_t verbose) {
  // Build a clone of the prototype executing all registered 'replace' rules
  // If verbose is set a message is printed for each leaf or branch node
  // modification. The returned head node owns all cloned branch nodes
  // that were created in the cloning proces

  // Execute build
  RooAbsArg *ret = doBuild(verbose);

  // Make root object own all cloned nodes

  // First make list of all objects that were created
  RooArgSet allOwned;
  if (_cloneNodeListOwned) {
    allOwned.add(*_cloneNodeListOwned);
  }
  allOwned.add(*_cloneBranchList);

  // Remove head node from list
  allOwned.remove(*ret);

  // If list with owned objects is not empty, assign
  // head node as owner
  if (allOwned.getSize() > 0) {
    ret->addOwnedComponents(allOwned);
  }

  return ret;
}

//_____________________________________________________________________________

RooAbsArg *RooFitUtils::RooCustomizerEnhanced::doBuild(Bool_t verbose) {
  // Back-end implementation of the p.d.f building functionality

  // Find nodes that must be replaced according to provided description, Clone
  // nodes, change their names
  RooArgSet masterNodesToBeReplaced("masterNodesToBeReplaced");
  RooArgSet masterReplacementNodes("masterReplacementNodes");
  RooArgSet clonedMasterNodes("clonedMasterNodes");

  masterNodesToBeReplaced.setHashTableSize(1000);
  masterReplacementNodes.setHashTableSize(1000);
  clonedMasterNodes.setHashTableSize(1000);

  _masterLeafListIter->Reset();
  RooAbsArg *node;

  RooArgSet nodeList(_masterLeafList);
  nodeList.setHashTableSize(1000);

  nodeList.add(_masterBranchList);

  RooArgSet *wilds = 0;
  if (_replaceWild.size()) {
    wilds = (RooArgSet *)_masterBranchList.selectByName(_replaceWild.c_str());
    if (verbose)
      coutI(ObjectHandling)
          << "RooCustomizerEnhanced::build(" << _masterPdf->GetName()
          << ") selected branches: " << wilds->contentsString() << std::endl;
  }

  //   cout << "loop over " << nodeList.getSize() << " nodes" << std::endl ;
  for(auto* node : nodeList){
    RooAbsArg *ReplaceArg =
        (RooAbsArg *)_replaceArgList.FindObject(node->GetName());
    if (ReplaceArg) {
      RooAbsArg *substArg =
          (RooAbsArg *)_replaceSubList.At(_replaceArgList.IndexOf(ReplaceArg));
      if (verbose) {
        coutI(ObjectHandling)
            << "RooCustomizerEnhanced::build(" << _masterPdf->GetName()
            << "): tree node " << node->GetName() << " will be replaced by "
            << substArg->GetName() << std::endl;
      }

      // Affix attribute with old name to support name changing server redirect
      TString nameAttrib("ORIGNAME:");
      nameAttrib.Append(node->GetName());
      substArg->setAttribute(nameAttrib);

      // Add to list
      masterNodesToBeReplaced.add(*node);
      masterReplacementNodes.add(*substArg);
    }
  }

  // Find branches that are affected and must be cloned
  RooArgSet masterBranchesToBeCloned("masterBranchesToBeCloned");
  masterBranchesToBeCloned.setHashTableSize(1000);
  _masterBranchListIter->Reset();
  RooAbsArg *branch;
  while ((branch = (RooAbsArg *)_masterBranchListIter->Next())) {

    if (masterNodesToBeReplaced.find(branch->GetName())) {
      if (verbose) {
        coutI(ObjectHandling)
            << "RooCustomizerEnhanced::build(" << _masterPdf->GetName()
            << ") Branch node " << branch->GetName() << " is already replaced"
            << std::endl;
      }
      continue;
    }

    if (branch->dependsOn(masterNodesToBeReplaced)) {
      if (wilds && !branch->dependsOn(*wilds))
        continue;
      if (verbose) {
        coutI(ObjectHandling)
            << "RooCustomizerEnhanced::build(" << _masterPdf->GetName()
            << ") Branch node " << branch->IsA()->GetName()
            << "::" << branch->GetName()
            << " cloned: depends on a replaced parameter" << std::endl;
      }
      masterBranchesToBeCloned.add(*branch);
    }
  }

  delete wilds;
  wilds = 0;

  // Clone branches, changes their names
  RooAbsArg *cloneTopPdf = 0;
  RooArgSet clonedMasterBranches("clonedMasterBranches");
  clonedMasterBranches.setHashTableSize(1000);
  for(auto branch : masterBranchesToBeCloned){
    // Affix attribute with old name to clone to support name changing server
    // redirect
    RooAbsArg *clone = (RooAbsArg *)branch->Clone(branch->GetName());
    clone->setStringAttribute("factory_tag", 0);
    TString nameAttrib("ORIGNAME:");
    nameAttrib.Append(branch->GetName());
    clone->setAttribute(nameAttrib);

    if (!clone->getStringAttribute("origName")) {
      clone->setStringAttribute("origName", branch->GetName());
    }

    clonedMasterBranches.add(*clone);

    // Save pointer to clone of top-level pdf
    if (branch == _masterPdf)
      cloneTopPdf = (RooAbsArg *)clone;
  }

  if (_owning) {
    _cloneBranchList->addOwned(clonedMasterBranches);
  } else {
    _cloneBranchList->add(clonedMasterBranches);
  }

  // Reconnect cloned branches to each other and to cloned nodess
  for(auto* branch : clonedMasterBranches){
    branch->redirectServers(clonedMasterBranches, kFALSE, kTRUE);
    branch->redirectServers(clonedMasterNodes, kFALSE, kTRUE);
    branch->redirectServers(masterReplacementNodes, kFALSE, kTRUE);
  }

  return cloneTopPdf ? cloneTopPdf : _masterPdf;
}

//_____________________________________________________________________________

void RooFitUtils::RooCustomizerEnhanced::printName(std::ostream &os) const {
  // Print name of customizer
  os << GetName();
}

//_____________________________________________________________________________

void RooFitUtils::RooCustomizerEnhanced::printTitle(std::ostream &os) const {
  // Print title of customizer
  os << GetTitle();
}

//_____________________________________________________________________________

void RooFitUtils::RooCustomizerEnhanced::printClassName(
    std::ostream &os) const {
  // Print class name of customizer
  os << IsA()->GetName();
}

//_____________________________________________________________________________

void RooFitUtils::RooCustomizerEnhanced::printArgs(std::ostream &os) const {
  // Print arguments of customizer, i.e. input p.d.f
  os << "[ masterPdf=" << _masterPdf->GetName();
  os << " ]";
}

//_____________________________________________________________________________

void RooFitUtils::RooCustomizerEnhanced::printMultiline(std::ostream &os,
                                                        Int_t /*content*/,
                                                        Bool_t /*verbose*/,
                                                        TString indent) const {
  // Print customizer configuration details

  os << indent << "RooCustomizerEnhanced for " << _masterPdf->GetName()
     << std::endl;

  Int_t i, nrepl = _replaceArgList.GetSize();
  if (nrepl > 0) {
    os << indent << "  Replacement rules:" << std::endl;
    for (i = 0; i < nrepl; i++) {
      os << indent << "   " << _replaceSubList.At(i)->GetName() << " replaces "
         << _replaceArgList.At(i)->GetName() << std::endl;
    }
  }

  return;
}

//_____________________________________________________________________________

void RooFitUtils::RooCustomizerEnhanced::setCloneBranchSet(
    RooArgSet &cloneBranchSet) {
  // Install the input RooArgSet as container in which all cloned branches
  // will be stored

  _cloneBranchList = &cloneBranchSet;
  _cloneBranchList->setHashTableSize(1000);
}

namespace {
std::string Join(const std::vector<std::string> &v,
                 const std::string &sep = ",", bool trim = true) {
  std::string s;
  size_t n = v.size();
  if (trim)
    for (; n > 0; n--)
      if (v[n - 1].size() > 0)
        break;
  for (size_t i = 0; i < n; i++) {
    if (i > 0)
      s += sep;
    s += v[i].c_str();
  }
  return s;
}
}

//_____________________________________________________________________________
std::string RooFitUtils::RooCustomizerEnhanced::CustIFace::create(
    RooFactoryWSTool &ft, const char *typeName, const char *instanceName,
    std::vector<std::string> args) {

  // Check number of arguments
  if (args.size() < 2) {
    throw std::runtime_error(
        Form("RooCustomizerEnhanced::CustIFace::create() ERROR: expect at "
             "least 2 arguments for SEDIT: the input object and at least one "
             "$Replace() rule"));
  }

  if (std::string(typeName) != "SEDIT") {
    throw std::runtime_error(Form("RooCustomizerEnhanced::CustIFace::create() "
                                  "ERROR: unknown type requested: %s",
                                  typeName));
  }

  // Check that first arg exists as RooAbsArg
  RooAbsArg *arg = ft.ws().arg(args[0].c_str());
  if (!arg) {
    throw std::runtime_error(Form("RooCustomizerEnhanced::CustIFace::create() "
                                  "ERROR: input RooAbsArg %s does not exist",
                                  args[0].c_str()));
  }

  // If name of new object is same as original, execute in sterile mode (i.e no
  // suffixes attached), and rename original nodes in workspace upon import
  if (args[0] == instanceName) {
    instanceName = 0;
  }

  size_t i = 1;
  Bool_t verbose = kFALSE;
  if (args[i] == "VERBOSE") {
    verbose = kTRUE;
    i++;
  }

  if (verbose)
    std::cout << typeName << "::" << instanceName << " (" << Join(args, ", ")
              << ")" << std::endl;

  size_t l = args[i].size();
  std::string wild;
  if (l >= 2 && ((args[i][0] == '{' && args[i][l - 1] == '}') ||
                 (args[i][0] == '\'' && args[i][l - 1] == '\''))) {
    wild = args[i].substr(1, l - 2);
    i++;
  }

  if (args.size() - i < 1) {
    throw std::runtime_error(
        Form("RooCustomizerEnhanced::CustIFace::create() ERROR: expect at "
             "least one $Replace() rule"));
  }

  // Create a customizer
  RooCustomizerEnhanced cust(*arg, instanceName, wild.c_str());

  for (; i < args.size(); i++) {
    char buf[1024];
    strlcpy(buf, args[i].c_str(), 1024);
    char *sep = strchr(buf, '=');
    if (!sep) {
      throw std::runtime_error(
          Form("RooCustomizerEnhanced::CustIFace::create() ERROR: unknown "
               "argument: %s, expect form orig=subst",
               args[i].c_str()));
    }
    *sep = 0;
    RooAbsArg *orig = ft.ws().arg(buf);
    RooAbsArg *subst(0);
    if (std::string(sep + 1).find("$REMOVE") == 0) {

      // Create a removal dummy ;
      subst = &RooRealConstant::removalDummy();

      // If removal instructed was annotated with target node, encode these in
      // removal dummy
      char *sep2 = strchr(sep + 1, '(');
      if (sep2) {
        char buf2[1024];
        strlcpy(buf2, sep2 + 1, 1024);
        char *saveptr;
        char *tok = strtok_r(buf2, ",)", &saveptr);
        while (tok) {
          // cout << "$REMOVE is restricted to " << tok << std::endl ;
          subst->setAttribute(Form("REMOVE_FROM_%s", tok));
          tok = strtok_r(0, ",)", &saveptr);
        }
      } else {
        // Otherwise mark as universal removal node
        subst->setAttribute("REMOVE_ALL");
      }

    } else {
      subst = ft.ws().arg(sep + 1);
    }
    if (orig && subst) {
      cust.replaceArg(*orig, *subst);
    } else {
      oocoutW((TObject *)0, ObjectHandling)
          << "RooCustomizerEnhanced::CustIFace::create() WARNING: input or "
             "replacement of a replacement operation not found, operation "
             "ignored"
          << std::endl;
    }
  }

  // Build the desired edited object
  RooAbsArg *targ = cust.build(verbose);
  if (!targ) {
    throw std::runtime_error(
        Form("RooCustomizerEnhanced::CustIFace::create() ERROR in customizer "
             "build, object %snot created",
             instanceName));
  }

  // Import the object into the workspace
  if (instanceName) {
    // Set the desired name of the top level node
    targ->SetName(instanceName);
    ft.ws().import(cust.cloneBranchList(), RooFit::Silence(),
                   RooFit::NoRecursion(kTRUE));
  } else {
    ft.ws().import(cust.cloneBranchList(), RooFit::Silence(),
                   RooFit::RenameConflictNodes("orig", 1),
                   RooFit::NoRecursion(kTRUE));
  }

  return instanceName ? instanceName : targ->GetName();
}
