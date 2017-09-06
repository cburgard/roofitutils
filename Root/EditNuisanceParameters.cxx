/*==============================================================================
$Id: EditNuisanceParameters.cxx 239512 2016-06-21 15:41:35Z adye $

Define additional workspace factory editing commands to manipulate nuisance
parameter terms.
Can be used with (or without) editws.cc to extend the editing commands provided
by standard
RooFit's RooFactoryWSTool.

Usage (see editws.cc):
  root -b -q INPUT_FILE npsplit.cc+ editws.cc+'(OUTPUT_FILE, EDIT_SPEC)'

Eg.
  root -b -q wsin.root npsplit.cc+ editws.cc+'("wsout.root",
"RemoveDuplicates()")'

The default script (run without parameters), registers custom RooFactoryWSTool
commands to
extend those documented in
http://root.cern.ch/root/html/RooFactoryWSTool.html#RooFactoryWSTool:process
with the following:

NP::name()
   -- Create nuisance parameter variable name, Gaussian constraint term
name_Pdf, and
      global observable name_In. Add to ModelConfig.
NP::name(constraint,global_obs,sigma)
   -- Specify constraint term, and optionally global observable and Gaussian
sigma

FlexibleInterpVar::name(np1,np2,np3,nominal,code,err1,err2,err3)
   -- Create FlexibleInterpVar with NPs np1,np2,np3, nominal value,
interpolation code, and errors err1,err2,err3
FlexibleInterpVar::name(np1,np2,np3,nominal,code,lo1,hi1,lo2,hi2,lo3,hi3)
   -- Create FlexibleInterpVar with NPs np1,np2,np3, nominal value,
interpolation code, and errors +hi1-lo1, +hi2-lo2, +hi3-lo3.

edit(origNode=substNode,...)
   -- Like EDIT::OLDPDF(orig,origNode=substNode,...), but always edits PDF
in-place (orig=OLDPDF) and works with QueueEdits

ReplaceConstraint(oldNP,newNP1,newNP2,...,{wildcards})
   -- edits PDF to replace oldNP's constraint term with new NPs. The new NPs can
be specified as the NP variable,
      or (less ambiguously) the constraint PDF. The optional wildcards restrict
which product terms to update.
      If wildcards include [] specs, specify them as 'wildcards' instead.

ReplaceResponse(oldNP,newFlexibleInterpVar,{wildcards})
   -- edits PDF to replace oldNP's FlexibleInterpVar term with new values (given
as a template FlexibleInterpVar).
      The optional wildcards restrict which terms to replace.
ReplaceResponse(oldNP,newTerm1,newTerm2,...,{wildcards})
   -- edits PDF to replace oldNP's response terms with new ones of the same
type.
      The optional wildcards restrict which terms to replace.
      If wildcards include [] specs, specify them as 'wildcards' instead.

RemoveNP(NP)
   -- Equivalent to: ReplaceConstraint(NP); ReplaceResponse(NP)

RemoveDuplicates()
   -- edits PDF to remove duplicate FlexibleInterpVar terms

QueueEdits()
   -- queue following edits (edit, ReplaceConstraint, ReplaceResponse, and
RemoveDuplicates) to perform all in one go
CommitEdits()
   -- execute queued edits

The ReplaceConstraint, ReplaceResponse, RemoveNP, and (optionally) NP use the
ModelConfig.
This is set automatically by editws.cc, or can be specified when registering
npsplit(mc).
All these editing commands (edit, ReplaceConstraint, ReplaceResponse, RemoveNP,
RemoveDuplicates, and CommitEdits)
operate by default on the ModelConfig's pdf. If ModelConfig is not defined, or a
different pdf should be edited,
specify as edit::MyPdf().

ReplaceConstraint and ReplaceResponse create new PDF terms (with names like
term_NEW or term_NEW2)
with updated values. These are replaced in the PDF with their new names. For
neatness' sake,
editws.cc will change names back to the originals once they are imported into
the new workspace.

Author: Tim Adye <T.J.Adye@rl.ac.uk>, 23 Oct 2014.

==============================================================================*/

#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooCustomizer.h"
#include "RooFactoryWSTool.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooRealVar.h"
#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "RooStats/ModelConfig.h"
#include "RooWorkspace.h"
#include "TError.h"
#include "TPRegexp.h"
#include "TRegexp.h"
#include "TString.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#if defined(__ACLIC__) || defined(__CINT__) || defined(__CLING__) ||           \
    defined(__MAKECINT__)
#include "RooStats/ModelConfig.h"
#include <string>
namespace RooFitUtils {
// Place for editws.cc to fill in current ModelConfig.
// This is used if CustIFace_npsplit is not initialised with a different
// ModelConfig.
#ifndef ROOFITUTILS_DEFINED
RooStats::ModelConfig *currentMC = 0;
std::string *renameSpec = 0;
#else
extern RooStats::ModelConfig *currentMC;
extern std::string *renameSpec;
#endif
}
#else
#include "RooFitUtils.h"
#endif

#ifdef ROOFITUTILS_DEFINED
#undef INFO
#undef ERROR
#undef WARN
#else
#define ROOFITUTILS_DEFINED
#endif

#define INFO(...) Info("npsplit", __VA_ARGS__)
#define ERROR(...) Error("npsplit", __VA_ARGS__)
#define WARN(...) Warning("npsplit", __VA_ARGS__)

namespace {

std::string Print(const RooAbsArg &a, Int_t level = 2) {
  std::ostringstream os;
  if (level >= 1) {
    a.printClassName(os);
    os << "::";
  }
  if (level >= 0)
    a.printName(os);
  if (level >= 2)
    a.printArgs(os);
  if (level >= 3) {
    os << " = ";
    a.printValue(os);
  }
  if (level >= 4) {
    os << " ";
    a.printExtras(os);
  }
  return os.str();
}

std::string Print(const RooAbsCollection &c, Int_t level = 0) {
  std::ostringstream os;
  int index = 0;
  for (RooFIter it = c.fwdIterator(); RooAbsArg *a = it.next();) {
    if (level <= 1) {
      if (index > 0)
        os << ", ";
      os << Print(*a, level);
    } else {
      os << Form("%3d) ", index + 1);
      os << Print(*a, level) << std::endl;
    }
    index++;
  }
  return os.str();
}

Bool_t isSubset(const RooAbsCollection &set, const RooAbsCollection &sub,
                Bool_t matchByNameOnly = kTRUE) {
  if (sub.getSize() <= 0)
    return kFALSE; // empty subset -> false (more useful, if not mathematically
                   // true)
  RooFIter it = sub.fwdIterator();
  const RooAbsArg *a;
  while ((a = it.next())) {
    if (!(matchByNameOnly ? set.contains(*a) : set.containsInstance(*a)))
      return kFALSE;
  }
  return kTRUE;
}

TString newName(const RooWorkspace *ws, const char *old,
                const char *suffix = "_NEW") {
  TString pre = old;
  Ssiz_t is = pre.Index(TRegexp(Form("%s[0-9]*$", suffix)));
  int num = 1;
  if (is != kNPOS) {
    is += strlen(suffix);
    if (is < pre.Length())
      num = atoi(TString(pre(is, pre.Length()))) + 1;
    else
      num = 2;
    pre.Remove(is);
  } else
    pre += suffix;
  for (;; num++) {
    TString name = pre;
    if (num >= 2)
      name += Form("%d", num);
    if (!ws->arg(name))
      return name;
  }
}

RooAbsArg *import(RooWorkspace *ws, const RooAbsArg &a,
                  Bool_t ifNotPresent = kFALSE, int verbose = 0) {
  if (ifNotPresent) {
    RooAbsArg *b = ws->arg(a.GetName());
    if (b)
      return b;
  }
  if (verbose >= 1)
    std::cout << "import(" << a.GetName() << ")"
              << (ws->arg(a.GetName()) ? " already in workspace" : "") << ""
              << std::endl;
  if (ws->import(a, RooFit::RecycleConflictNodes(), RooFit::Silence())) {
    if (verbose >= 1)
      std::cout << "-> import error" << std::endl;
    return 0;
  }
  RooAbsArg *b = ws->arg(a.GetName());
  if (!b) {
    std::cout << "Could not import into workspace '" << ws << "': " << Print(a)
              << std::endl;
  } else if (verbose >= 2)
    std::cout << "-> import OK" << std::endl;
  return b;
}

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

inline bool Strtod(const std::string &s, double &d) {
  char *p = 0;
  d = strtod(s.c_str(), &p);
  return p && !*p;
}
}

RooProdPdf *editTerm(const RooProdPdf &pdf, const RooArgList &oldterms,
                     const RooArgList &newterms, const char *name = 0,
                     const char *title = 0) {
  if (!name)
    name = Form("%s_NEW", pdf.GetName());
  if (!title)
    title = pdf.GetTitle();
  RooArgList pdfList = pdf.pdfList();
  if (!isSubset(pdfList, oldterms)) {
    WARN("%s does not contain %s", Print(pdf).c_str(), Print(oldterms).c_str());
  }
  pdfList.remove(oldterms, kTRUE);
  pdfList.add(newterms);
  return new RooProdPdf(name, title, pdfList);
}

RooProduct *editTerm(const RooProduct &pdf, const RooArgList &oldterms,
                     const RooArgList &newterms, const char *name = 0,
                     const char *title = 0) {
  if (!name)
    name = Form("%s_NEW", pdf.GetName());
  if (!title)
    title = pdf.GetTitle();
  RooArgList pdfList = const_cast<RooProduct &>(pdf).components();
  if (!isSubset(pdfList, oldterms)) {
    WARN("%s does not contain %s", Print(pdf).c_str(), Print(oldterms).c_str());
  }
  pdfList.remove(oldterms, kTRUE);
  pdfList.add(newterms);
  return new RooProduct(name, title, pdfList);
}

class FlexibleInterpVarInspector
    : public RooStats::HistFactory::FlexibleInterpVar {
public:
  static const FlexibleInterpVarInspector &
  inspect(const RooStats::HistFactory::FlexibleInterpVar &f) {
    return (const FlexibleInterpVarInspector &)f;
  }
  static const FlexibleInterpVarInspector *inspect(const TObject *o) {
    return (
        const FlexibleInterpVarInspector
            *)dynamic_cast<const RooStats::HistFactory::FlexibleInterpVar *>(o);
  }

  Double_t nominal() const { return _nominal; }
  Int_t index(const RooAbsArg &a) const { return _paramList.index(&a); }
  double low(size_t i) const { return _low[i]; }
  double high(size_t i) const { return _high[i]; }
  int code(size_t i) const { return _interpCode[i]; }
  std::vector<double> &low() { return _low; }
  std::vector<double> &high() { return _high; }
  std::vector<int> &code() { return _interpCode; }
  const RooArgList &paramList() const { return _paramList; }
  void print() const {
    std::cout << ::Print(*this, 4) << std::endl;
    for (Int_t i = 0, n = _paramList.getSize(); i < n; i++)
      std::cout << "  " << _paramList[i].GetName() << " = "
                << printParam(_paramList[i]) << std::endl;
  }

  std::string printParam(const RooAbsArg &a) const {
    Int_t i = index(a);
    if (i < 0) {
      return "no parameter";
    }
    double hi = _high[i] - _nominal, lo = _low[i] - _nominal;
    if (_nominal == 0.0 || _interpCode[i] == 0 ||
        _interpCode[i] == 2) { // linear
      if (fabs(hi + lo) > 1e-7)
        return Form("%g %+g %+g (%d)", _nominal, hi, lo, _interpCode[i]);
      else if (lo > 0.0)
        return Form("%g -/+ %g (%d)", _nominal, lo, _interpCode[i]);
      else
        return Form("%g +/- %g (%d)", _nominal, hi, _interpCode[i]);
    } else { // log
      double fhi = _high[i] / _nominal, flo = _low[i] / _nominal;
      if (fabs(fhi * flo - 1.0) <= 1e-7) {
        if (flo > 1.0)
          return Form("%g -/+ %g%% (%d)", _nominal, 100.0 * log(flo),
                      _interpCode[i]);
        else
          return Form("%g +/- %g%% (%d)", _nominal, 100.0 * log(fhi),
                      _interpCode[i]);
      } else if (fabs(hi + lo) <= 1e-7) {
        if (lo > 0.0)
          return Form("%g -/+ %g (%d)", _nominal, lo, _interpCode[i]);
        else
          return Form("%g +/- %g (%d)", _nominal, hi, _interpCode[i]);
      } else {
        return Form("%g %+g%% %+g%% (%d)", _nominal, 100.0 * log(fhi),
                    100.0 * log(flo), _interpCode[i]);
      }
    }
  }

  RooStats::HistFactory::FlexibleInterpVar *
  edit(const RooArgList &oldnp,
       const RooStats::HistFactory::FlexibleInterpVar &newterm,
       Double_t scale = 1.0, Bool_t copyCode = kFALSE, const char *name = 0,
       const char *title = 0) const {
    const FlexibleInterpVarInspector &fi = inspect(newterm);
    return edit(fi._paramList, fi._low, fi._high,
                (copyCode ? fi._interpCode : std::vector<int>()), oldnp, scale,
                name, title);
  }

  RooStats::HistFactory::FlexibleInterpVar *
  edit(const RooArgList &paramList_, const std::vector<double> &low_,
       const std::vector<double> &high_,
       const std::vector<int> &code_ = std::vector<int>(),
       const RooArgList &remove = RooArgList(), Double_t scale = 1.0,
       const char *name = 0, const char *title = 0) const {
    if (!name)
      name = Form("%s_NEW", GetName());
    if (!title)
      title = GetTitle();
    size_t n = _low.size(), m = low_.size(), nc = code_.size();
    if (m != (size_t)paramList_.getSize() || m != high_.size() || m < nc) {
      ERROR("FlexibleInterpVarInspector::edit(%s) ERROR: input list sizes "
            "don't match",
            name);
      return 0;
    }
    if (n != (size_t)_paramList.getSize() || n != _high.size() ||
        n != _interpCode.size()) {
      ERROR("FlexibleInterpVarInspector::edit(%s) ERROR: existing "
            "FlexibleInterpVar list sizes don't match",
            name);
      return 0;
    }
    RooArgList p = _paramList;
    std::vector<double> l = _low, h = _high;
    std::vector<int> c = _interpCode;

    std::vector<bool> rem(n);
    for (RooFIter it = remove.fwdIterator(); RooAbsArg *a = it.next();) {
      Int_t i = p.index(a);
      if (i < 0)
        ERROR("FlexibleInterpVarInspector::edit(%s) ERROR: parameter %s not "
              "present for removal",
              name, a->GetName());
      else
        rem[i] = true;
    }
    int defcode = n > 0 ? c[0] : 0;
    for (int i = n - 1; i >= 0; i--) {
      if (c[i] != defcode)
        defcode = 0;
      if (rem[i]) {
        p.remove(*p.at(i));
        l.erase(l.begin() + i);
        h.erase(h.begin() + i);
        c.erase(c.begin() + i);
      }
    }
    for (size_t i = 0; i < m; i++) {
      RooAbsArg *a = paramList_.at(i);
      Int_t j = p.index(a);
      Double_t lo = low_[i], hi = high_[i];
      if (scale != 1.0) {
        lo = _nominal - scale * (_nominal - lo);
        hi = _nominal + scale * (hi - _nominal);
      }
      if (j < 0) {
        p.add(*a);
        l.push_back(lo);
        h.push_back(hi);
        c.push_back(i < nc ? code_[i] : defcode);
      } else {
        l[j] = lo;
        h[j] = hi;
        c[j] = i < nc ? code_[i] : defcode;
      }
    }
    RooStats::HistFactory::FlexibleInterpVar *f =
        new RooStats::HistFactory::FlexibleInterpVar(name, title, p, _nominal,
                                                     l, h, c);
    f->setGlobalBoundary(_interpBoundary);
    return f;
  }

  Bool_t equals(const RooStats::HistFactory::FlexibleInterpVar &other,
                const RooAbsCollection *params = 0, int test = 2) const {
    // test=0: just params, 1=+values, 2=+code
    const FlexibleInterpVarInspector &o =
        FlexibleInterpVarInspector::inspect(other);
    if (params) {
      if (!isSubset(_paramList, *params) || !isSubset(o._paramList, *params))
        return kFALSE;
    } else {
      if (!isSubset(_paramList, o._paramList))
        return kFALSE;
      params = &o._paramList;
    }
    if (test <= 0)
      return kTRUE;
    if (_nominal != o._nominal)
      return kFALSE;
    Int_t j = 0;
    for (RooFIter it = params->fwdIterator(); RooAbsArg *a = it.next(); j++) {
      Int_t i = index(*a);
      if (i < 0)
        return kFALSE; // just in case
      if (_low[i] != o._low[j] || _high[i] != o._high[j])
        return kFALSE;
      if (test >= 2 && _interpCode[i] != o._interpCode[j])
        return kFALSE;
    }
    return kTRUE;
  }

  ClassDef(FlexibleInterpVarInspector, 0);
};

int flexmerge(RooCustomizer &cust, const RooAbsPdf &pdf,
              const RooArgSet *params = 0, int verbose = 0) {
  // optimise by combining FlexibleInterpVar terms
  const RooArgSet *comp = pdf.getComponents();
  RooArgList fiv;
  for (RooFIter it = comp->fwdIterator(); RooAbsArg *a = it.next();) {
    if (dynamic_cast<const RooStats::HistFactory::FlexibleInterpVar *>(a))
      fiv.add(*a);
  }
  delete comp;

  Int_t nfiv = fiv.getSize();
  std::vector<bool> isdup(nfiv);
  int nuniq = 0, nrep = 0;
  for (Int_t ia = 0; ia < nfiv; ia++) {
    if (isdup[ia])
      continue;
    nuniq++;
    const FlexibleInterpVarInspector *fa =
        FlexibleInterpVarInspector::inspect(fiv.at(ia));
    RooArgList dup;
    for (Int_t ib = ia + 1; ib < nfiv; ib++) {
      const RooStats::HistFactory::FlexibleInterpVar *fb =
          dynamic_cast<const RooStats::HistFactory::FlexibleInterpVar *>(
              fiv.at(ib));
      if (fa->equals(*fb, params)) {
        dup.add(*fb);
        isdup[ib] = true;
      }
    }
    Int_t ndup = dup.getSize();
    if (!ndup)
      continue;
    if (verbose >= 1) {
      std::cout << "Duplicate terms: " << fa->GetName();
      for (Int_t i = 0; i < ndup; i++)
        std::cout << ' ' << dup[i].GetName();
      std::cout << std::endl;
    }
    for (Int_t i = 0; i < ndup; i++)
      cust.replaceArg(dup[i], *fa);
    nrep += ndup;
  }
  std::cout << "There are " << nuniq << "/" << nfiv
            << " unique FlexibleInterpVar terms." << std::endl;
  return nrep;
}

static bool sortdup(std::pair<std::pair<double, double>, int> a,
                    std::pair<std::pair<double, double>, int> b) {
  return (a.second > b.second); // reverse order in count
}

int npsplit_response(RooWorkspace *ws, RooCustomizer &cust,
                     const RooAbsPdf &pdf, const RooArgList &oldnp,
                     const RooArgList &newterms, const RooArgSet *allnp,
                     const RooArgSet *observables, const char *select = 0,
                     int verbose = 1, const char *suffix = "_NEW") {
  const TClass *newtermClass = 0;
  if (newterms.getSize()) {
    newtermClass = newterms[0].IsA();
    for (Int_t i = 1, n = newterms.getSize(); i < n; i++) {
      if (newterms[i].IsA() != newtermClass) {
        WARN("new terms have different types: %s", Print(newterms, 2).c_str());
        break;
      }
    }
  }

  RooArgSet excl = *allnp;
  excl.remove(oldnp, kTRUE);
  excl.add(*observables);
  Int_t noldnp = oldnp.getSize();

  const RooArgSet *comp = pdf.getComponents();
  if (select && *select) {
    const RooArgSet *sel =
        dynamic_cast<const RooArgSet *>(comp->selectByName(select));
    delete comp;
    comp = sel;
  }
  if (!comp->getSize()) {
    if (select && *select)
      WARN("no PDF terms matching %s", select);
    else
      WARN("no PDF terms");
    delete comp;
    return 0;
  }

  // Find all terms of the right type that depend on oldnp, but not on any of
  // the other NPs,
  // nor any of the global observables (constraint terms).
  RooArgList terms, clients;
  std::vector<RooArgList> clientTerms;
  for (RooFIter it = comp->fwdIterator(); RooAbsArg *a = it.next();) {
    RooArgSet fp;
    if (const RooProdPdf *p1 = dynamic_cast<const RooProdPdf *>(a))
      fp.add(p1->pdfList());
    else if (RooProduct *p2 = dynamic_cast<RooProduct *>(a))
      fp.add(p2->components());
    else
      continue;
    RooArgSet ft, ftnew;
    RooArgList oldnp1 = oldnp;
    for (RooFIter ip = fp.fwdIterator(); RooAbsArg *t = ip.next();) {
      if (newtermClass && t->IsA() != newtermClass)
        continue;
      int iv, ndep = 0;
      for (iv = 0; iv < noldnp; iv++) {
        if (!t->dependsOnValue(oldnp[iv]))
          continue;
        ndep++;
        if (!oldnp1.remove(oldnp[iv], kTRUE))
          break;
      }
      if (iv < noldnp || ndep != 1)
        continue;
      if (t->dependsOnValue(excl))
        continue;
      ft.add(*t);
      if (!terms.find(*t))
        ftnew.add(*t);
    }
    if (oldnp1.getSize() > 0)
      continue; // all oldnp included in separate terms of the product
    if (!ft.getSize())
      continue;
    terms.add(ftnew);
    clients.add(*a);
    clientTerms.push_back(ft);
  }

  // Select the top-level ones of this type
  RooArgList dterms;
  for (Int_t i = 0, nterms = terms.getSize(); i < nterms; i++) {
    Int_t j;
    for (j = i + 1; j < nterms; j++)
      if (terms[j].dependsOnValue(terms[i]))
        break;
    if (j < nterms)
      dterms.add(terms[i]);
  }
  if (dterms.getSize() > 0) {
    if (verbose >= 1) {
      std::cout << "Remove dependent response terms:" << std::endl
                << Print(dterms, 2);
    }
    terms.remove(dterms);
    for (Int_t ic = clients.getSize() - 1; ic >= 0; ic--) {
      clientTerms[ic].remove(dterms, kTRUE);
      if (!clientTerms[ic].getSize()) {
        clients.remove(clients[ic]);
        clientTerms.erase(clientTerms.begin() + ic);
      }
    }
  }
  if (verbose >= 2) {
    std::cout << "Client terms:" << std::endl << Print(clients, 2);
  }
  if (verbose >= 1) {
    std::cout << "Response terms:" << std::endl << Print(terms, 2);
  }
  assert(clients.getSize() == (Int_t)clientTerms.size());

  int nrep = 0;
  for (Int_t ic = 0, nclients = clients.getSize(); ic < nclients; ic++) {
    const RooAbsArg *client = clients.at(ic);
    if (verbose >= 2)
      std::cout << "Client " << Print(*client) << ", client terms "
                << Print(clientTerms[ic]) << std::endl;
    RooAbsArg *fedit = 0;
    if (const RooProdPdf *term1 = dynamic_cast<const RooProdPdf *>(client)) {
      fedit = editTerm(*term1, clientTerms[ic], newterms,
                       newName(ws, client->GetName(), suffix));
    } else if (const RooProduct *term2 =
                   dynamic_cast<const RooProduct *>(client)) {
      fedit = editTerm(*term2, clientTerms[ic], newterms,
                       newName(ws, client->GetName(), suffix));
    }
    if (!fedit) {
      WARN("Unable to replace response terms %s in %s",
           Print(clientTerms[ic], 1).c_str(), Print(*client).c_str());
      continue;
    }
    if (verbose >= 1) {
      std::cout << "Replace " << Print(*client) << std::endl;
      std::cout << "with    " << Print(*fedit) << std::endl;
    }
    RooAbsArg *fimp = import(ws, *fedit);
    if (fimp) {
      if (RooFitUtils::renameSpec)
        *RooFitUtils::renameSpec =
            Form("%s=%s%s", fimp->GetName(), client->GetName(),
                 (RooFitUtils::renameSpec->empty() ? "" : ",")) +
            *RooFitUtils::renameSpec;
      cust.replaceArg(*client, *fimp);
      nrep++;
    }
    delete fedit;
  }

  return nrep;
}

int npsplit_response(RooWorkspace *ws, RooCustomizer &cust,
                     const RooAbsPdf &pdf, const RooArgList &oldnp,
                     const RooStats::HistFactory::FlexibleInterpVar *newterm_,
                     Bool_t copyCode = kFALSE, const char *select = 0,
                     int verbose = 1, const char *suffix = "_NEW") {
  const FlexibleInterpVarInspector &newterm =
      FlexibleInterpVarInspector::inspect(*newterm_);
  RooArgList newnp = newterm.paramList();

  int nrep = 0, nf = 0;
  Int_t np = oldnp.getSize(), nnewnp = newnp.getSize();
  const RooArgSet *comp = pdf.getComponents();
  if (select && *select) {
    const RooArgSet *sel =
        dynamic_cast<const RooArgSet *>(comp->selectByName(select));
    delete comp;
    comp = sel;
  }
  if (!comp->getSize()) {
    if (select && *select)
      WARN("no PDF terms matching %s", select);
    else
      WARN("no PDF terms");
    delete comp;
    return 0;
  }
  std::map<double, int> ndup;
  std::vector<std::map<std::pair<double, double>, int>> dup(np);
  for (RooFIter it = comp->fwdIterator(); RooAbsArg *a = it.next();) {
    const FlexibleInterpVarInspector *f =
        FlexibleInterpVarInspector::inspect(a);
    if (!f)
      continue;
    if (!isSubset(f->paramList(), oldnp))
      continue;
    Double_t nominal = f->nominal();
    ndup[nominal]++;
    for (Int_t i = 0; i < np; i++) {
      Int_t j = f->index(oldnp[i]);
      dup[i][std::make_pair(nominal - f->low(j), f->high(j) - nominal)]++;
    }
    RooStats::HistFactory::FlexibleInterpVar *fedit = f->edit(
        oldnp, newterm, 1.0, copyCode, newName(ws, f->GetName(), suffix));
    if (verbose >= 1) {
      std::cout << fedit->GetName();
      for (Int_t i = 0; i < np; i++) {
        const RooAbsArg *v = oldnp.at(i);
        std::cout << ' ' << v->GetName() << " = " << f->printParam(*v);
      }
      std::cout << " ->";
      const FlexibleInterpVarInspector &fn =
          FlexibleInterpVarInspector::inspect(*fedit);
      for (Int_t i = 0; i < nnewnp; i++) {
        const RooAbsArg *v = newnp.at(i);
        std::cout << ' ' << v->GetName() << " = " << fn.printParam(*v);
      }
      std::cout << std::endl;
    }
    nf++;
    RooAbsArg *fimp = import(ws, *fedit);
    if (fimp) {
      if (RooFitUtils::renameSpec)
        *RooFitUtils::renameSpec =
            Form("%s=%s%s", fimp->GetName(), f->GetName(),
                 (RooFitUtils::renameSpec->empty() ? "" : ",")) +
            *RooFitUtils::renameSpec;
      cust.replaceArg(*f, *fimp);
      nrep++;
    }
    delete fedit;
  }
  delete comp;

  if (!nf) {
    WARN("no matching FlexibleInterpVar terms for NPs %s",
         Print(oldnp).c_str());
  }
  if (ndup.size() > 1) {
    std::stringstream ss;
    for (std::map<double, int>::iterator it = ndup.begin(); it != ndup.end();
         ++it)
      ss << ' ' << it->first << " (x" << it->second << ")";
    WARN("different nominal values:%s", ss.str().c_str());
  }
  for (Int_t i = 0; i < np; i++) {
    size_t nd = dup[i].size();
    if (nd <= 1)
      continue;
    std::vector<std::pair<std::pair<double, double>, int>> vdup(dup[i].begin(),
                                                                dup[i].end());
    std::sort(vdup.begin(), vdup.end(), sortdup);
    std::stringstream ss;
    for (size_t j = 0; j < nd; j++) {
      double lo = vdup[j].first.first, hi = vdup[j].first.second;
      if (hi != lo)
        std::cout << " +" << hi << "-" << lo;
      else if (hi < 0.0)
        std::cout << " -/+" << fabs(hi);
      else
        std::cout << " +/-" << hi;
      ss << " (x" << vdup[j].second << ")";
    }
    WARN("%s different error ranges:%s", oldnp[i].GetName(), ss.str().c_str());
  }
  return nrep;
}

RooArgSet *np2constraint(const RooAbsPdf &pdf, const RooAbsArg &np,
                         const RooArgSet *allnp, const RooArgSet *observables,
                         int verbose = 1) {
  RooArgSet excl = *allnp;
  excl.remove(np, kTRUE);
  excl.add(*observables);
  // constraint term depends on the NP we'll replace, and on no observables or
  // other NPs.
  RooArgSet ptmp = np; // getAllConstraints updates its second arg
  RooArgSet *constraints = pdf.getAllConstraints(excl, ptmp, kFALSE);
  if (constraints->getSize() <= 0)
    WARN("no constraints for NP %s", Print(np).c_str());
  else if (constraints->getSize() > 1)
    WARN("multiple constraint terms for NP %s: %s", Print(np).c_str(),
         Print(*constraints).c_str());
  else if (verbose >= 1)
    std::cout << "Constraint term: " << Print(*constraints->first(), 2)
              << std::endl;
  return constraints;
}

int npsplit_constraint(RooWorkspace *ws, RooCustomizer &cust,
                       const RooAbsPdf &pdf, const RooArgList &oldconstraints,
                       const RooArgSet &newconstraints,
                       const char *selectClients = 0, int verbose = 1,
                       const char *suffix = "_NEW") {
  int nrep = 0;
  for (RooFIter it = oldconstraints.fwdIterator(); RooAbsArg *a = it.next();) {
    RooArgSet clients;
    TIterator *ia = a->clientIterator();
    while (TObject *o = ia->Next()) {
      const RooAbsArg *client = dynamic_cast<RooAbsArg *>(o);
      if (!pdf.dependsOn(*client))
        continue; // could be hanging around from previous edit
      if (client)
        clients.add(
            *client,
            kTRUE); // work round duplicates returned by clientIterator()
    }
    delete ia;
    if (selectClients && *selectClients) {
      TString sel = selectClients, suf = suffix;
      TPRegexp(",").Substitute(sel, suf + "*"); // also include recreated terms
      sel = selectClients + ("," + sel + suf + "*");
      RooAbsCollection *selClients = clients.selectByName(sel);
      if (verbose >= 2)
        std::cout << "Select " << Print(*a, 1) << " '" << sel
                  << "' clients: " << Print(clients) << " -> "
                  << Print(*selClients) << std::endl;
      clients.removeAll();
      clients.add(*selClients);
      delete selClients;
    }
    int nc = 0;
    for (RooFIter ic = clients.fwdIterator(); RooAbsArg *client = ic.next();) {
      if (!(client && pdf.dependsOnValue(*client)))
        continue;
      nc++;
      RooAbsArg *fedit = 0;
      if (const RooProdPdf *term1 = dynamic_cast<const RooProdPdf *>(client)) {
        fedit = editTerm(*term1, oldconstraints, newconstraints,
                         newName(ws, client->GetName(), suffix));
      } else if (const RooProduct *term2 =
                     dynamic_cast<const RooProduct *>(client)) {
        fedit = editTerm(*term2, oldconstraints, newconstraints,
                         newName(ws, client->GetName(), suffix));
      }
      if (!fedit) {
        WARN("Unable to replace constraint %s in %s", Print(*a, 1).c_str(),
             Print(*client).c_str());
        continue;
      }
      if (verbose >= 1) {
        std::cout << "Replace " << Print(*client) << std::endl;
        std::cout << "with    " << Print(*fedit) << std::endl;
      }
      RooAbsArg *fimp = import(ws, *fedit);
      if (fimp) {
        if (RooFitUtils::renameSpec)
          *RooFitUtils::renameSpec =
              Form("%s=%s%s", fimp->GetName(), client->GetName(),
                   (RooFitUtils::renameSpec->empty() ? "" : ",")) +
              *RooFitUtils::renameSpec;
        cust.replaceArg(*client, *fimp);
        nrep++;
      }
      delete fedit;
    }
    if (!nc) {
      WARN("no clients for constraint %s", Print(*a, 1).c_str());
    };
  }
  return nrep;
}

int editWorkspace(RooWorkspace *ws, RooCustomizer &cust, const RooAbsPdf *&pdf,
                  int nr, const char *termtype = "", int verbose = 1,
                  const char *orig = "orig") {
  // Edit workspace using RooCustomizer settings. Updates pdf pointer.
  TString suffix;
  for (int num = 1, nsel = -1; nsel; num++) {
    suffix = (num <= 1) ? orig : Form("%s%d", orig, num);
    RooAbsCollection *sel = ws->components().selectByName("*_" + suffix);
    nsel = sel->getSize();
    delete sel;
  }
  std::cout << "Replace " << nr << " " << termtype << " in " << Print(*pdf, 1)
            << ". Save original in " << pdf->GetName() << "_" << suffix << "."
            << std::endl;
  if (verbose >= 2)
    cust.Print("v");
  if (nr <= 0)
    return 0;
  const RooAbsArg *updpdf = cust.build(verbose >= 3);
  if (!updpdf) {
    ERROR("Could not build new PDF from %s", Print(*pdf, 1).c_str());
    return 0;
  }
  TString newname = updpdf->GetName();
  const RooAbsPdf *newpdf = 0;
  if (!ws->import(cust.cloneBranchList(), RooFit::Silence(),
                  RooFit::RenameConflictNodes(suffix, kTRUE),
                  RooFit::NoRecursion(kTRUE)))
    newpdf = ws->pdf(newname);
  if (!newpdf) {
    ERROR("Could not import new PDF %s into workspace '%s'",
          Print(*updpdf, 1).c_str(), ws->GetName());
    return 0;
  }
  pdf = newpdf;
  return nr;
}

// Factory interface
class CustIFace_npsplit : public RooFactoryWSTool::IFace {
public:
  CustIFace_npsplit(int verbose = 0)
      : _pdf(0), _cust(0), _verbose(verbose), _replaceCount(0),
        _queueEdits(false) {
    RooFactoryWSTool::registerSpecial("NP", this);
    RooFactoryWSTool::registerSpecial("FlexibleInterpVar", this);
    RooFactoryWSTool::registerSpecial("edit", this);
    RooFactoryWSTool::registerSpecial("ReplaceConstraint", this);
    RooFactoryWSTool::registerSpecial("ReplaceResponse", this);
    RooFactoryWSTool::registerSpecial("RemoveNP", this);
    RooFactoryWSTool::registerSpecial("RemoveDuplicates", this);
    RooFactoryWSTool::registerSpecial("QueueEdits", this);
    RooFactoryWSTool::registerSpecial("CommitEdits", this);
  }

  std::string create(RooFactoryWSTool &ft, const char *typeName,
                     const char *instanceName, std::vector<std::string> args) {
    if (_verbose >= 2) {
      std::cout << typeName << "::" << instanceName << "(" << Join(args, ",")
                << ")";
      if (RooFitUtils::currentMC)
        std::cout << " for '" << RooFitUtils::currentMC->GetName() << "'";
      std::cout << std::endl;
    }
    if (!RooFitUtils::renameSpec)
      RooFitUtils::renameSpec = new std::string;
    std::string type = typeName;
    if (type == "NP")
      return NP(ft, instanceName, args);
    if (type == "FlexibleInterpVar")
      return FlexibleInterpVar(ft, instanceName, args);
    if (type == "edit")
      return edit(ft, instanceName, args);
    if (type == "ReplaceConstraint")
      return ReplaceConstraint(ft, instanceName, args);
    if (type == "ReplaceResponse")
      return ReplaceResponse(ft, instanceName, args);
    if (type == "RemoveNP")
      return RemoveNP(ft, instanceName, args);
    if (type == "RemoveDuplicates")
      return RemoveDuplicates(ft, instanceName, args);
    if (type == "QueueEdits")
      return QueueEdits(ft, instanceName, args);
    if (type == "CommitEdits")
      return CommitEdits(ft, instanceName, args);
    throw std::string(
        Form("CustIFace_npsplit::create() ERROR: unknown type requested: %s",
             typeName));
  }

  std::string NP(RooFactoryWSTool &ft, const char *instanceName,
                 const std::vector<std::string> args);
  std::string FlexibleInterpVar(RooFactoryWSTool &ft, const char *instanceName,
                                const std::vector<std::string> args);
  std::string edit(RooFactoryWSTool &ft, const char *instanceName,
                   const std::vector<std::string> args);
  std::string ReplaceConstraint(RooFactoryWSTool &ft, const char *instanceName,
                                const std::vector<std::string> args);
  std::string ReplaceResponse(RooFactoryWSTool &ft, const char *instanceName,
                              const std::vector<std::string> args);
  std::string RemoveDuplicates(RooFactoryWSTool &ft, const char *instanceName,
                               const std::vector<std::string> args);
  std::string RemoveNP(RooFactoryWSTool &ft, const char *instanceName,
                       const std::vector<std::string> args);
  std::string QueueEdits(RooFactoryWSTool &ft, const char *instanceName,
                         const std::vector<std::string> args);
  std::string CommitEdits(RooFactoryWSTool &ft, const char *instanceName,
                          const std::vector<std::string> args);
  std::string CommitEdits(RooFactoryWSTool &ft,
                          const char *type = "queued edits");
  virtual ~CustIFace_npsplit() {
    if (_cust) {
      if (_replaceCount > 0)
        ERROR(
            "CustIFace_npsplit destructor ERROR: %d queued edits not committed",
            _replaceCount);
      delete _cust;
    }
  }

private:
  const RooAbsPdf *_pdf;
  RooCustomizer *_cust;
  int _verbose, _replaceCount;
  bool _queueEdits;
};

std::string CustIFace_npsplit::NP(RooFactoryWSTool &ft,
                                  const char *instanceName,
                                  const std::vector<std::string> args) {
  RooWorkspace *ws = &ft.ws();
  size_t nargs = args.size();
  std::string npName = instanceName;
  std::string pdfName =
      nargs >= 1 && args[0].size() ? args[0] : npName + "_Pdf";
  std::string gobName = nargs >= 2 && args[1].size() ? args[1] : npName + "_In";
  std::string sigName =
      nargs >= 3 && args[2].size() ? args[2] : std::string("1");
  if (nargs >= 4)
    throw std::string(
        Form("NuisanceParameter ERROR: too many arguments: %d", (int)nargs));

  RooAbsArg *np =
      import(ws, RooRealVar(instanceName, instanceName, 0.0, -5.0, 5.0));
  if (!np)
    throw std::string(
        Form("NuisanceParameter ERROR: importing NP %s", instanceName));
  RooAbsArg *gob =
      import(ws, RooRealVar(gobName.c_str(), gobName.c_str(), 0.0));
  if (!gob)
    throw std::string(
        Form("NuisanceParameter ERROR: importing global observable %s",
             gobName.c_str()));
  if (!import(ws, ft.asARG(sigName.c_str()), kTRUE))
    throw std::string(
        Form("NuisanceParameter ERROR: importing sigma %s", sigName.c_str()));
  if (!import(ws, RooGaussian(
                      pdfName.c_str(), pdfName.c_str(),
                      *dynamic_cast<RooAbsReal *>(ws->arg(instanceName)),
                      *dynamic_cast<RooAbsReal *>(ws->arg(gobName.c_str())),
                      *dynamic_cast<RooAbsReal *>(ws->arg(sigName.c_str())))))
    throw std::string(
        Form("NuisanceParameter ERROR: importing PDF %s", pdfName.c_str()));
  if (RooFitUtils::currentMC) {
    RooFitUtils::currentMC->SetNuisanceParameters(
        RooArgSet(*RooFitUtils::currentMC->GetNuisanceParameters(), *np));
    RooFitUtils::currentMC->SetGlobalObservables(
        RooArgSet(*RooFitUtils::currentMC->GetGlobalObservables(), *gob));
  }
  return instanceName;
}

std::string
CustIFace_npsplit::FlexibleInterpVar(RooFactoryWSTool &ft,
                                     const char *instanceName,
                                     const std::vector<std::string> args) {
  size_t nargs = args.size();
  RooArgList np;
  std::vector<double> vx;
  for (size_t i = 0; i < nargs; i++) {
    double x;
    if (Strtod(args[i], x))
      vx.push_back(x);
    else if (vx.size())
      throw std::string(
          Form("FlexibleInterpVar ERROR: bad number '%s'", args[i].c_str()));
    else {
      RooAbsArg *a = ft.ws().arg(args[i].c_str());
      if (!a)
        throw std::string(
            Form("FlexibleInterpVar ERROR: input RooAbsArg %s does not exist",
                 args[i].c_str()));
      np.add(*a);
    }
  }
  Int_t nnp = np.getSize(), nvx = vx.size();
  if (nvx < 2)
    throw std::string(Form(
        "FlexibleInterpVar ERROR: insufficient numeric arguments: %d", nvx));
  double nominal = vx[0];
  int code = vx[1] + 0.5;
  if (fabs(vx[1] - code) > 1e-6)
    throw std::string(
        Form("FlexibleInterpVar ERROR: bad code %s", args[nnp + 1].c_str()));
  std::vector<double> lo, hi;
  if (nnp + 2 == nvx) {
    for (Int_t i = 0; i < nnp; i++) {
      lo.push_back(nominal - vx[i + 2]);
      hi.push_back(nominal + vx[i + 2]);
    }
  } else if (2 * nnp + 2 == nvx) {
    for (Int_t i = 0; i < nnp; i++) {
      lo.push_back(nominal - vx[2 * i + 2]);
      hi.push_back(nominal + vx[2 * i + 3]);
    }
  } else {
    throw std::string(Form("FlexibleInterpVar ERROR: number of NPs (%d) and "
                           "numeric arguments (%d) don't match",
                           nnp, nvx));
  }
  RooAbsArg *f = new RooStats::HistFactory::FlexibleInterpVar(
      instanceName, instanceName, np, nominal, lo, hi,
      std::vector<int>(nnp, code));
  RooAbsArg *fi = import(&ft.ws(), *f);
  return fi ? fi->GetName() : instanceName;
}

std::string CustIFace_npsplit::QueueEdits(RooFactoryWSTool & /*ft*/,
                                          const char * /*instanceName*/,
                                          const std::vector<std::string> args) {
  if (args.size() > 1 || (args.size() == 1 && args[0].size() > 0))
    throw std::string(Form("QueueEdits ERROR: does no take arguments, have %d",
                           (int)args.size()));
  _queueEdits = true;
  if (_verbose >= 1)
    std::cout << "PDF edits will be queued until CommitEdits command is issued"
              << std::endl;
  return "";
}

std::string CustIFace_npsplit::CommitEdits(RooFactoryWSTool &ft,
                                           const char *type) {
  if (!_cust || !_pdf)
    throw std::runtime_error("CommitEdits ERROR: no edits to commit");
  const RooAbsPdf *newpdf = 0;
  if (_replaceCount <= 0)
    WARN("in CommitEdits: no queued edits to commit");
  else {
    newpdf = _pdf;
    editWorkspace(&ft.ws(), *_cust, newpdf, _replaceCount, type, _verbose);
    if (RooFitUtils::currentMC && newpdf &&
        _pdf == RooFitUtils::currentMC->GetPdf() && newpdf != _pdf)
      RooFitUtils::currentMC->SetPdf(*newpdf);
  }
  delete _cust;
  _cust = 0;
  _replaceCount = 0;
  _pdf = 0;
  return newpdf ? newpdf->GetName() : "";
}

std::string
CustIFace_npsplit::CommitEdits(RooFactoryWSTool &ft, const char *instanceName,
                               const std::vector<std::string> args) {
  if (args.size() > 1 || (args.size() == 1 && args[0].size() > 0))
    throw std::string(Form("CommitEdits ERROR: does no take arguments, have %d",
                           (int)args.size()));
  if (_replaceCount > 0 && strcmp(_pdf->GetName(), instanceName) != 0 &&
      !TPRegexp("^gobj\\d+$").MatchB(instanceName))
    WARN("CommitEdits::%s WARNING: queued edits are for PDF '%s'", instanceName,
         _pdf->GetName());
  std::string res = CommitEdits(ft);
  _queueEdits = false;
  return res;
}

std::string CustIFace_npsplit::edit(RooFactoryWSTool &ft,
                                    const char *instanceName,
                                    const std::vector<std::string> args) {
  RooWorkspace *ws = &ft.ws();
  const RooAbsPdf *pdf = ws->pdf(instanceName);
  if (!pdf && RooFitUtils::currentMC &&
      TPRegexp("^gobj\\d+$").MatchB(instanceName))
    pdf = RooFitUtils::currentMC->GetPdf();
  if (!pdf)
    throw std::string(
        Form("edit ERROR: PDF '%s' does not exist", instanceName));
  if (!(args.size() > 0 && args[0].size() > 0))
    throw std::runtime_error("edit ERROR: need at least one argument");
  if (_pdf && pdf != _pdf)
    CommitEdits(ft);
  _pdf = pdf;

  if (!_cust)
    _cust = new RooCustomizer(*pdf, 0);
  for (size_t i = 0, nargs = args.size(); i < nargs; i++) {
    const std::string &arg = args[i];
    size_t ic = arg.find('=');
    if (ic == std::string::npos)
      throw std::string(
          Form("edit ERROR: unknown argument: %s, expect form orig=subst",
               arg.c_str()));
    std::string origName = arg.substr(0, ic);
    const char *substName = arg.c_str() + ic + 1;
    RooAbsArg *orig = ws->arg(origName.c_str());
    RooAbsArg *subst = ws->arg(substName);
    if (!orig)
      throw std::string(Form("edit ERROR: input RooAbsArg %s does not exist",
                             origName.c_str()));
    if (!subst)
      throw std::string(
          Form("edit ERROR: input RooAbsArg %s does not exist", substName));
    _cust->replaceArg(*orig, *subst);
    _replaceCount++;
  }
  if (_queueEdits)
    return "";
  else
    return CommitEdits(ft, "edits");
}

std::string
CustIFace_npsplit::ReplaceConstraint(RooFactoryWSTool &ft,
                                     const char *instanceName,
                                     const std::vector<std::string> args) {
  // ReplaceConstraint (oldNP, newNP1, newNP2,...)
  // The new NPs can be specified as the NP variable, or (less ambiguously) the
  // constraint PDF.
  if (!RooFitUtils::currentMC)
    throw std::runtime_error("ReplaceConstraint ERROR: ModelConfig not set");
  RooWorkspace *ws = &ft.ws();
  const RooAbsPdf *pdf = ws->pdf(instanceName);
  if (!pdf && TPRegexp("^gobj\\d+$").MatchB(instanceName))
    pdf = RooFitUtils::currentMC->GetPdf();
  if (!pdf)
    throw std::string(
        Form("ReplaceConstraint ERROR: PDF '%s' does not exist", instanceName));
  if (_pdf && pdf != _pdf)
    CommitEdits(ft, _pdf->GetName());
  _pdf = pdf;

  int nargs = args.size(), iarg = nargs - 1;
  std::string select;
  if (nargs > 0) {
    size_t l = args[iarg].size();
    if (l >= 2 && ((args[iarg][0] == '{' && args[iarg][l - 1] == '}') ||
                   (args[iarg][0] == '\'' && args[iarg][l - 1] == '\''))) {
      select = args[iarg].substr(1, l - 2);
      nargs--;
    }
  }
  RooArgSet allPdfs = ws->allPdfs();
  RooArgList oldconstraints, newconstraints;
  for (iarg = 0; iarg < nargs; iarg++) {
    RooAbsArg *a = ws->arg(args[iarg].c_str());
    if (!a)
      throw std::string(
          Form("ReplaceConstraint ERROR: input RooAbsArg %s does not exist",
               args[iarg].c_str()));
    if (dynamic_cast<RooAbsRealLValue *>(a)) {
      RooArgSet dep;
      if (iarg <= 0) { // constraint term depends on the NP we'll replace, and
                       // on no observables or other NPs.
        RooArgSet *c = np2constraint(
            *pdf, *a, RooFitUtils::currentMC->GetNuisanceParameters(),
            RooFitUtils::currentMC->GetObservables(), _verbose);
        dep.add(*c);
        delete c;
      } else {
        for (RooFIter it = allPdfs.fwdIterator(); RooAbsArg *p = it.next();) {
          std::string c = p->ClassName();
          if ((c == "RooGaussian" || c == "RooLognormal" || c == "RooGamma" ||
               c == "RooPoisson" || c == "RooBifurGauss") &&
              p->dependsOnValue(*a))
            dep.add(*p);
        }
      }
      if (dep.getSize() <= 0)
        throw std::string(
            Form("ReplaceConstraint ERROR: no constraint term for NP %s",
                 a->GetName()));
      if (dep.getSize() > 1)
        throw std::string(Form(
            "ReplaceConstraint ERROR: multiple constraint terms for NP %s: %s",
            a->GetName(), dep.contentsString().c_str()));
      if (_verbose >= 1)
        std::cout << "NP " << a->GetName() << " has constraint term "
                  << Print(*dep.first()) << std::endl;
      a = dep.first();
    }
    if (iarg <= 0)
      oldconstraints.add(*a);
    else
      newconstraints.add(*a);
  }
  if (!_cust)
    _cust = new RooCustomizer(*pdf, 0);
  _replaceCount += npsplit_constraint(ws, *_cust, *pdf, oldconstraints,
                                      newconstraints, select.c_str(), _verbose);
  if (_queueEdits)
    return "";
  else
    return CommitEdits(ft, "constraint terms");
}

std::string
CustIFace_npsplit::ReplaceResponse(RooFactoryWSTool &ft,
                                   const char *instanceName,
                                   const std::vector<std::string> args) {
  if (!RooFitUtils::currentMC)
    throw std::runtime_error("ReplaceResponse ERROR: ModelConfig not set");
  RooWorkspace *ws = &ft.ws();
  const RooAbsPdf *pdf = ws->pdf(instanceName);
  if (!pdf && TPRegexp("^gobj\\d+$").MatchB(instanceName))
    pdf = RooFitUtils::currentMC->GetPdf();
  if (!pdf)
    throw std::string(
        Form("ReplaceResponse ERROR: PDF '%s' does not exist", instanceName));
  if (args.size() < 1)
    throw std::string(
        Form("NuisanceParameter ERROR: need at least 1 arguments, have %d",
             (int)args.size()));
  if (_pdf && pdf != _pdf)
    CommitEdits(ft, _pdf->GetName());
  _pdf = pdf;

  RooArgList argList;
  int nargs = args.size(), iarg = nargs - 1;
  std::string select;
  bool checkWild = true;
  if (nargs > 0) {
    size_t l = args[iarg].size();
    if (l >= 2 && ((args[iarg][0] == '{' && args[iarg][l - 1] == '}') ||
                   (args[iarg][0] == '\'' && args[iarg][l - 1] == '\''))) {
      select = args[iarg].substr(1, l - 2);
      nargs--;
      checkWild = false;
    }
  }
  for (iarg = 0; iarg < nargs; iarg++) {
    // compatibility: also allow wildcards at the end of the list without {} or
    // ''
    if (checkWild && args[iarg].find_first_of("*?[]") != std::string::npos)
      break;
    RooAbsArg *a = ws->arg(args[iarg].c_str());
    if (!a)
      throw std::string(
          Form("ReplaceResponse ERROR: input RooAbsArg %s does not exist",
               args[iarg].c_str()));
    argList.add(*a);
  }
  if (iarg < 1)
    throw std::string(Form("NuisanceParameter ERROR: need at least 1 "
                           "non-wildcard arguments, have %d",
                           iarg));
  RooArgList oldnp = argList[0];
  argList.remove(argList[0]);
  if (checkWild)
    select = Join(std::vector<std::string>(args.begin() + iarg, args.end()));

  if (!_cust)
    _cust = new RooCustomizer(*pdf, 0);
  if (iarg <= 1) {
    RooStats::HistFactory::FlexibleInterpVar fiv("empty_template",
                                                 "empty_template");
    int nrep = npsplit_response(ws, *_cust, *pdf, oldnp, &fiv, kFALSE,
                                select.c_str(), _verbose);
    if (nrep)
      _replaceCount += nrep;
    else
      _replaceCount += npsplit_response(
          ws, *_cust, *pdf, oldnp, argList,
          RooFitUtils::currentMC->GetNuisanceParameters(),
          RooFitUtils::currentMC->GetObservables(), select.c_str(), _verbose);
  } else if (const RooStats::HistFactory::FlexibleInterpVar *fiv =
                 dynamic_cast<const RooStats::HistFactory::FlexibleInterpVar *>(
                     argList.at(0))) {
    if (argList.getSize() != 1)
      throw std::string(Form("NuisanceParameter ERROR: only specify one "
                             "FlexibleInterpVar replacement, have %d",
                             argList.getSize()));
    _replaceCount += npsplit_response(ws, *_cust, *pdf, oldnp, fiv, kTRUE,
                                      select.c_str(), _verbose);
  } else
    _replaceCount += npsplit_response(
        ws, *_cust, *pdf, oldnp, argList,
        RooFitUtils::currentMC->GetNuisanceParameters(),
        RooFitUtils::currentMC->GetObservables(), select.c_str(), _verbose);
  if (_queueEdits)
    return "";
  else
    return CommitEdits(ft, "response terms");
}

std::string CustIFace_npsplit::RemoveNP(RooFactoryWSTool &ft,
                                        const char *instanceName,
                                        const std::vector<std::string> args) {
  // RemoveNP(NP)
  // Equivalent to:
  //   ReplaceConstraint(NP)
  //   ReplaceResponse(NP)
  bool separateEdits = true; // why doesn't response change work when included
                             // in the same edit as the constraint term
  if (args.size() == 2 && args[1] == "ONE")
    separateEdits = false;
  else if (args.size() != 1 || args[0].size() == 0)
    throw std::string(
        Form("RemoveNP ERROR: need 1 argument, have %d", (int)args.size()));
  if (!RooFitUtils::currentMC)
    throw std::runtime_error("RemoveNP ERROR: ModelConfig not set");
  RooWorkspace *ws = &ft.ws();
  const RooAbsPdf *pdf = ws->pdf(instanceName);
  if (!pdf && TPRegexp("^gobj\\d+$").MatchB(instanceName))
    pdf = RooFitUtils::currentMC->GetPdf();
  if (!pdf)
    throw std::string(
        Form("RemoveNP ERROR: PDF '%s' does not exist", instanceName));
  if (_pdf && pdf != _pdf)
    CommitEdits(ft, _pdf->GetName());
  _pdf = pdf;

  RooAbsArg *np = ws->arg(args[0].c_str());
  if (!np)
    throw std::string(Form("RemoveNP ERROR: input RooAbsArg %s does not exist",
                           args[0].c_str()));
  RooArgSet *constraint =
      np2constraint(*pdf, *np, RooFitUtils::currentMC->GetNuisanceParameters(),
                    RooFitUtils::currentMC->GetObservables(), _verbose);

  if (!_cust)
    _cust = new RooCustomizer(*pdf, 0);
  if (constraint->getSize() > 0) {
    _replaceCount += npsplit_constraint(ws, *_cust, *pdf, *constraint,
                                        RooArgList(), 0, _verbose);
    delete constraint;
    if (!_queueEdits && separateEdits) {
      std::string pdfname = CommitEdits(ft, "constraint terms");
      const RooAbsArg *a = pdfname.size() ? ws->arg(pdfname.c_str()) : 0;
      if (a && _verbose) {
        std::cout << "-> " << Print(*a) << std::endl;
      }
      const RooAbsPdf *newpdf = dynamic_cast<const RooAbsPdf *>(a);
      if (newpdf)
        _pdf = pdf = newpdf;
    }
  }

  if (!_cust)
    _cust = new RooCustomizer(*pdf, 0);
  RooStats::HistFactory::FlexibleInterpVar fiv("empty_template",
                                               "empty_template");
  int nrep = npsplit_response(ws, *_cust, *pdf, RooArgSet(*np), &fiv, kFALSE, 0,
                              _verbose);
  if (nrep)
    _replaceCount += nrep;
  else
    _replaceCount +=
        npsplit_response(ws, *_cust, *pdf, RooArgSet(*np), RooArgList(),
                         RooFitUtils::currentMC->GetNuisanceParameters(),
                         RooFitUtils::currentMC->GetObservables(), 0, _verbose);
  if (_queueEdits)
    return "";
  else
    return CommitEdits(ft, separateEdits ? "response terms" : "NP terms");
}

std::string
CustIFace_npsplit::RemoveDuplicates(RooFactoryWSTool &ft,
                                    const char *instanceName,
                                    const std::vector<std::string> args) {
  if (args.size() > 1 || (args.size() == 1 && args[0].size() > 0))
    throw std::string(Form("CommitEdits ERROR: does no take arguments, have %d",
                           (int)args.size()));
  RooWorkspace *ws = &ft.ws();
  const RooAbsPdf *pdf = ws->pdf(instanceName);
  if (!pdf && RooFitUtils::currentMC &&
      TPRegexp("^gobj\\d+$").MatchB(instanceName))
    pdf = RooFitUtils::currentMC->GetPdf();
  if (!pdf)
    throw std::string(
        Form("RemoveDuplicates ERROR: PDF '%s' does not exist", instanceName));
  if (_pdf && pdf != _pdf)
    CommitEdits(ft, _pdf->GetName());
  _pdf = pdf;

  if (!_cust)
    _cust = new RooCustomizer(*pdf, 0);
  _replaceCount += flexmerge(*_cust, *pdf, 0, _verbose);
  if (_queueEdits)
    return "";
  else
    return CommitEdits(ft, "duplicate terms");
}

#if defined(__ACLIC__) || defined(__CINT__) || defined(__CLING__) ||           \
    defined(__MAKECINT__)
void npsplit(int verbose = 1) {
  new CustIFace_npsplit(verbose); // register interface with RooFactoryWSTool
}
#else
namespace {
static CustIFace_npsplit *dummy =
    new CustIFace_npsplit(1); // register interface with RooFactoryWSTool
}
#endif
