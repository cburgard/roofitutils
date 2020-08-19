/*==============================================================================
$Id: EditWorkspaces.cxx 239521 2016-06-29 16:33:13Z cburgard $

Copy workspace to another file, with renaming and edits, removing disconnected
items.

Usage:
  root -b -q INPUT_FILE editws.cc+'(OUTPUT_FILE, EDIT_SPEC)'
Or:
  root -b -q editws.cc+'(OUTPUT_FILE, EDIT_SPEC, INPUT_FILE)'
Or called from within ROOT:
  editws (RooWorkspace*, OUTPUT_FILE, EDIT_SPEC);

Examples:
  root -b -q wsin.root editws.cc+'("wsout.root", "edit_params.cfg")'
  root -b -q wsin.root editws.cc+'("wsout.root",
"EDIT::OLDPDF(OLDPDF,relM=sum::(mH,-125)); oldvar=newvar,var*=Var*")'

EDIT_SPEC gives the names of one or more files (separated by commas) containing
workspace factory editing commands and rename specs.
EDIT_SPEC can also include the commands (separated by ";") directly.

See http://root.cern.ch/root/html/RooFactoryWSTool.html#RooFactoryWSTool:process
for
details of the standard RooFit editing commands.
See npsplit.cc for additional editing commands.

"OLDPDF" in the EDIT_SPEC can be used for a PDF name - the string is replaced by
the ModelConfig's PDF name.

In addition, the edit spec can include rename specs of the form OLDVAR=NEWVAR,
where OLDVAR can include wildcards (*, ?, or [a-z]), and the matched text
replaces *s in NEWVAR. Note that all the rename specs are performed together
(serially)
when the workspace is copied following the other edit commands.

Author: Tim Adye <T.J.Adye@rl.ac.uk>, 23 Oct 2014.

==============================================================================*/

#include "TClass.h"
#include "TDataMember.h"
#include "TFile.h"
#include "TIterator.h"
#include "TKey.h"
#include "TError.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TPRegexp.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TRealData.h"
#include <algorithm>
#include <fstream>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooMsgService.h"
#include "RooStats/ModelConfig.h"
#include "RooTObjWrap.h"
#include "RooWorkspace.h"

#if defined(__ACLIC__) || defined(__CINT__) || defined(__CLING__) ||           \
    defined(__MAKECINT__)
#include "RooStats/ModelConfig.h"
#include <string>
namespace RooFitUtils {
#ifndef ROOFITUTILS_DEFINED
RooStats::ModelConfig *currentMC = 0;
std::string *renameSpec = 0;
#else
extern RooStats::ModelConfig *currentMC;
extern std::string *renameSpec;
#endif
}
int editws(const std::vector<std::string> &lines, RooWorkspace *w,
           const char *mcname = 0);
int editws(const RooWorkspace *w, RooWorkspace &wout, const char *mcname = 0,
           const char *dataname = 0);
int editws(const char *edit, RooWorkspace *w, const char *mcname = 0);
int editws(RooWorkspace *w, TDirectory *fout, const char *edit = 0,
           const char *mcname = 0, const char *dataname = 0);
int editws(RooWorkspace *w, const char *outfile, const char *edit = 0,
           const char *mcname = 0, const char *dataname = 0);
int editws(TDirectory *f, TDirectory *fout, const char *edit = 0,
           const char *wsname = 0, const char *mcname = 0,
           const char *dataname = 0);
int editws(TDirectory *f, const char *outfile, const char *edit = 0,
           const char *wsname = 0, const char *mcname = 0,
           const char *dataname = 0);
int editws(const char *outfile, const char *edit = 0, const char *wsfile = 0,
           const char *wsname = 0, const char *mcname = 0,
           const char *dataname = 0);
int editws(RooWorkspace *w, RooWorkspace *wout, const char *edit,
           const char *mcname = 0, const char *dataname = 0);
int editws(RooWorkspace *w, RooWorkspace *wout,
           const std::vector<std::string> &edits, const char *mcname = 0,
           const char *dataname = 0);
#else
#include "RooFitUtils.h"
#include <RooFitUtils/EditWorkspaces.h>
#endif

#ifdef ROOFITUTILS_DEFINED
#undef INFO
#undef ERROR
#undef WARN
#else
#define ROOFITUTILS_DEFINED
#endif

#define INFO(...) Info("editws", __VA_ARGS__)
#define ERROR(...) Error("editws", __VA_ARGS__)
#define WARN(...) Warning("editws", __VA_ARGS__)

#ifndef __CINT__
namespace {
#endif
// local definitions encapsulated in anonymous namespace
// Maximum string length for RooFit import renaming.
// This can be increased to 64000 in HSG7-patched ROOT 5.34.32.
int RooFit_buflen = 10240;

static int nocopy = 0;

class Rename {
  std::vector<TString> _spec, _exact, _re, _sub;
  mutable std::vector<int> _nmod;

public:
  Rename(const char *p = 0) {
    if (!p || !p[0])
      return;
#ifdef __CINT__
#pragma optimize 0 // turn off CINT bytecode optimization, which doesn't work
                   // for the following (ROOT bug #75875 ?)
#endif
    for (;;) {
      const char *q = strchr(p, ',');
      if (!q)
        q = p + strlen(p);
      TString exact, re, sub, spec;
      if (wild2subs(TString(p, q - p), exact, re, sub, spec)) {
        _exact.push_back(exact);
        _re.push_back(re);
        _sub.push_back(sub);
        _spec.push_back(spec);
        _nmod.push_back(0);
      }
      if (!*q)
        break;
      p = q + 1;
    }
#ifdef __CINT__
#pragma optimize 4 // restore default CINT optimization
#endif
  }

  size_t size() const { return _re.size(); }

  int checkUsed() const {
    int n = 0;
    for (size_t i = 0, nre = _re.size(); i < nre; i++) {
      if (_nmod[i])
        continue;
      WARN("Nothing to rename %s", _spec[i].Data());
      n++;
    }
    return n;
  }

  int rename(TString &name, bool verbose = true, bool canRemove = false,
             bool track = true) const {
    // Rename Note that substitutions occur serially, so can undo several edits
    int nmod = 0;
    TString oldname = name;
    for (size_t i = 0, nre = _re.size(); i < nre; i++) {
      if (!_re[i].Length()) {
        if (name != _exact[i])
          continue;
        name = _sub[i];
      } else {
        if (!TPRegexp(_re[i]).Substitute(name, _sub[i]))
          continue;
      }
      if (!canRemove && !name.Length()) {
        name = oldname;
        return 0;
      }
      nmod++;
      if (track)
        _nmod[i]++;
    }
    if (verbose && nmod) {
      if (name.Length())
        INFO("Rename %s -> %s", oldname.Data(), name.Data());
      else
        INFO("Do not copy %s", oldname.Data());
    }
    return nmod;
  }

  TString renamed(const char *name, bool verbose = true, bool canRemove = false,
                  bool track = true) const {
    TString n = name;
    rename(n, verbose, canRemove, track);
    return n;
  }

  static bool wild2subs(const TString &s, TString &exact, TString &re,
                        TString &sub, TString &spec) {
    sub = "";
    Int_t ie = s.Index('=');
    if (ie == kNPOS || ie == 0 || s.Index("=", ie + 1) != kNPOS ||
        s.Index("$", ie + 1) != kNPOS) {
      ERROR("Bad replacement spec '%s'", s.Data());
      re = "";
      return false;
    }
    re = TString(s(0, ie)).Strip(TString::kBoth);
    const TString s2 = TString(s(ie + 1, s.Length())).Strip(TString::kBoth);
    spec = re + " -> " + s2;
    TPRegexp("[.^$\\\\()]").Substitute(re, "\\$&", "g");
    re = "^" + re + "$";
    Int_t n1 =
        TPRegexp("(?:\\*+|\\?+|\\[[^]]+\\])+").Substitute(re, "($&)", "g");
    if (!n1) {
      re = "";
      exact = TString(s(0, ie)).Strip(TString::kBoth);
      sub = s2;
      return true;
    }
    exact = "";
    INFO("Wildcard rename %s", spec.Data());
    re.ReplaceAll("*", ".*");
    re.ReplaceAll("?", ".");
    Int_t iw = 0, n2 = 0;
    for (;; n2++) {
      Int_t jw = s2.Index('*', iw);
      if (jw == kNPOS)
        break;
      sub += s2(iw, jw - iw);
      sub += Form("$%d", n2 + 1);
      iw = jw + 1;
    }
    sub += s2(iw, s2.Length());
    if (n2 > n1 || n2 > 10) {
      ERROR("Replacement spec wildcard mismatch '%s'", s.Data());
      return false;
    }
    if (gDebug > 0)
      INFO("Replace s/%s/%s/g", re.Data(), sub.Data());
    return true;
  }
};

const char *Scan(const char *a, std::vector<std::string> &v, char sep = ',',
                 bool append = false) {
  if (!append)
    v.clear();
  if (!*a)
    return a;
  for (;; a++) {
    const char *b = strchr(a, sep);
    if (!b)
      b = a + strlen(a);
    v.push_back(std::string(a, b - a));
    if (!*b)
      return b;
    a = b;
  }
}

TObject *GetDataMember(const TObject *o, const char *name,
                       const TClass *cls = 0) {
  // Get access to a private data member. Checks type if cls is given.
  if (!o)
    return 0;
  const TClass *c = o->IsA();
  TRealData *rd = c->GetRealData(name);
  if (!rd)
    return 0;
  TDataMember *dm = rd->GetDataMember();
  if (!dm)
    return 0;
  if (cls && strcmp(dm->GetTypeName(), cls->GetName()) != 0)
    return 0;
  if (dm->IsaPointer())
    return *((TObject **)(((char *)o) + rd->GetThisOffset()));
  else
    return (TObject *)(((char *)o) + rd->GetThisOffset());
}

std::string PrintArg(const RooAbsArg &a, Int_t level = 2) {
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
    std::cout << " = ";
    a.printValue(os);
  }
  if (level >= 4) {
    std::cout << " ";
    a.printExtras(os);
  }
  return os.str();
}

RooAbsArg *
import(RooWorkspace &ws, const RooAbsArg &a, const char *newname = 0,
       const RooCmdArg &arg1 = RooCmdArg(), const RooCmdArg &arg2 = RooCmdArg(),
       const RooCmdArg &arg3 = RooCmdArg(), const RooCmdArg &arg4 = RooCmdArg(),
       const RooCmdArg &arg5 = RooCmdArg(), const RooCmdArg &arg6 = RooCmdArg(),
       const RooCmdArg &arg7 = RooCmdArg(),
       const RooCmdArg &arg8 = RooCmdArg()) {
  TString name = newname ? newname : a.GetName();
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  RooMsgService::instance().setGlobalKillBelow(
      RooFit::PROGRESS); // import(RooAbsData), which may be called internally,
                         // doesn't respect RooFit::Silence()
  ws.import(a, RooFit::Silence(), arg1, arg2, arg3, arg4, arg5, arg6, arg7,
            arg8);
  RooMsgService::instance().setGlobalKillBelow(msglevel);
  RooAbsArg *b = ws.arg(name);
  if (!b)
    WARN("Imported %s::%s not found", a.ClassName(), name.Data());
  return b;
}

RooAbsData *
import(RooWorkspace &ws, RooAbsData &a, const char *newname = 0,
       const RooCmdArg &arg1 = RooCmdArg(), const RooCmdArg &arg2 = RooCmdArg(),
       const RooCmdArg &arg3 = RooCmdArg(), const RooCmdArg &arg4 = RooCmdArg(),
       const RooCmdArg &arg5 = RooCmdArg(), const RooCmdArg &arg6 = RooCmdArg(),
       const RooCmdArg &arg7 = RooCmdArg(), const RooCmdArg &arg8 = RooCmdArg(),
       const RooCmdArg &arg9 = RooCmdArg()) {
  TString name = newname ? newname : a.GetName();
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  RooMsgService::instance().setGlobalKillBelow(
      RooFit::PROGRESS); // import(RooAbsData), which may be called internally,
                         // doesn't respect RooFit::Silence()
  ws.import(a, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
  RooMsgService::instance().setGlobalKillBelow(msglevel);
  RooAbsData *b = ws.data(name);
  if (!b)
    WARN("Imported %s::%s not found", a.ClassName(), name.Data());
  return b;
}

size_t ReadFile(const char *filename, std::vector<std::string> &v,
                bool strip = true, char comment = '#', bool append = true) {
  if (!append)
    v.clear();
  std::ifstream f(filename);
  if (!f) {
    ERROR("Could not read %s", filename);
    return 0;
  }
  size_t n = 1; // offset by 1 so 0=fail
  while (f) {
    std::string s;
    std::getline(f, s);
    if (strip) {
      if (comment != '\0') {
        size_t ic = s.find('#');
        if (ic != std::string::npos)
          s.erase(ic);
      }
      size_t i1 = s.find_first_not_of(' ');
      if (i1 == std::string::npos)
        continue; // blank line
      size_t i2 = s.find_last_not_of(' ');
      if (i1 > 0 || i2 + 1 < s.size())
        s = s.substr(i1, i2 - i1 + 1);
    }
    v.push_back(s);
    n++;
  }
  return n;
}

RooArgSet copyset(RooWorkspace &w, const RooArgSet *in, const Rename &rename,
                  const char *listtype = "args", bool renameInput = false) {
  RooArgSet out;
  if (!in || in->getSize() <= 0)
    return out;
  INFO("Copy %d %s", in->getSize(), listtype);
  RooFIter it = in->fwdIterator();
  RooAbsArg *a;
  while ((a = it.next())) {
    TString newname = a->GetName();
    bool renamed = rename.rename(newname, false, false, false);
    RooAbsArg *b = w.arg(newname);
    if (!b) {
      WARN("%s::%s not copied to %s", a->ClassName(), a->GetName(), listtype);
      continue;
    }
    out.add(*b);
    // Should be OK to rename arg in the input set if it is our own snapshot
    if (renameInput && renamed)
      a->SetName(newname);
  }
  return out;
}

int importSnaphots(const RooWorkspace *w, RooWorkspace &wout,
                   const Rename &rename) {
  const RooLinkedList *snaps = dynamic_cast<const RooLinkedList *>(
      GetDataMember(w, "_snapshots", RooLinkedList::Class()));
  if (!snaps)
    return 1;
  RooLinkedListIter it = snaps->iterator();
  TObject *o;
  while ((o = it.Next())) {
    const RooArgSet *snap = dynamic_cast<const RooArgSet *>(o);
    if (!snap)
      continue;
    TString newname = rename.renamed(snap->GetName(), true, true);
    if (!newname.Length())
      continue; // rename to "" to remove
    RooArgSet *newsnap = dynamic_cast<RooArgSet *>(snap->snapshot());
    newsnap->setName(newname);
    RooArgSet orig = copyset(wout, newsnap, rename,
                             Form("snapshot '%s' args", newname.Data()), true);
    RooArgSet *origsnap = dynamic_cast<RooArgSet *>(orig.snapshot());
    wout.saveSnapshot(newname, *newsnap, true);
    orig = *origsnap;
    delete origsnap;
    delete newsnap;
  }
  return 0;
}

RooCmdArg renameVariable(const RooAbsCollection *components,
                         const Rename &rename, int *cont = 0,
                         std::vector<std::string> *newTitle = 0) {
  if (!rename.size())
    return RooCmdArg();
  RooFIter it = components->fwdIterator();
  RooAbsArg *a;
  TString from = ",", to = ",";
  bool docont = false;
  for (int ind = 0; (a = it.next()); ind++) {
    if (cont && ind < *cont)
      continue; // skip to where we left off
    TString oldname = a->GetName(), newname = oldname;
    if (!rename.rename(newname, false))
      continue;
    if (!newname.Length() || components->find(newname) ||
        to.Index("," + newname + ",") != kNPOS) {
      WARN("Cannot rename %s to %s - duplicate name", oldname.Data(),
           newname.Data());
      continue;
    }
    if (cont && ind > 0 &&
        (from.Length() + oldname.Length() - 1 > RooFit_buflen ||
         to.Length() + newname.Length() - 1 > RooFit_buflen)) {
      *cont = ind;
      docont = true;
      break;
    }
    INFO("Rename %s::%s -> %s", a->ClassName(), oldname.Data(), newname.Data());
    from += oldname;
    to += newname;
    from += ",";
    to += ",";
    if (newTitle && a->GetTitle() == oldname)
      newTitle->push_back(newname.Data());
  }
  if (cont && !docont)
    *cont = 0;
  if (from.Length() <= 2)
    return RooCmdArg();
  from = from(1, from.Length() - 2);
  to = to(1, to.Length() - 2);
  if (from.Length() > RooFit_buflen || to.Length() > RooFit_buflen)
    WARN("RooFit::RenameVariable(from,to) string lengths (%d,%d) may be too "
         "long for RooWorkspace::import()",
         from.Length(), to.Length());
  return RooFit::RenameVariable(from, to);
}

Bool_t importPdf(RooWorkspace &w, RooAbsPdf *pdf, const Rename &rename,
                 TString &newname, const char *pdftype = "PDF") {
  if (!pdf) {
    newname = "";
    return false;
  }
  newname = pdf->GetName();
  RooCmdArg rentop;
  if (rename.rename(newname))
    rentop = RooFit::RenameVariable(pdf->GetName(),
                                    newname); // rename top-level pdf first
  if (w.pdf(newname)) {
    WARN("%s %s::%s already in workspace '%s'", pdftype, pdf->ClassName(),
         newname.Data(), w.GetName());
    return false;
  }
  Bool_t err = false;
  RooAbsCollection *components =
      RooArgSet(*pdf).snapshot(); // easy (if not so quick) way of getting names
                                  // of all dependents
  components->remove(*pdf, true, true);
  RooAbsArg *imppdf = pdf;
  RooWorkspace *savews = 0;
  std::vector<std::string> newTitle;
  int cont = 0;
  for (int itmp = 0;; itmp++) {
    RooCmdArg ren = renameVariable(components, rename, &cont, &newTitle);
    if (cont > 0) { // work round stupid limit on length of rename args
      RooWorkspace *tmpws = new RooWorkspace(w);
      tmpws->SetName(Form("%s_TEMP%d", w.GetName(), itmp));
      INFO("Import %s %s::%s to temporary workspace '%s' after testing %d/%d "
           "components",
           pdftype, pdf->ClassName(), newname.Data(), tmpws->GetName(), cont,
           components->getSize());
      imppdf = import(*tmpws, *imppdf, newname, RooFit::RecycleConflictNodes(),
                      rentop, ren);
      rentop = RooCmdArg();
      err = !imppdf;
      delete savews;
      savews = tmpws;
      if (err)
        break;
    } else {
      INFO("Import %s %s::%s to workspace '%s'", pdftype, pdf->ClassName(),
           newname.Data(), w.GetName());
      err = !import(w, *imppdf, newname, RooFit::RecycleConflictNodes(),
                    RooFit::Silence(), rentop, ren);
      break;
    }
  }
  delete savews;
  delete components;
  for (size_t i = 0, n = newTitle.size(); i < n; i++) {
    RooAbsArg *a = w.arg(newTitle[i].c_str());
    if (a)
      a->SetTitle(a->GetName());
  }
  return err;
}

int importModel(RooWorkspace &w, const RooStats::ModelConfig &mc,
                const Rename &rename, const char *mcname = 0) {
  int nerr = 0;
  TString newPdf, newPriorPdf;
  nerr += importPdf(w, mc.GetPdf(), rename, newPdf, "PDF");
  nerr += importPdf(w, mc.GetPriorPdf(), rename, newPriorPdf, "prior PDF");

  if (!mcname)
    mcname = mc.GetName(); // otherwise may have been stored by an alias
  RooStats::ModelConfig mcout(rename.renamed(mcname), &w);
  if (newPdf.Length())
    mcout.SetPdf(newPdf);
  if (newPriorPdf.Length())
    mcout.SetPriorPdf(newPriorPdf);
  const RooAbsData *pdata;
  if ((pdata = mc.GetProtoData()))
    mcout.SetProtoData(pdata->GetName());
  const RooArgSet *set;
  if ((set = mc.GetParametersOfInterest()))
    mcout.SetParametersOfInterest(copyset(w, set, rename, "POIs"));
  if ((set = mc.GetNuisanceParameters()))
    mcout.SetNuisanceParameters(copyset(w, set, rename, "NPs"));
  if ((set = mc.GetGlobalObservables()))
    mcout.SetGlobalObservables(copyset(w, set, rename, "global observables"));
  if ((set = mc.GetObservables()))
    mcout.SetObservables(copyset(w, set, rename, "observables"));
  if ((set = mc.GetConditionalObservables()))
    mcout.SetConditionalObservables(
        copyset(w, set, rename, "conditional observables"));
  if ((set = mc.GetConstraintParameters()))
    mcout.SetConstraintParameters(copyset(w, set, rename, "constraints"));
  INFO("Import model '%s' to workspace '%s'", mcout.GetName(), w.GetName());
  nerr += w.import(mcout);
  return nerr;
}

int importData(RooWorkspace &w, RooAbsData *data, const Rename &rename) {
  TString dname = data->GetName();
  RooCmdArg opt;
  if (rename.rename(dname, true, true))
    opt = RooFit::Rename(dname);
  if (!dname.Length())
    return 0;
  RooCmdArg ren = renameVariable(data->get(), rename);
  INFO("Import %s::%s to workspace '%s'", data->ClassName(), dname.Data(),
       w.GetName());
  return !import(w, *data, dname, opt, ren);
}

#ifndef __CINT__
} // end anonymous namespace
#endif

int editws(const RooWorkspace *w, RooWorkspace &wout, const char *mcname,
           const char *dataname) {
  // Copy workspace, all models (or just mcname), all datasets (or just
  // dataname), all snapshots,
  // and all generic objects. Everything is imported with renaming applied.
  // Does not include nodes that are not used by any of the models' pdfs.
  // Does not copy named sets, apart from those used by the models.
  Rename ren = (RooFitUtils::renameSpec ? RooFitUtils::renameSpec->c_str() : 0);
  delete RooFitUtils::renameSpec;
  RooFitUtils::renameSpec = 0;
  TString wsname = wout.GetName();
  if (ren.rename(wsname))
    wout.SetName(wsname);

  int nerr = 0;
  if (dataname && dataname[0]) {
    RooAbsData *data = w->data(dataname);
    if (data)
      nerr += importData(wout, data, ren);
  } else {
    std::list<RooAbsData *> dl = w->allData();
    while (!dl.empty()) {
      RooAbsData *data = dl.front();
      dl.pop_front();
      if (data)
        nerr += importData(wout, data, ren);
    }
  }

  if (mcname && mcname[0]) {
    RooStats::ModelConfig *mc =
        dynamic_cast<RooStats::ModelConfig *>(w->genobj(mcname));
    if (mc)
      nerr += importModel(wout, *mc, ren, mcname);
  } else {
    // Access _genObjects list directly (instead of via allGenericObjects()) to
    // see alias names
    const RooLinkedList *objs = dynamic_cast<const RooLinkedList *>(
        GetDataMember(w, "_genObjects", RooLinkedList::Class()));
    RooLinkedListIter ito = objs->iterator();
    TObject *obj;
    while ((obj = ito.Next())) {
      TString name = obj->GetName();
      if (obj->IsA() == RooTObjWrap::Class()) {
        RooTObjWrap *wobj = dynamic_cast<RooTObjWrap *>(obj);
        if (wobj)
          obj = wobj->obj();
      }
      if (obj->IsA() == RooStats::ModelConfig::Class()) {
        RooStats::ModelConfig *mc = dynamic_cast<RooStats::ModelConfig *>(obj);
        if (mc)
          nerr += importModel(wout, *mc, ren, name);
      } else {
        ren.rename(name);
        INFO("Import %s::%s to workspace '%s'", obj->ClassName(), name.Data(),
             wout.GetName());
        if (name == obj->GetName())
          nerr += wout.import(*obj);
        else if (obj->InheritsFrom(TNamed::Class())) {
          TObject *clone = obj->Clone(name);
          nerr += wout.import(*clone);
          delete clone;
        } else
          nerr += wout.import(*obj, name);
      }
    }
  }

  // Warn about unreferenced components.
  const RooArgSet &aout = wout.components();
  RooArgSet ain = w->components();
  ain.sort();
  RooFIter ita = ain.fwdIterator();
  RooAbsArg *a;
  TPRegexp re("_orig\\d*(?:_\\d+)?$");
  std::map<std::string, int> orig;
  while ((a = ita.next())) {
    if (!aout.find(ren.renamed(a->GetName(), false, false, false))) {
      TObjArray *m = re.MatchS(a->GetName());
      if (m->GetEntries() == 1)
        orig[((TObjString *)m->At(0))->GetString().Data()]++;
      else
        INFO("Did not copy %s::%s", a->ClassName(), a->GetName());
      delete m;
    }
  }
  for (std::map<std::string, int>::const_iterator i = orig.begin(),
                                                  e = orig.end();
       i != e; i++)
    INFO("Did not copy %4d '*%s' items", i->second, i->first.c_str());

  nerr += importSnaphots(w, wout, ren);
  nerr += ren.checkUsed();
  return nerr;
}

int editws(const std::vector<std::string> &lines, RooWorkspace *w,
           const char *mcname) {
  int nerr = 0;
  TObjArray allmc;
  if (mcname && mcname[0]) {
    TObject *obj = w->genobj(mcname);
    if (obj && obj->IsA() == RooStats::ModelConfig::Class())
      allmc.Add(obj);
  } else {
    std::list<TObject *> objs = w->allGenericObjects();
    while (!objs.empty()) {
      TObject *obj = objs.front();
      objs.pop_front();
      if (obj && obj->IsA() == RooStats::ModelConfig::Class())
        allmc.Add(obj);
    }
  }
  if (!allmc.GetEntries()) {
    WARN("No ModelConfig found in workspace '%s'", w->GetName());
    nerr++;
  }
  for (size_t il = 0, nl = lines.size(); il < nl; il++) {
    if (lines[il] == "NOCOPY") {
      nocopy++;
    }
    if (allmc.GetEntries() > 0 &&
        lines[il].find("OLDPDF") != std::string::npos) {
      TObjArray pdfNames(allmc.GetEntries()); // Use TNamed to map old->new PDF
                                              // nnames with name->title
      pdfNames.SetOwner();
      Int_t i, n = allmc.GetEntries();
      for (i = 0; i < n; i++) {
        RooStats::ModelConfig *mc =
            dynamic_cast<RooStats::ModelConfig *>(allmc.At(i));
        if (!mc)
          continue;
        RooAbsPdf *pdf = mc->GetPdf();
        if (pdf)
          pdfNames.AddAt(new TNamed(pdf->GetName(), ""), i);
      }
      for (i = 0; i < n; i++) {
        TNamed *pdfName = dynamic_cast<TNamed *>(pdfNames.At(i));
        if (!pdfName)
          continue;
        TString newPdfName, oldPdfName = pdfName->GetName();
        TNamed *done = dynamic_cast<TNamed *>(pdfNames.FindObject(oldPdfName));
        if (done == pdfName) { // didn't do this pdf before
          TString es = lines[il].c_str();
          es.ReplaceAll("OLDPDF", oldPdfName);
          std::cout << es << std::endl;
          RooFitUtils::currentMC =
              dynamic_cast<RooStats::ModelConfig *>(allmc.At(i));
          RooAbsArg *a = w->factory(es);
          RooFitUtils::currentMC = 0;
          if (!a)
            continue;
          std::cout << "-> " << PrintArg(*a) << std::endl;
          if (a->InheritsFrom(RooAbsPdf::Class()))
            newPdfName = a->GetName();
          else
            nerr++;
        } else
          newPdfName = done->GetTitle();
        pdfName->SetTitle(newPdfName);
        if (newPdfName.Length() > 0 && newPdfName != oldPdfName) {
          RooStats::ModelConfig *mc =
              dynamic_cast<RooStats::ModelConfig *>(allmc.At(i));
          mc->SetPdf(newPdfName);
          INFO("Set '%s' as PDF in model '%s'", newPdfName.Data(),
               mc->GetName());
        }
      }
    } else if (lines[il].find("=") != std::string::npos &&
               lines[il].find_first_of("({})") == std::string::npos) {
      if (!RooFitUtils::renameSpec)
        RooFitUtils::renameSpec = new std::string;
      if (RooFitUtils::renameSpec->size())
        *RooFitUtils::renameSpec += ",";
      *RooFitUtils::renameSpec += lines[il];
    } else {
      const char *e = lines[il].c_str();
      std::cout << e << std::endl;
      if (allmc.GetEntries() > 0)
        RooFitUtils::currentMC =
            dynamic_cast<RooStats::ModelConfig *>(allmc.At(0));
      RooAbsArg *a = w->factory(e);
      RooFitUtils::currentMC = 0;
      if (a) {
        std::cout << "-> " << PrintArg(*a) << std::endl;
      }
    }
  }
  return nerr;
}

int editws(const char *edit, RooWorkspace *w, const char *mcname) {
  // Perform RooFit factory edit of workspace with string or file.
  // Any lines containing the string "OLDPDF" will be run with the PDF name from
  // each ModelConfig.
  if (!edit || !edit[0])
    return 0;
  int nerr = 0;
  std::vector<std::string> lines;
  if (strpbrk("([{}]);=", edit)) {
    Scan(edit, lines, ';');
  } else {
    std::vector<std::string> files;
    Scan(edit, files);
    for (size_t jf = 0, nf = files.size(); jf < nf; jf++) {
      if (!ReadFile(files[jf].c_str(), lines))
        nerr++;
    }
  }
  return nerr + editws(lines, w, mcname);
}

int editws(RooWorkspace *w, TDirectory *fout, const char *edit,
           const char *mcname, const char *dataname) {
  int nerr = 0;
  if (edit)
    nerr += editws(edit, w, mcname);
  if (nocopy) {
    INFO("Write workspace '%s' to file %s", w->GetName(), fout->GetName());
    fout->cd();
    if (!w->Write())
      nerr++;
  } else {
    RooWorkspace wout(w->GetName(), w->GetTitle());
    nerr += editws(w, wout, mcname, dataname);
    INFO("Write workspace '%s' to file %s", wout.GetName(), fout->GetName());
    fout->cd();
    if (!wout.Write())
      nerr++;
  }
  gROOT->cd();
  return nerr;
}

int editws(RooWorkspace *w, RooWorkspace *wout, const char *edit,
           const char *mcname, const char *dataname) {
  int nerr = 0;
  if (!wout)
    return -1;
  if (edit)
    nerr += editws(edit, w, mcname);
  nerr += editws(w, *wout, mcname, dataname);
  return nerr;
}

int editws(RooWorkspace *w, RooWorkspace *wout,
           const std::vector<std::string> &edits, const char *mcname,
           const char *dataname) {
  int nerr = 0;
  if (!wout)
    return -1;
  if (edits.size() > 0)
    nerr += editws(edits, w, mcname);
  nerr += editws(w, *wout, mcname, dataname);
  return nerr;
}

int editws(RooWorkspace *w, const char *outfile, const char *edit,
           const char *mcname, const char *dataname) {
  TFile *fout = TFile::Open(outfile, "recreate");
  if (!fout)
    return 1;
  INFO("Copy%s workspace '%s' to file %s", ((edit && *edit) ? " and edit" : ""),
       w->GetName(), fout->GetName());
  int nerr = editws(w, fout, edit, mcname, dataname);
  delete fout;
  return nerr;
}

int editws(TDirectory *f, TDirectory *fout, const char *edit,
           const char *wsname, const char *mcname, const char *dataname) {
  INFO("Copy%s workspaces from file %s to file %s",
       ((edit && *edit) ? " and edit" : ""), f->GetName(), fout->GetName());
  int nws = 0, nerr = 0;
  if (wsname && wsname[0]) {
    RooWorkspace *w = 0;
    f->GetObject(wsname, w);
    if (w) {
      nerr += editws(w, fout, edit, mcname, dataname);
      nws++;
      delete w;
    }
  } else {
    TIter nextKey = f->GetListOfKeys();
    TKey *key;
    while ((key = dynamic_cast<TKey *>(nextKey()))) {
      const TClass *c = TClass::GetClass(key->GetClassName());
      if (c && c->InheritsFrom(RooWorkspace::Class())) {
        const TKey *ktop = dynamic_cast<TKey *>(
            f->GetListOfKeys()->FindObject(key->GetName()));
        if (ktop && key->GetCycle() != ktop->GetCycle())
          continue; // check this is the highest cycle
        RooWorkspace *w = dynamic_cast<RooWorkspace *>(key->ReadObj());
        if (w) {
          nerr += editws(w, fout, edit, mcname, dataname);
          nws++;
          delete w;
        }
      }
    }
  }
  if (!nws) {
    WARN("No workspaces found in file %s", f->GetName());
    nerr++;
  }
  return nerr;
}

int editws(TDirectory *f, const char *outfile, const char *edit,
           const char *wsname, const char *mcname, const char *dataname) {
  TFile *fout = TFile::Open(outfile, "recreate");
  if (!fout)
    return 1;
  int nerr = editws(f, fout, edit, wsname, mcname, dataname);
  delete fout;
  return nerr;
}

int editws(const char *outfile, const char *edit, const char *wsfile,
           const char *wsname, const char *mcname, const char *dataname) {
  TFile *f, *fopened = 0;
  if (wsfile && wsfile[0]) {
    f = fopened = TFile::Open(wsfile);
    gROOT->cd();
  } else {
    f = dynamic_cast<TFile *>(gROOT->GetListOfFiles()->At(0));
    if (!f)
      ERROR("No input workspace file specified");
  }
  if (!f)
    return 1;
  int nerr = editws(f, outfile, edit, wsname, mcname, dataname);
  delete fopened;
  return nerr;
}
