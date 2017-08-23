// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2016-03-17
// Description : Common helper functions

#include <chrono>

#include "TROOT.h"
#include "TTime.h"
#include "TSystem.h"
#include "TMath.h"

#include "RooProdPdf.h"
#include "RooArgSet.h"
#include "RooWorkspace.h"

#include "RooFitUtils/Log.h"

#include "RooFitUtils/Utils.h"
namespace RooFitUtils {
  bool RooStarMomentMorphFix = true;
  bool RooMultiPdfFix = true;
}

namespace { 
  union MyFloat_t {
    MyFloat_t(float num = 0.0f) : f(num) {}
    // Portable extraction of components.
    bool Negative() const { return (i >> 31) != 0; }
    Int_t RawMantissa() const { return i & ((1 << 23) - 1); }
  Int_t RawExponent() const { return (i >> 23) & 0xFF; }
    
    Int_t i;
    float f;
  };
}

#include "definitions.h"
#ifdef HAS_ROOSTARMOMENTMORPH
#include "RooStarMomentMorph.h"
#endif
int fixRooStarMomentMorph(RooWorkspace* workspace){
  int retval = 0;
#ifdef HAS_ROOSTARMOMENTMORPH
  RooFIter iter(workspace->components().fwdIterator());
  RooAbsArg* arg;
  while ((arg = iter.next())) {
    RooStarMomentMorph* pdf = dynamic_cast<RooStarMomentMorph*>(arg);
    if (pdf){
      retval++;
      pdf->fixCache();
    }
  }
#endif
  return retval;
}


bool AlmostEqualUlpsAndAbs(float A, float B, float maxDiff, int maxUlpsDiff) {
  // Check if the numbers are really close -- needed  when comparing numbers
  // near zero.
  float absDiff = fabs(A - B);
  if (absDiff <= maxDiff) return true;

  MyFloat_t uA(A);
  MyFloat_t uB(B);

  // Different signs means they do not match.
  if (uA.Negative() != uB.Negative()) return false;

  // Find the difference in ULPs.
  int ulpsDiff = abs(uA.i - uB.i);
  if (ulpsDiff <= maxUlpsDiff) return true;

  return false;
}

// _____________________________________________________________________________
// Print used resources
// Courtesy of Tim Adye <T.J.Adye@rl.ac.uk>.
void PrintResourcesUsed(const TTime& progStart)
{
  ProcInfo_t info;
  if (gSystem->GetProcInfo(&info)<0) return;
  Long_t cput= TMath::CeilNint(info.fCpuUser);
  Long_t wall= Long64_t(gSystem->Now()-progStart+TTime(500))/Long64_t(1000);
  LOG(logINFO) << Form("resources used: cput=%02ld:%02ld:%02ld, mem=%ldkb, vmem=%ldkb, walltime=%02ld:%02ld:%02ld",
                       cput/3600, (cput/60)%60, cput%60,
                       info.fMemResident, info.fMemVirtual,
                       wall/3600, (wall/60)%60, wall%60);
}

// _____________________________________________________________________________
// Split strings according to separator
std::vector<std::string> parseString(const std::string& str, const std::string& sep)
{
  std::vector<std::string> parsed;
  int pos = 0;
  bool first = true;
  if (str.size() == 0) return parsed;
  if (str.find(sep) == std::string::npos) {
    parsed.push_back(str);
    return parsed;
  }

  while (true) {
    int newPos = str.find(sep, pos);
    if (str.find(sep, pos) == std::string::npos) {
      if (!first) parsed.push_back(str.substr(pos, newPos-pos));
      break;
    }

    std::string sub = str.substr(pos, newPos-pos);
    parsed.push_back(sub);
    pos = newPos+1;
    first = false;
  }

  return parsed;
}

// _____________________________________________________________________________
// Split a RooProdPdf into its components
void FindUniqueProdComponents( RooProdPdf* Pdf, RooArgSet& Components )
{
  static int counter = 0;
  counter++;

  if (counter > 50) {
    LOG(logERROR) << "FindUniqueProdComponents detected infinite loop. Please check.";
    exit(1);
  }

  RooArgList pdfList = Pdf->pdfList();
  if (pdfList.getSize() == 1) {
    LOG(logINFO) << "FindUniqueProdComponents " << pdfList.at(0)->GetName() << " is fundamental.";
    Components.add(pdfList);
  } else {
    TIterator* pdfItr = pdfList.createIterator();
    RooAbsArg* nextArg;
    while ((nextArg = (RooAbsArg*)pdfItr->Next())) {
      RooProdPdf* Pdf = (RooProdPdf*)nextArg;
      if (std::string(Pdf->ClassName()) != "RooProdPdf") {
        LOG(logINFO) << "FindUniqueProdComponents " << Pdf->GetName() << " is no RooProdPdf. Adding it.";
        Components.add(*Pdf);
        continue;
      }
      FindUniqueProdComponents(Pdf, Components);
    }
    delete pdfItr;
  }
  counter = 0;
}

bool ensureDirectory(const TString& path) {
  // ensure that the directory with the given path exists
  // check if directory <path> exists
  Long_t flags = 0;
  gSystem->GetPathInfo(path.Data(), (Long_t*)0, (Long_t*)0, &flags, (Long_t*)0);
  if (flags & 2) {
    // directory exists
    return true;
  } 
  //create directory
  if (0 == gSystem->mkdir(path.Data(),true)) return true;
  else return false;
}

bool ensureDirectoryForFile(const TString& file) {
  // ensure that the directory for the given file exists
  Ssiz_t pos = file.Last('/');
  if(pos == kNPOS) return false;
  return ensureDirectory(file(0,pos));
}
