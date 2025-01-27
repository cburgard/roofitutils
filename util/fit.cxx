// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2016-03-17
// Description : Perform quick unconditional fits
#include <iomanip>
#include <stdlib.h>
#include <list>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <chrono>
#include <regex>

#include "TFile.h"
#include "TH1D.h"
#include "TTime.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH2.h"
#include "Math/MinimizerOptions.h"

#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/RooStatsUtils.h"
#include "RooNLLVar.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooRealSumPdf.h"
#include "RooGaussian.h"
#include "RooMinimizer.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooLognormal.h"
#include "RooGamma.h"
#include "RooPoisson.h"
#include "RooBifurGauss.h"

#include "RooFitUtils/ExtendedModel.h"
#include "RooFitUtils/ExtendedMinimizer.h"
#include "RooFitUtils/Log.h"
#include "RooFitUtils/Utils.h"

#ifdef ROOTCORE
#include <RootCore/Packages.h>
#endif

#ifdef ROOTCORE_PACKAGE_Asg_Boost
#include "boost/program_options.hpp"
#include "boost/program_options/cmdline.hpp"
#include "boost/program_options/options_description.hpp"
#include "boost/program_options/parsers.hpp"
#include "boost/program_options/variables_map.hpp"

using namespace std;
using namespace RooFit;
using namespace RooStats;
using namespace RooFitUtils;

#include <chrono>
#include "TROOT.h"
class MyTimer {
public:
  MyTimer() : beg_(clock_::now()) {}
  void reset() { beg_ = clock_::now(); }
  double elapsed() const {
    return std::chrono::duration_cast<second_>
      (clock_::now() - beg_).count(); }

private:
  typedef std::chrono::high_resolution_clock clock_;
  typedef std::chrono::duration<double, std::ratio<1> > second_;
  std::chrono::time_point<clock_> beg_;
};

// _____________________________________________________________________________
// Declarations of functions used in this file

// _____________________________________________________________________________
// Main routine
int main(int argc, char** argv)
{
  TTime thistime = gSystem->Now();

  // Model information
  string inFileName      = "path/to/workspace.root";
  string wsName          = "combined";
  string modelConfigName = "ModelConfig";
  string dataName        = "combData";
  string snapshot        = "nominalNuis";
  string folder          = "test";

  // Parameter settings
  string poiName         = "mu";
  string profileName     = "";
  string fixName         = "";
  int fixAllNP           = 0;
  bool makeParameterSnapshots = true;

  // Fit settings
  string minimizerType   = "Minuit2";
  string minimizerAlgo   = "Migrad";
  int defaultStrategy    = 0;
  int binnedLikelihood   = 1;
  int offsetting         = 1;
  int constOpt           = 2;
  double eps             = 1.0;
  int numCPU             = 1;
  double precision       = 0.001;
  bool setInitialError   = false;

  // Misc settings
  int fixCache           = 1;
  int fixMulti           = 1;
  int eigendecomposition = 0;
  string loglevel        = "DEBUG";

  // Bookkeeping

  using namespace boost;
  namespace po = boost::program_options;
  po::options_description desc( "Program options" );
  desc.add_options()
    ( "help           , h"                                                                         , "Print this help message")
    ( "input"         , po::value<string>( &inFileName )->default_value( inFileName )              , "File to run over." )
    ( "poi"           , po::value<string>( &poiName )->default_value( poiName )                    , "POIs to measure." )
    ( "snapshot"      , po::value<string>( &snapshot )->default_value( snapshot )                  , "Initial snapshot." )
    ( "folder"        , po::value<string>( &folder )->default_value( folder )                      , "Output folder." )
    ( "profile"       , po::value<string>( &profileName )->default_value( profileName )            , "Parameters to profile." )
    ( "fix"           , po::value<string>( &fixName )->default_value( fixName )                    , "Parameters to fix." )
    ( "workspace"     , po::value<string>( &wsName )->default_value( wsName )                      , "WS to grab." )
    ( "modelconfig"   , po::value<string>( &modelConfigName )->default_value( modelConfigName )    , "MC to load." )
    ( "data"          , po::value<string>( &dataName )->default_value( dataName )                  , "Data to use." )
    ( "minimizerType" , po::value<string>( &minimizerType )->default_value( minimizerType )        , "Minimizer type." )
    ( "minimizerAlgo" , po::value<string>( &minimizerAlgo )->default_value( minimizerAlgo )        , "Minimizer algorithm." )
    ( "strategy"      , po::value<int>( &defaultStrategy )->default_value( defaultStrategy )       , "Default strategy." )
    ( "numCPU"        , po::value<int>( &numCPU )->default_value( numCPU )                         , "Number of CPUs." )
    ( "binned"        , po::value<int>( &binnedLikelihood )->default_value( binnedLikelihood )     , "Binned likelihood." )
    ( "starfix"       , po::value<int>( &fixCache )->default_value( fixCache )                     , "Fix StarMomentMorph cache." )
    ( "multifix"      , po::value<int>( &fixMulti )->default_value( fixMulti )                     , "Fix MultiPdf level 2." )
    ( "precision"     , po::value<double>( &precision )->default_value( precision )                , "Precision for scan." )
    ( "eps"           , po::value<double>( &eps )->default_value( eps )                            , "Convergence criterium." )
    ( "eigen"         , po::value<int>( &eigendecomposition )->default_value( eigendecomposition ) , "Eigenvalues and vectors." )
    ( "offset"        , po::value<int>( &offsetting )->default_value( offsetting )                 , "Offset likelihood." )
    ( "optimize"      , po::value<int>( &constOpt )->default_value( constOpt )                     , "Optimize constant terms." )
    ( "loglevel"      , po::value<string>( &loglevel )->default_value( loglevel )                  , "POIs to use." )
    ( "fixAllNP"      , po::value<int>( &fixAllNP )->default_value( fixAllNP )                     , "Fix all NP." )
    ;

  po::variables_map vm0;

  try {
    po::store( po::command_line_parser( argc, argv ).options( desc ).run(), vm0 );
    po::notify( vm0 );
  }

  catch ( std::exception& ex ) {
    cerr << "Invalid options: " << ex.what() << endl;
    cout << "Use ./a.out --help to print the allowed program options" << endl;
    return -1;
  }

  catch ( ... ) {
    cerr << "Unidentified error parsing program options." << endl;
    return -1;
  }

  // if help, print help
  if ( vm0.count( "help" ) ) {
    cout << "Usage: ./a.out [PROGRAMOPTIONS]\n";
    cout << desc;
    return 0;
  }

  // DEBUG OUTPUT
  // - ERROR
  // - WARNING
  // - INFO
  // - DEBUG
  Log::ReportingLevel() = Log::FromString(loglevel);

  // Configuration of minimizer
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minimizerType.c_str(), minimizerAlgo.c_str());
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(defaultStrategy);
  if (loglevel == "DEBUG") {
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);
  } else {
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  }
  int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
  if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  // Parse options
  vector<string> parsed = parseString(poiName, ",");

  // Load the model
  ExtendedModel* model = new ExtendedModel("model", inFileName, wsName,
                                            modelConfigName, dataName, snapshot,
                                            binnedLikelihood, "pdf_", fixCache,
                                            fixMulti);

  RooWorkspace* ws = model->GetWorkspace();
  ModelConfig* mc = model->GetModelConfig();
  RooAbsPdf* pdf = model->GetPdf();
  RooAbsData* data = model->GetData();
  RooArgSet* nuis = model->GetNuisanceParameters();
  RooArgSet* globs = model->GetGlobalObservables();
  RooArgSet* pois = model->GetParametersOfInterest();
  RooArgSet* obs = model->GetObservables();

  system(("mkdir -vp root-files/" + folder).c_str());

  if (fixAllNP) {
    model->fixNuisanceParameters();
  }

  if (setInitialError) {
    model->setInitialErrors();
  }

  if (fixName != "") {
    model->fixNuisanceParameters(fixName);
  }

  model->fixParametersOfInterest();

















  // Collect POIs
  vector<RooRealVar*> scan_poi_vector;
  RooArgSet scan_poi_set;
  vector<double> poiVals;
  for (size_t i = 0; i < parsed.size(); i++) {
    TString thisName = parsed[i];
    TString range;
    TString boundary;
    int sign = 0;

    bool useRange    = kFALSE;
    bool useBoundary = kFALSE;

    // Get ranges
    if (thisName.Contains("[")) {
      assert(thisName.Contains("]"));
      TObjArray* thisNameArray = thisName.Tokenize("[");
      thisName = ((TObjString*)thisNameArray->At(0))->GetString();
      range = ((TObjString*)thisNameArray->At(1))->GetString();
      range.ReplaceAll("]","");
      assert(range.Contains(":"));
      useRange = kTRUE;
    }

    // Get sign
    if (thisName.Contains("+")) {
      thisName.ReplaceAll("+",">0");
    } else if (thisName.Contains("-")) {
      thisName.ReplaceAll("-","<0");
    }

    // Get boundaries
    if (thisName.Contains(">")) {
      TObjArray* thisNameArray = thisName.Tokenize(">");
      thisName = ((TObjString*)thisNameArray->At(0))->GetString();
      boundary = ((TObjString*)thisNameArray->At(1))->GetString();
      sign = +1;
      useBoundary = kTRUE;
    } else if (thisName.Contains("<")) {
      TObjArray* thisNameArray = thisName.Tokenize("<");
      thisName = ((TObjString*)thisNameArray->At(0))->GetString();
      boundary = ((TObjString*)thisNameArray->At(1))->GetString();
      sign = -1;
      useBoundary = kTRUE;
    }

    RooRealVar* thisPoi = (RooRealVar*)ws->var(thisName);
    if (!thisPoi) {
      Log(logERROR) << "POI: " << thisName << " doesn't exist!";
      exit(-1);
    }

    thisPoi->removeRange();

    if (useRange) {
      double origVal = thisPoi->getVal();
      TObjArray* rangeArray = range.Tokenize(":");
      TString s_lo = ((TObjString*)rangeArray->At(0))->GetString();
      TString s_hi = ((TObjString*)rangeArray->At(1))->GetString();
      double lo = atof(s_lo.Data());
      double hi = atof(s_hi.Data());
      thisPoi->setRange(lo, hi);
      if ((origVal < lo) || (origVal > hi)) {
        double newVal = (hi - lo) / 2;
        thisPoi->setVal(newVal);
        cout << "setting " << newVal << endl;
      }
    }

    if (useBoundary) {
      double tmpBoundary = atof(boundary.Data());
      double origVal = thisPoi->getVal();
      double forigVal = fabs(thisPoi->getVal());
      bool boundaryIsZero = AlmostEqualUlpsAndAbs(tmpBoundary, 0.0, 0.0001, 4);

      if (sign > 0) {
        thisPoi->setMin(tmpBoundary);
        if (origVal < tmpBoundary) {
          thisPoi->setVal(tmpBoundary);
        }
        if (boundaryIsZero && origVal < 0) {
          thisPoi->setVal(forigVal);
        }
      } else if (sign < 0) {
        thisPoi->setMax(tmpBoundary);
        if (origVal > tmpBoundary) {
          thisPoi->setVal(tmpBoundary);
        }
        if (boundaryIsZero && origVal > 0) {
          thisPoi->setVal(-forigVal);
        }
      }
    }
    thisPoi->setError(0.2);
    thisPoi->setConstant(0);
    double val = thisPoi->getVal();
    Log(logINFO) << "Getting POI " << thisPoi->GetName() << " and set value to " << val;

    scan_poi_vector.push_back(thisPoi);
    scan_poi_set.add(*thisPoi);
    poiVals.push_back(val);
  }







  model->profileParameters(profileName);



  if (makeParameterSnapshots) {
    // Save the snapshots of nominal parameters
    Log(logINFO) << "MakeSnapshots() Saving nominal snapshots.";
    ws->saveSnapshot("nominalGlobs", *mc->GetGlobalObservables());
    ws->saveSnapshot("nominalNuis", *mc->GetNuisanceParameters());
    ws->saveSnapshot("nominalPois", *mc->GetParametersOfInterest());
  }

  MyTimer timer;
  ExtendedMinimizer minimizer("minimizer", pdf, data);
  minimizer.minimize(Minimizer(minimizerType.c_str(), minimizerAlgo.c_str()), Strategy(defaultStrategy), ExtendedMinimizer::Eps(eps),
                     Constrain(*nuis), GlobalObservables(*globs),
                     NumCPU(numCPU, 3), Offset(offsetting), Optimize(constOpt),
                     ExtendedMinimizer::Scan(scan_poi_set), Precision(precision), Hesse(), Save());//stars around extended minimizer 
  double time = timer.elapsed();
  Log(logINFO) << "Fitting time: " << setprecision(9) << time << " seconds";

  double minNll = minimizer.GetMinNll();
  Log(logINFO) << "NLL after minimisation: " << setprecision(15) << minNll;

  // Save fit result
  RooFitResult* result = (minimizer.GetFitResult());
  if(result){
    stringstream filename;
    filename << "cov_" <<  inFileName;
    const char* fname = filename.str().c_str();
    ensureDirectoryForFile(fname);
    TFile fout(fname, "recreate");
    if(fout.IsOpen()){
      result->Write("",TObject::kOverwrite);
      fout.Close();
      Log(logINFO) << "Saved fitresult as " << fname;
    } else {
      Log(logERROR) << "Unable to open file " << fname;
    }
    if (makeParameterSnapshots) {
      TString filename_with_snapshot(inFileName);
      filename_with_snapshot.ReplaceAll(".root", "_snapshots.root");
      
      RooArgSet floatingParameters(result->floatParsFinal());
      ws->saveSnapshot("ucmles", floatingParameters);
      
      // Bring us back to nominal for exporting
      Log(logINFO) << "MakeSnapshots() Return to nominal parameter values.";
      ws->loadSnapshot("nominalNuis");
      ws->loadSnapshot("nominalGlobs");
      ws->loadSnapshot("nominalPois");
      
      ws->writeToFile(filename_with_snapshot.Data());
    }
  } else {
    Log(logERROR) << "fit failed, no fit result obtained!";
  }  
  
  PrintResourcesUsed(thistime);
}

#else

int main(){
  std::cout << "this executable requires compilation with Asg_Boost" << std::endl;
  std::cout << "you can acquire it by typing" << std::endl;
  std::cout << "  rc checkout_pkg atlasoff/AsgExternal/Asg_Boost/tags && rc find_packages" << std::endl;
  std::cout << "and recompiling this package" << std::endl;
}

#endif
