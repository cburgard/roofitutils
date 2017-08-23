#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <map>

#include "TFile.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TTree.h"

#include "Math/MinimizerOptions.h"

#include "RooWorkspace.h"
#include "RooNLLVar.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooSimultaneous.h"
#include "RooRealSumPdf.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooPlot.h"

#include "RooStats/ModelConfig.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

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


RooFitResult* minimize(RooAbsReal* fcn, RooArgSet* minosSet = NULL) {
    TStopwatch timer;
    timer.Start();

    int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
    RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
    if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    int strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
    int save_strat = strat;

    RooMinimizer minim(*fcn);
    minim.optimizeConst(2);
    minim.setStrategy(strat);
    minim.setPrintLevel(printLevel);
    minim.setProfile(1);

    int status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());

    // up the strategy
    if (status != 0 && status != 1 && strat < 2) {
        strat++;
        cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
        minim.setStrategy(strat);
        status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    }

    if (status != 0 && status != 1) {
        cout << "Fit failed with status " << status << endl;
    }

    if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(msglevel);
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(save_strat);

    string name = Form("fitresult_%s",fcn->GetName());
    string title = Form("fitresult_%s",fcn->GetName());

    RooFitResult* fitresult = minim.save(name.c_str(), title.c_str());

    if (minosSet != NULL) {
        minim.minos(*minosSet);
    }

    timer.Stop();
    timer.Print();

    return fitresult;
}


double runMuhatQuick(string inFileName,
    const char* mypoi,
    string outFolder = "test",
    string wsName = "combined",
    string mcName = "ModelConfig",
    string dataName = "combData",
    string conditionalSnapshot = "")
{
  TStopwatch timer;
  timer.Start();

  TFile f(inFileName.c_str());
  RooWorkspace* ws = (RooWorkspace*)(f.Get(wsName.c_str()));
  if (!ws) {
    cout << "ERROR::Workspace: " << wsName << " doesn't exist!" << endl;
    return 0;
  }

  // Activate binned likelihood calculation for binned models
  // RooFIter iter = ws->components().fwdIterator();
  // RooAbsArg* arg;
  // while ((arg = iter.next())) {
  //   if (arg->IsA() == RooRealSumPdf::Class()) {
  //     arg->setAttribute("BinnedLikelihood");
  //   }
  // }

  ModelConfig* mc = (ModelConfig*)(ws->obj(mcName.c_str()));
  if (!mc) {
    cout << "ERROR::ModelConfig: " << mcName << " doesn't exist!" << endl;
    return 0;
  }

  RooDataSet* data = (RooDataSet*)(ws->data(dataName.c_str()));
  if (!data) {
    cout << "ERROR::Dataset: " << dataName << " doesn't exist!" << endl;
    return 0;
  }

  if (conditionalSnapshot != "") ws->loadSnapshot(conditionalSnapshot.c_str());
  ws->loadSnapshot("nominalGlobs");

  TIterator* nitr = mc->GetParametersOfInterest()->createIterator();
  RooRealVar* var;
  nitr->Reset();
  while ((var = (RooRealVar*)(nitr->Next()))) {
    cout << "Setting " << var->GetName() << " constant" << endl;
    var->setConstant(1);
  }

  RooRealVar* firstPOI = (RooRealVar*)(mc->GetParametersOfInterest()->find(mypoi));
  if (!firstPOI) {
    cout << "Trying to get the POI directly from the workspace." << endl;
    firstPOI = (RooRealVar*)(ws->var(mypoi));
  }
	if(!firstPOI){
		cout << "unable to find poi " << mypoi << endl;
		exit(1);
	}
  firstPOI->setVal(1.0);
  firstPOI->setConstant(0);

  int strategy = 0;
  //RooNLLVar::SetIgnoreZeroEntries(1);//
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(strategy);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);

  cout << "creating nll " << endl;

  RooNLLVar* nll = (RooNLLVar*)mc->GetPdf()->createNLL(*data, RooFit::Constrain(*mc->GetNuisanceParameters()), RooFit::GlobalObservables(*mc->GetGlobalObservables()), RooFit::Offset(1), RooFit::NumCPU(1, RooFit::Hybrid));

  cout.precision(15);
  cout << "NLL value before minimisation: " << nll->getVal() << endl;

  RooFitResult* result = minimize(nll);
  nll->enableOffsetting(0);

  f.Close();

  stringstream filename;
  filename << "root-files/" << "fit_" << outFolder << ".root";
  const char* fname(filename.str().c_str());
  ensureDirectoryForFile(fname);
  TFile of(fname, "recreate");
  
  TTree* resultTree = new TTree("result", "result");
  resultTree->SetDirectory(0);
  resultTree->Branch("result", &result);
  resultTree->Fill();
  resultTree->ResetBranchAddresses();
  resultTree->Write("", TObject::kOverwrite);
  of.Close();

  result->Print("v");
  double minnll = nll->getVal();
  cout.precision(15);
  cout << "NLL value after minimisation: " << minnll << endl;
  timer.Print();

  cout << "POI '" << firstPOI->GetName() << "' value is " << firstPOI->getVal() << " +/-"  << firstPOI->getError() << endl;

  return firstPOI->getVal();
}


int main(int argc, const char* argv[]){
  if(argc < 5){
    std::cout << "usage: runMuhatQuick inFileName mypoi outFolder wsName [mcName] [dataName] [conditionalSnapshot]" << std::endl;
    return 1;
  }
  std::string inFileName( argv[1] );
  std::string mypoi( argv[2] );
  std::string outFolder( argv[3] );
  std::string wsName( argv[4] );
  std::string mcName( argc > 5 ? argv[5] : "ModelConfig" );
  std::string dataName( argc > 6 ? argv[6] : "combData" );
  std::string conditionalSnapshot( argc > 7 ? argv[7] : "" );
  
  runMuhatQuick(inFileName,mypoi.c_str(),outFolder,wsName,mcName,dataName,conditionalSnapshot);
  return 0;
}
