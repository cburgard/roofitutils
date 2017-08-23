// Performs a test workspace edit and print the edited workspace
{
  gROOT->ProcessLine (".x $ROOTCOREDIR/scripts/load_packages.C");
  editws ("editout.root", "RooFitUtils/test/atlas_hgamgam.cfg", "/afs/cern.ch/atlas/project/HSG7/Spring2014/HSG1/v14/Comb2011p2012_eEPSptt70_unblinded_v5.root");
  TFile* f= TFile::Open("editout.root");
  combination->Print();
  delete f;
}
