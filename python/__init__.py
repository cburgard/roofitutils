try:
    # this dummy code is needed to make sure the library is loaded
    import ROOT
    ROOT.gSystem.Load("libRooFitUtils.so")
    ExtendedModel = ROOT.RooFitUtils.ExtendedModel
    ExtendedMinimizer = ROOT.RooFitUtils.ExtendedMinimizer
except ImportError as err:
    pass
