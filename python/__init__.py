try:
    # this dummy code is needed to make sure the library is loaded
    import ROOT
    ROOT.RooFitUtils.ExtendedModel
    for k in dir(ROOT.RooFitUtils):
        locals()[k] = getattr(ROOT.RooFitUtils,k)
except AttributeError or ImportError as err:
    raise ImportError("unable to import RooFitUtils: "+err)


