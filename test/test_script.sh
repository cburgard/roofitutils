# fit likelihood and obtain correlation matrix at minimum of likelihood (observed and expected)
#python ../scripts/fit.py --input test_workspace.root --workspace Test --data obsData --output test_results/fit_obsData.txt --findSigma --hesse
#python ../scripts/fit.py --input test_workspace.root --workspace Test --data asimovData --output test_results/fit_asimovData.txt --findSigma --hesse

# perform one dimension likelihood scans (observed and expected)
#python ../scripts/fit.py --input test_workspace.root --workspace Test --data obsData --output test_results/nllscan_mu_obsData.txt --scan "mu" "50" "-0.3" "3" --no-findSigma --no-hesse
#python ../scripts/fit.py --input test_workspace.root --workspace Test --data asimovData --output test_results/nllscan_mu_asimovData.txt --scan "mu" "50" "-0.3" "3" --no-findSigma --no-hesse
#python ../scripts/fit.py --input test_workspace.root --workspace Test --data obsData --output test_results/nllscan_norm_bkg_obsData.txt --scan "norm_bkg" "50" "0.9" "1.1" --no-findSigma --no-hesse
#python ../scripts/fit.py --input test_workspace.root --workspace Test --data asimovData --output test_results/nllscan_norm_bkg_asimovData.txt --scan "norm_bkg" "50" "0.9" "1.1" --no-findSigma --no-hesse

# perform two dimension likelihood scan (observed and expected)
#python ../scripts/fit.py --input test_workspace.root --workspace Test --data obsData --output test_results/nllscan_mu-norm_bkg_obsData.txt --scan "mu" "100" "-1" "3.5" --scan "norm_bkg" "50" "0.8" "1.2" --no-findSigma --no-hesse  
#python ../scripts/fit.py --input test_workspace.root --workspace Test --data asimovData --output test_results/nllscan_mu-norm_bkg_asimovData.txt --scan "mu" "100" "-1" "3.5" --scan "norm_bkg" "50" "0.8" "1.2" --no-findSigma --no-hesse  

# plot the correlation matrix (observed and expected)

#python ../scripts/plotcorrelation.py --input "color=black" test_results/fit_obsData.txt --output test_plots/corrmat_obsData.tex
#python ../scripts/plotcorrelation.py --input "color=black" test_results/fit_asimovData.txt --output test_plots/corrmat_asimovData.tex

# plot the one dimension likelihood scans (observed and expected)
#python ../scripts/plotscan.py --input "color=black" test_results/nllscan_mu_obsData.txt --labels "\\mu_{s}" --poi "mu" --output test_plots/nllscan_mu_obsData.tex
#python ../scripts/plotscan.py --input "color=black" test_results/nllscan_mu_asimovData.txt --labels "\\mu_{S}" --poi "mu" --output test_plots/nllscan_mu_asimovData.tex
#python ../scripts/plotscan.py --input "color=black" test_results/nllscan_norm_bkg_obsData.txt --labels "\\mu_{B}" --poi "norm_bkg" --output test_plots/nllscan_norm_bkg_obsData.tex
#python ../scripts/plotscan.py --input "color=black" test_results/nllscan_norm_bkg_asimovData.txt --labels "\\mu_{B}" --poi "norm_bkg" --output test_plots/nllscan_norm_bkg_asimovData.tex

# plot two dimension likelihood scans (observed and expected)
#python ../scripts/plotscan.py --input "color=black" test_results/nllscan_mu-norm_bkg_asimovData.txt --labels "\\mu_{S}" "\\mu_{B}" --output test_plots/nllscan_mu-norm_bkg_asimovData.tex
#python ../scripts/plotscan.py --input "color=black" test_results/nllscan_mu-norm_bkg_obsData.txt --labels "\\mu_{S}" "\\mu_{B}" --output test_plots/nllscan_mu-norm_bkg_obsData.tex

