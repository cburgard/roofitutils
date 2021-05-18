python scripts/fit.py -i test/test_workspace.root --data asimovData --poi mu --fit --write-workspace postfit.root
python scripts/fit.py -i postfit.root --poi mu --data asimovData --impacts 'norm_bkg' 'gamma_stat_*' --writeSubmit impacts 
bash impacts/jobs.txt
scripts/plotpulls.py -i red "impacts/*nominal*.json" --impacts mu "impacts/*.json" --range -2 2 --scaleimpacts 5 --numbers --atlas True --output impacts.tex
