PYTHONVER=$(python -c "import sys; print(str(sys.version_info.major)+'.'+str(sys.version_info.minor))")
rm -rf get-pip.py
wget -nv https://bootstrap.pypa.io/pip/$PYTHONVER/get-pip.py
python get-pip.py --user
python -m pip install --user scipy
python -m pip install --user scikit-image
