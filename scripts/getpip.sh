if ! command -v pip &> /dev/null; then
    PYTHONVER=$(python -c "import sys; print(str(sys.version_info.major)+'.'+str(sys.version_info.minor))")
    python -c "import sys; from RooFitUtils.util import version_greater_equal ; exit(version_greater_equal(sys.argv[1],sys.argv[2]))" $PYTHONVER 3.7
    PYTHONVER_GE_37=$?
    rm -rf get-pip.py
    if [  ${PYTHONVER_GE_37} -eq 1 ]; then
	wget -nv https://bootstrap.pypa.io/pip/get-pip.py
    else
	wget -nv https://bootstrap.pypa.io/pip/$PYTHONVER/get-pip.py
    fi
    python get-pip.py --user   
fi

PYTHON_ENV=$(python -c "import sys; sys.stdout.write('1') if hasattr(sys, 'real_prefix') else sys.stdout.write('0')")

if [ $PYTHON_ENV -eq 1 ]; then
    USER="--user"
else
    USER=""
fi

python -c "import scipy"   &> /dev/null || python -m pip install --user scipy         || python -m pip install scipy
python -c "import skimage" &> /dev/null || python -m pip install --user scikit-image  || python -m pip install scikit-image
   
