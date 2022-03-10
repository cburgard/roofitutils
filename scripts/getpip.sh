if ! command -v pip &> /dev/null; then
    rm -rf get-pip.py
    wget -nv https://bootstrap.pypa.io/pip/get-pip.py
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
   
