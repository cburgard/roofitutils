if ! command -v pip &> /dev/null; then
    PYTHONVER=$(python -c "import sys; print(str(sys.version_info.major)+'.'+str(sys.version_info.minor))")
    rm -rf get-pip.py
    wget -nv https://bootstrap.pypa.io/pip/$PYTHONVER/get-pip.py
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
   
