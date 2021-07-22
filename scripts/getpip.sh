if ! command -v pip &> /dev/null; then
    PYTHONVER=$(python -c "import sys; print(str(sys.version_info.major)+'.'+str(sys.version_info.minor))")
    rm -rf get-pip.py
    wget -nv https://bootstrap.pypa.io/pip/$PYTHONVER/get-pip.py
    python get-pip.py --user   
fi

PYTHON_ENV=$(python -c "import sys; sys.stdout.write('1') if hasattr(sys, 'real_prefix') else sys.stdout.write('0')")

if $PYTHON_ENV; then
    USER="--user"
else
    USER=""
fi

python -c "import scipy"   || python -m pip install $USER scipy
python -c "import skimage" || python -m pip install $USER scikit-image
   
