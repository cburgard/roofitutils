if ! command -v pip &> /dev/null; then
    PYTHONVER=$(python -c "import sys; print(str(sys.version_info.major)+'.'+str(sys.version_info.minor))")
    rm -rf get-pip.py
    wget -nv https://bootstrap.pypa.io/pip/$PYTHONVER/get-pip.py
    python get-pip.py --user   
fi

USER="--user"
if [[ $(whoami) == "root" ]]; then
   USER=""
fi

python -c "import scipy"   || python -m pip install $USER scipy
python -c "import skimage" || python -m pip install $USER scikit-image
   
