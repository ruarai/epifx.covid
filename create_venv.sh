python3 -m venv venv_dev_local

source venv_dev_local/bin/activate

python3 -m pip install -r requirements.txt

python3 -m pip install --editable local_pypfilt
python3 -m pip install --editable local_epifx