#!/bin/sh

VENV_NAME="venv"

virtualenv -p python3 "${VENV_NAME}"

. "./${VENV_NAME}/bin/activate"

pip install -r requirements.txt
