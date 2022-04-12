FROM python:3.8-slim

#Install SRRS python package
ADD setup.py setup.cfg pyproject.toml requirements_test.txt MANIFEST.in SRRS/
ADD src SRRS/src/
RUN pip install ./SRRS

