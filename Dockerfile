FROM python:3.8

#Install sprawl python package into the docker image
ADD setup.py setup.cfg pyproject.toml requirements_test.txt MANIFEST.in sprawl/
ADD src sprawl/src/
RUN pip install ./sprawl

