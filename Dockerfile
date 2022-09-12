FROM python:3.8

#Install sprawl python package into the docker image
ADD package sprawl/
RUN pip install ./sprawl

