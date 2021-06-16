# syntax=docker/dockerfile:1

FROM python:3.8-slim-buster
COPY . /app
WORKDIR /app
CMD python question2.1.py --reference=data/reference.csv.gz --input=data/input.tped --output=data/output.txt