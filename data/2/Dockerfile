# syntax=docker/dockerfile:1

FROM python:3.8-slim-buster
COPY . /app
WORKDIR /app
CMD python question2.1.py --reference=reference.csv.gz --input=input.tped --output=output.txt
 