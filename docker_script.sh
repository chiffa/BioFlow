#!/usr/bin/env bash

# To build container:
docker build -t "bioflow" .

# To run the container:
docker run bioflow

# create the volumes to which the docker containers will bind on the disk:
mkdir -p $BIOFLOWHOME/input
mkdir -p $BIOFLOWHOME/source
mkdir -p $BIOFLOWHOME/.internal/docker-mongo/db-data
mkdir -p $BIOFLOWHOME/.internal/docker-neo4j/db-data

# To build the container after the compose:
docker-compose build

# To run the container after the compose: (in the compose directory)
docker-compose up -d

# to publish to Docker, first tag the branch of interest:
docker tag XXXXXXXX chiffa/bioflow:latest

# then push to dockerhub
sudo docker push chiffa/bioflow


# basically the idea is to create a volume that would be mounted for every new instance and
# filled with data for the relevant organisms. Since the data is unlikely to change, the instance
# of neo4j can actually kept on it.

# this is however a non-essential development and will be skipped for now.

docker volume create my-vol

docker run -d \
    -it \
    --name devtest \
    --mount source=myvol2, target=/app \
    nginx:latest

# acutally, with my version:
docker run -d \
  -it \
  --name devtest \
  -v myvol2:/app \
  nginx:latest

docker container stop devtest
docker container rm devtest
docker volume rm myvol2

