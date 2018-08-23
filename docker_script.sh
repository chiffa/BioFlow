#!/usr/bin/env bash

# To build container:
docker build -t "bioflow:dockerfile" .

# To run the container:
docker run -p biflow:dockerfile

# To build the container after the compose:
docker-compose build

# To run the container after the compose: (in the compose directory)
docker-compose up -d


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

