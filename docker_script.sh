#!/usr/bin/env bash

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

