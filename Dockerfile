FROM ubuntu:14.04

# create and set-up home directory:
RUN cd /home \
     && mkdir ank
     && cd ank

# Mount data volumes that contain pre-build neo4j data