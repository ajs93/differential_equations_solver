#!/bin/bash

DOCKER_IMAGE="numeric_calculation:test"
INCLUDE_PATH_DEST=".include"

rm -rf $INCLUDE_PATH_DEST
DOCKER_CONTAINER_ID=`docker run --name builder -d $DOCKER_IMAGE`
echo "Container ID: $DOCKER_CONTAINER_ID"
echo "Retrieving include files to path: $INCLUDE_PATH_DEST"
docker cp $DOCKER_CONTAINER_ID:/usr/include/. $INCLUDE_PATH_DEST
docker cp $DOCKER_CONTAINER_ID:/usr/local/include/. $INCLUDE_PATH_DEST
docker stop $DOCKER_CONTAINER_ID
docker container rm $DOCKER_CONTAINER_ID