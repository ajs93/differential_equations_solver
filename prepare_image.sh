#!/bin/bash

DOCKER_IMAGE=`cat docker_image_name.txt`

docker build -t "$DOCKER_IMAGE" .