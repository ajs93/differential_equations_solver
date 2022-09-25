#!/bin/bash
DOCKER_IMAGE=`cat docker_image_name.txt`

docker run -v `pwd`:/home/user/project -w /home/user/project "$DOCKER_IMAGE"