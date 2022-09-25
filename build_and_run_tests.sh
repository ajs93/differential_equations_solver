#!/bin/bash
DOCKER_IMAGE="numeric_calculation:test"

docker run -v `pwd`:/home/user/project -w /home/user/project "$DOCKER_IMAGE"