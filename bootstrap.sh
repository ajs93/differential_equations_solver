#!/bin/bash

BOOTSTRAP_FILE="bootstrap.sh"

if ! [[ -f "$BOOTSTRAP_FILE" ]]; then
  echo 'Please run this script from the EOL top path. Aborting...'
  return 1
fi

if [[ -z "$EOL_IMAGE_NAME" ]]; then
  echo 'Please export the EOL_IMAGE_NAME variable first with the desired docker image to be used. Aborting...'
  return 1
fi

if [[ -z "$FLASHING_PORT" ]]; then
  echo 'Please export the FLASHING_PORT variable first with the desired port to be used to flash the devices. Aborting...'
  return 1
fi

echo 'This script does a basic checking that you are in the top of the EOL project, please make sure that this is correct.'
echo 'Creating necessary environment variables...'

# Top path of the EOL project
export EOL_TOP_PATH="`pwd`"

# Auth mac for display needs
XAUTH_MAC=`xauth list 2>&1 | head -n 1`

XAUTH_MAC_F1=`echo $XAUTH_MAC | cut -d' ' -f1`
XAUTH_MAC_F2=`echo $XAUTH_MAC | cut -d' ' -f2`
XAUTH_MAC_F3=`echo $XAUTH_MAC | cut -d' ' -f3`

DISPLAY_NUMBER=`echo ${DISPLAY:1:1}`

FINAL_XAUTH_MAC="$XAUTH_MAC_F1$DISPLAY_NUMBER  $XAUTH_MAC_F2  $XAUTH_MAC_F3"

export XAUTH_MAC=$FINAL_XAUTH_MAC

# Docker GID to allow user in container to connect with host docker daemon
DOCKER_GID=`cut -d: -f3 < <(getent group docker)`

echo 'Environment ready!'

return 0