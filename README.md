# Differential equations solver in C++

This module is intended to provide a base for differential equation solvers written in C++.

## How to use

### Build testing image

This module provides a dockerfile image to allow unit tests to be run. In order to do so, first prepare the docker image in your local environment:

```bash
./prepare_image.sh
```

Then simply build and run the tests:

```bash
./build_and_run_tests.sh
```

### Fetching includes to work properly in the local environment

To fetch the include files of the docker image to develop in the local environment, the user may use a helper script:

```bash
./copy_includes.sh
```

It will fetch all the necessary header files from the docker container into the local environment to the ".include" directory.

### Library API

In order to check how to use this library, check the unit tests located in the tests/ directory of this project.