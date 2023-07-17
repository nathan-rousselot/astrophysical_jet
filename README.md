# Astrophysical Jet simulations using MPI-AMRVAC

This is a numerical simulation of an astrophysical jet. Those jets are outflows of ionised matter emitted from the axis of rotation of galaxies. Most notable example of astrophysical jets are Quasars.

## Prerequisites

`apt install git`

`apt install perl`

`apt install gfortran`

`apt install openmpi-bin libopenmpi-dev`

`git clone https://github.com/amrvac/amrvac.git`

Note: In this repository, we use MPI-AMRVAC version 3.1 which is the latest version at the time of writing this README.md file. At the time of writing this README.md file, the latest release is 3.0. To switch to 3.1, use `git checkout amrvac3.1` after cloning the repository.

## Contents

This repository contains the following files:

## Running The Simulation

To run the simulation, first you need to create an output directory (if not created already).

`mkdir output`

Then, you need to compile the code. Go to root directory of MPI-AMRVAC 3.1 and run

`.\setup.pl -d=2`

Then, go to the directory of this repository and run

`make`

Finally, run the simulation using

`mpirun -np 112 ./amrvac -i jet.par`

NOTE : The number of processors used in the simulation is 112. You can change (and probably want) to change it. For the provided parameter file, you need a lot of CPUs or a lot of time. All simulations in the paper were run on 112 cores from AMD EPYC 7713 CPU.
