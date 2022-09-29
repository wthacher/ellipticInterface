# ellipticInterface

This code solves, in two dimensions, the elliptic equation -\beta u + div (eta F(u)) = f, where eta, beta, and f may be disconintinuous across some internal interface. We employ a conservative cut cell, finite volume scheme which uses least-squares interpolation to acheive high-order accuracy.

## Dependencies:

You will need the following libraries:

1. [Chombo](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations)
2. [PETSC](https://petsc.org/release/download/)

The [Chombo Design Document](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations?preview=/73468344/169279490/chomboDesign.pdf) has the information you will need to build the library on your system. It is highly recommended, but not necesarry, to build with MPI.

## Running an example

Once you have built the Chombo library, you are most of the way there. In one of the directories in the 'test' folder, edit the makefile to relfect your own directory layout; the PETSC_DIR and CHOMBO_HOME variables need to point to those respective libraries. Type 'make all' - if all goes well, this should produce an executable. The inputs file allows to you control options for each example.
