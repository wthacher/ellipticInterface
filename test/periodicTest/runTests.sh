#!/bin/bash

rm sout
mpirun -np 8 ./periodicTest2d.Darwin.64.mpic++.mpif90.DEBUG.MPI.PETSC.ex inputs test.nCells=32 test.convergence = false >> sout
mv pout.0 pout.32
mv time.table.0 time.table.32
mpirun -np 8 ./periodicTest2d.Darwin.64.mpic++.mpif90.DEBUG.MPI.PETSC.ex inputs test.nCells=64 test.convergence = false >> sout
mv pout.0 pout.64
mv time.table.0 time.table.64
mpirun -np 8 ./periodicTest2d.Darwin.64.mpic++.mpif90.DEBUG.MPI.PETSC.ex inputs test.nCells=128 test.convergence = false >> sout
mv pout.0 pout.128
mv time.table.0 time.table.128
mpirun -np 8 ./periodicTest2d.Darwin.64.mpic++.mpif90.DEBUG.MPI.PETSC.ex inputs test.nCells=256 test.convergence = false >> sout
mv pout.0 pout.256
mv time.table.0 time.table.256
mpirun -np 8 ./periodicTest2d.Darwin.64.mpic++.mpif90.DEBUG.MPI.PETSC.ex inputs test.nCells=512 test.convergence = false >> sout
mv pout.0 pout.512
mv time.table.0 time.table.512
mpirun -np 8 ./periodicTest2d.Darwin.64.mpic++.mpif90.DEBUG.MPI.PETSC.ex inputs test.nCells=1024 test.convergence = false test.noSolve = true >> sout
mv pout.0 pout.1024
mv time.table.0 time.table.1024

rm sout2
./periodicTest2d.Darwin.64.mpic++.mpif90.DEBUG.MPI.PETSC.ex inputs test.nCells=32 test.convergence = true test.nCellsFine = 1024 >> sout2
./periodicTest2d.Darwin.64.mpic++.mpif90.DEBUG.MPI.PETSC.ex inputs test.nCells=64 test.convergence = true test.nCellsFine = 1024 >> sout2
./periodicTest2d.Darwin.64.mpic++.mpif90.DEBUG.MPI.PETSC.ex inputs test.nCells=128 test.convergence = true test.nCellsFine = 1024 >> sout2
./periodicTest2d.Darwin.64.mpic++.mpif90.DEBUG.MPI.PETSC.ex inputs test.nCells=256 test.convergence = true test.nCellsFine = 1024 >> sout2
./periodicTest2d.Darwin.64.mpic++.mpif90.DEBUG.MPI.PETSC.ex inputs test.nCells=512 test.convergence = true test.nCellsFine = 1024 >> sout2

rm pout.0 pout.1 pout.2 pout.3 pout.4 pout.5 pout.6 pout.7
rm time.table.0 time.table.1 time.table.2 time.table.3 time.table.4 time.table.5 time.table.6 time.table.7
