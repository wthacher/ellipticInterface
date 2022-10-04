#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef __INTERSECTIONUTILS_H_
#define __INTERSECTIONUTILS_H_


#include  <iostream>
#include "FArrayBox.H"
#include "FluxBox.H"
#include "NodeFArrayBox.H"
#include "REAL.H"
#include "Vector.H"
#include "IntVect.H"
#include "LevelData.H"
#include "IVSFAB.H"
#include "RefCountedPtr.H"

#include "BoxIterator.H"
#include "Vector.H"
#include <CH_HDF5.H>
#include "LevelData.H"
#include "Lapack.H"
#include "AMRIO.H"
#include "LAPACKMatrix.H"

#include <assert.h> 

#include <PetscCompGrid.H>
#include "PetscCompGridPois.H"

#include "Stencil.H"

#include "MomentPolygon.H"
#include "Romberg.H"

#include "NamespaceHeader.H"


//helper functions
void get_points(Vector<Real>& x, Vector<Real>& y, Vector<Real>& x_int, Vector<Real>& y_int, Real dx=1);

//for some reason this is not in AMRIO
void writeLevelIVSFAB(const LevelData<IVSFAB<Real>>* a_dataPtr);

void print_matrix(LAPACKMatrix& A, int digits = 10);

void writeIVSFAB(const IVSFAB<Real>* a_dataPtr);

void writeIVS(const IntVectSet* a_dataPtr);
  
void shiftMoments(Vector<Real>& moments,int Q, Real sx, Real sy);                                      

void writeLDFAB(const LevelData<FArrayBox>* dataPtr);

Real faceMom(int px,int py, int d, int s);

Real faceMom(int px,int py, int d, int s,Vector<Real> bds);

void writePETSC(Vec u, const LevelData<FArrayBox>& cellIDs , LevelData<FArrayBox>& cellU);

int pseudoInvert(LAPACKMatrix& A);


#include "NamespaceFooter.H"

#endif
