#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef __EBRECONSTRUCTION_H_
#define __EBRECONSTRUCTION_H_

    
//#include "slu_ddefs.h"

#include "IntersectionUtils.H"


#include "NamespaceHeader.H" 
//get valid nodal values, interpolate valid nodes to invalid ones using 2D polynomials
void initializePsiNodes(ImplicitFunction& psi,
                    Real xl, Real yl, Real dx,
                    LevelData<NodeFArrayBox>& psiNodes,
                    LayoutData<IntVectSet>& EBCells,
                    Vector<LayoutData<IntVectSet>* > faceIntersections);


void tagEBCells(const LevelData<NodeFArrayBox>& psiNodes,
                LayoutData<IntVectSet>& EBCells,
                Vector<LayoutData<IntVectSet>* > faceIntersections);


void findIntersections(ImplicitFunction& psi,
                Real xl, Real yl, Real dx,
                const LevelData<NodeFArrayBox>& psiNodes,
                LayoutData<IntVectSet>& EBCells,
                Vector<LayoutData<IntVectSet>* > faceIntersections,
                LevelData<FluxBox>& faceIntVals);


void getMoments(ImplicitFunction& psi,
                Real xl, Real yl, Real dx,
                const LevelData<NodeFArrayBox>& psiNodes,
                const LevelData<FluxBox>& faceIntVals,
                LevelData<IVSFAB<Real>>& moments,
                int Q);

void tagCellsAndFaces(const LevelData<NodeFArrayBox>& psiNodes,
                const LayoutData<IntVectSet>& EBCells,
                RefCountedPtr<LayoutData<IntVectSet> >& taggedFluxCellsLD,
                Vector< RefCountedPtr< LayoutData<IntVectSet> > >& taggedFacesLD,
                int order,
                int velOrder);

#include "NamespaceFooter.H"

#endif