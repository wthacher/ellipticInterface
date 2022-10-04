#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef __EBSOLVEUTILITIES_H_
#define __EBSOLVEUTILITIES_H_

#include "IntersectionUtils.H"

#include "NamespaceHeader.H"



void convergenceTest(const LevelData<FArrayBox>& solCoarse,
                     const LevelData<FArrayBox>& solFine,
                     Real dxC, Real dxF, int Q,
                     LAPACKMatrix& err);
                     
void makeRegularCellStencil(LAPACKMatrix& SBeta,
                            Vector <LAPACKMatrix>&,
                            int order);



#include "NamespaceFooter.H"

#endif
