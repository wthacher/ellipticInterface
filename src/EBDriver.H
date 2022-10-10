#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef __EBDRIVER_H_
#define __EBDRIVER_H_


#include "EBSolveUtilities.H"
#include "EBReconstruction.H"
#include "EBEllipticSetUp.H"
#ifdef CH_USE_PETSC
#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscviewer.h"
#endif 

#include "NamespaceHeader.H"

class EBDriver
{
public:
    Real m_weight ;
    Real m_dx ;
    int m_recOrder;
    int m_velOrder;
    bool m_verbose;
    bool useTagged; //use irreg cells?
    bool useGradJump;
    bool factorEta;
    bool HOCut;
    bool conserve;
    bool noEB;
    PetscInt id0, patchid0;
    LAPACKMatrix betaStencil;
    Vector<LAPACKMatrix> etaStencil;

    ImplicitFunction* m_IFPtr;
    RefCountedPtr<LevelData<IVSFAB<Real> > > m_moments; //i want to accesss this to do convergence tests
    RefCountedPtr<LevelData<NodeFArrayBox> > m_psiNodes; //if we want to pass an initial guess
    RefCountedPtr<EBEllipticSetUp> m_EBEllipticSU;
    RefCountedPtr< LevelData<FluxBox> > m_faceIntVals;
    RefCountedPtr<BoxLayoutData<IVSFAB<EBStencil> > > EBOpStencilCut;


    LevelData<FArrayBox> cellIDs;

    Vec  x, b;      /* approx solution, RHS */
    Mat  A;         /* linear system matrix */
    //KSP  ksp;

    int DOF;
    
    EBDriver()
    {
        m_weight = 3;
        m_dx  =1;
        m_recOrder = 2;
        m_verbose = true;
        m_velOrder = 2;
        useTagged = true;
        useGradJump = true;
        factorEta = true;
        HOCut = false;
        conserve = true;
        noEB = false;
        

    }

    void define()
    {
        
        RefCountedPtr<EBEllipticSetUp> m_EBEllipticSUPTr(new EBEllipticSetUp) ;
        m_EBEllipticSU = m_EBEllipticSUPTr;

        RefCountedPtr<LevelData<IVSFAB<Real> > > m_momentsPtr(new LevelData<IVSFAB<Real> > );
        m_moments = m_momentsPtr;
        RefCountedPtr<LevelData<NodeFArrayBox> > m_psiNodesPtr(new LevelData<NodeFArrayBox>);
        m_psiNodes = m_psiNodesPtr;

        LAPACKMatrix SBeta,empty;
        Vector<LAPACKMatrix> SEta(1, empty);
        makeRegularCellStencil(SBeta,SEta,m_velOrder);

        betaStencil = SBeta;
        etaStencil = SEta;


    }

    //uses the fast algorithm
    void reconstructEB(LevelData<FArrayBox>& H, 
                        ImplicitFunction& IF,
                         Real xl, Real yl, Real dx);
    //set up BC object
    int solver(const RefCountedPtr<LevelData<FArrayBox> >& a_cellBeta,
              const RefCountedPtr<LevelData<FArrayBox> >& a_cellEta,
              const RefCountedPtr<LevelData<IVSFAB<Real> > > uJumps,
              LevelData<FArrayBox>& a_phi,
              LevelData<IVSFAB<Real> >& a_phiCut,
              LevelData<FArrayBox>& a_rhs,
              LevelData<IVSFAB<Real> >& a_rhsCut);

    #ifdef CH_USE_PETSC

    // void MFPetscSolver(LevelData<FArrayBox>& a_phi,
    //                     const LevelData<FArrayBox>& a_rhs);
    
    PetscErrorCode labelCells(const LevelData<FArrayBox>& a_rhs,
                              const LevelData<IVSFAB<Real>>& a_rhsCut);

    PetscErrorCode assembleSystem(const LevelData<FArrayBox>& a_phi,
                        const LevelData<IVSFAB<Real>>& a_phiCut,
                        const LevelData<FArrayBox>& a_rhs,
                        const LevelData<IVSFAB<Real>>& a_rhsCut,
                        const LevelData<FArrayBox>& a_cellBeta,
                        const LevelData<FArrayBox>& a_cellEta);
    
    PetscErrorCode solveSystem(LevelData<FArrayBox>& a_phi,
                        LevelData<IVSFAB<Real>>& a_phiCut);
    #endif

    

private:


};


#include "NamespaceFooter.H"

#endif