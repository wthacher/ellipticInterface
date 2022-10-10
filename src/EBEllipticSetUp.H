#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef __EBELLIPTICSETUP_H_
#define __EBELLIPTICSETUP_H_

#include "EBStencil.H"


#include "NamespaceHeader.H" 

/*whats the big idea here:


*/
class EBEllipticSetUp
{
public:
    //input data
    Vector< RefCountedPtr< LayoutData<IntVectSet > > > m_taggedFaces; //faces that need new fluxes
    RefCountedPtr< LayoutData<IntVectSet > > m_taggedFluxCells; //cells that need a new operator
    RefCountedPtr<LevelData<NodeFArrayBox> > m_psiNodes; //has geom info
    RefCountedPtr<LevelData<FluxBox> > m_faceIntVals; //has face intersection info
    int m_Q;
    int n;
    int loVol; 
    int hiVol ;
    int area ;
    int nX ;
    int nY ;
    int r;
    int nStenCells; 
    Real m_dx;
    Real m_weight;
    bool useGradJump;
    bool factorEta;
    bool conserve;
    
   
    RefCountedPtr<LevelData<IVSFAB<Real> > > m_moments;
    RefCountedPtr<LevelData<IVSFAB<Real> > > m_jumps;
    
    RefCountedPtr<BoxLayoutData<IVSFAB<EBStencil> > > EBOpStencil; //contains stencils in tagged Cells and EB cells for Elliptic calc
    RefCountedPtr<BoxLayoutData<IVSFAB<EBStencil> > > EBOpStencilCut;

    virtual ~EBEllipticSetUp(){}

    EBEllipticSetUp(){};

    EBEllipticSetUp(const Vector< RefCountedPtr< LayoutData<IntVectSet > > > a_taggedFaces,
                const RefCountedPtr< LayoutData<IntVectSet > > a_taggedFluxCells,
                const RefCountedPtr<LevelData<NodeFArrayBox> >& a_psiNodes,
                const RefCountedPtr<LevelData<IVSFAB<Real> > >& a_moments,
                int a_Q,
                Real a_dx,
                Real a_weight)
    {
        useGradJump = true;
        define(a_taggedFaces, a_taggedFluxCells, a_psiNodes, a_moments, a_Q, a_dx, a_weight);
    }

    void define(const Vector< RefCountedPtr< LayoutData<IntVectSet > > > a_taggedFaces,
                const RefCountedPtr< LayoutData<IntVectSet > > a_taggedFluxCells,
                const RefCountedPtr<LevelData<NodeFArrayBox> >& a_psiNodes,
                const RefCountedPtr<LevelData<IVSFAB<Real> > >& a_moments,
                int a_Q,
                Real a_dx,
                Real a_weight)
    {
        m_taggedFaces.push_back(a_taggedFaces[0]);
        m_taggedFaces.push_back(a_taggedFaces[1]);
        m_taggedFluxCells = a_taggedFluxCells;
        m_psiNodes = a_psiNodes;
        m_moments = a_moments;
        m_Q = a_Q;
        n = ((m_Q+1)*(m_Q+2))/2;
        //accessors in m_moments. to get the moment you want do (loVol + ith moment
        loVol = 0;
        hiVol = 1*n;
        area = 2*n;
        nX = 3*n;
        nY = 4*n;
        m_dx = a_dx;
        m_weight = a_weight;
        //r = (*m_psiNodes).ghostVect()[0]; 
        r = m_Q;
        nStenCells = (2*r+1)*(2*r+1);

        //set up LevelDatas
        DisjointBoxLayout dbl((*m_psiNodes).disjointBoxLayout());
        LayoutData<IntVectSet> EBFact(dbl);
        IVSFABFactory<Real> ivf(EBFact);
        IVSFABFactory<EBStencil> ivfEB(EBFact);
        
        //not sure why i have to do this
        RefCountedPtr<BoxLayoutData<IVSFAB<EBStencil> > > EBOpStencilPtr( new BoxLayoutData<IVSFAB<EBStencil>>(dbl,2,ivfEB) );
        EBOpStencil = EBOpStencilPtr;

        RefCountedPtr<BoxLayoutData<IVSFAB<EBStencil> > > EBOpStencilCutPtr( new BoxLayoutData<IVSFAB<EBStencil>>(dbl,2,ivfEB) );
        EBOpStencilCut = EBOpStencilCutPtr;

        RefCountedPtr<LevelData<IVSFAB<Real> > > m_jumpsPtr (new LevelData<IVSFAB<Real> >(dbl,2,m_moments->ghostVect(), ivf) );
        m_jumps = m_jumpsPtr;

        // for(DataIterator dit(dbl);dit.ok();++dit)
        // {
        //     const IVSFAB<Real>& momentsFab = (*m_moments)[dit()];
        //     const IntVectSet& EBCellsFab = momentsFab.getIVS();
        //     IVSFAB<Real>& jumpsFab = (*m_jumps)[dit()];
        //     jumpsFab.define(EBCellsFab,4);
        //     jumpsFab.setVal(0);
        // }
    }

    
    
    void getOperator(const LevelData<FArrayBox>& cellBeta, const LevelData<FArrayBox>& cellEta); 

    void getCoefficients(const IntVect center, const FArrayBox& etaFab, const FArrayBox& betaFab,
                        LAPACKMatrix& betaCoef, LAPACKMatrix& etaCoef,
                        const DataIndex& dit);
    void getCoefficientsFull(const IntVect center, const FArrayBox& etaFab, const FArrayBox& betaFab,
                        LAPACKMatrix& betaCoef, LAPACKMatrix& etaCoef,
                        const DataIndex& dit);
    //use faceIntersections that I pass you
    void getFluxFaceCell(const LAPACKMatrix& WMpW, const FArrayBox& cellNodes, 
                         const FluxBox& faceIntValsFab, const LAPACKMatrix& etaCoef, 
                         Vector<Real>& fluxFaceCell,
                         const LAPACKMatrix& CJump,Vector<Real>& fluxFaceCellCJump);

    void getFluxEBFaceCell(const LAPACKMatrix& WMpW, const IntVect iv, const DataIndex& dit,
                            const LAPACKMatrix& etaCoef, Vector<Real>& fluxEBFaceCell,
                             const LAPACKMatrix& CJump, Vector<Real>& fluxEBFaceCellCJump);

    void getEBFluxMoments(const Vector<Real>& nxMoments, const Vector<Real>& nyMoments,
                            const LAPACKMatrix& etaCoef, LAPACKMatrix& MFGrd,LAPACKMatrix& MFFlt);
    
    void getFrictionEBCell(const LAPACKMatrix& WMpW, IntVect iv, const DataIndex& dit,
                            const LAPACKMatrix& betaCoef, Vector<Real>& frictionEBCell,
                            const LAPACKMatrix& CJump,Vector<Real>& frictionEBCellCJump);

    //get the coefficients for one stencil
    void getPinvMSten(const IntVect center, const LAPACKMatrix& betaCoef,
                const LAPACKMatrix& etaCoef, const DataIndex& dit,
                LAPACKMatrix& WMpW, LAPACKMatrix& CJump);

    void getPinvMStenFull(const IntVect center, const LAPACKMatrix& betaCoef,
                        const LAPACKMatrix& etaCoef, const DataIndex& dit,
                        LAPACKMatrix& WMpW);


    void getFluxFaceCellFull(const LAPACKMatrix& WMpW, const LAPACKMatrix& etaCoef, 
                             const FArrayBox& cellNodes, Vector<Real>& fluxFaceCell);

    void getFrictionCell(const LAPACKMatrix& WMpW, IntVect iv, const DataIndex& dit,
                            const LAPACKMatrix& betaCoef, Vector<Real>& frictionCell);

    void getRegFaceStencil(int s, int d, const FArrayBox& etaFab, 
                            const IntVect iv, const DataIndex& dit,Vector<Real>& regFaceFlux );

    //eta will have 2 comps - grd then floating
  
};

#include "NamespaceFooter.H"

#endif