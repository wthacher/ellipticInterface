#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef __BCSETUP_H_
#define __BCSETUP_H_

#include "EBStencil.H"

#include "NamespaceHeader.H" 

/*whats the big idea here:


*/
class BCSetUp
{
public:
    //input data
    Vector< RefCountedPtr<LevelData<IVSFAB<Real> > > > m_taggedFaces; //faces that need new fluxers
    RefCountedPtr<LevelData<IVSFAB<Real> > > m_taggedFluxCells; //cells that need a new operator
    RefCountedPtr<LevelData<NodeFArrayBox> > m_psiNodes; //has geom info
    RefCountedPtr<LevelData<FluxBox> > m_faceIntVals; //has face intersection info
    RefCountedPtr<LevelData<IVSFAB<Real> > > m_BCData; //[BC0, BC1, data0, data1]. Bc \in {0,1,2}={per, no slip, calving front}


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
    bool conserve;
    Real constMu;
    
   
    RefCountedPtr<LevelData<IVSFAB<Real> > > m_moments;
    
    RefCountedPtr<BoxLayoutData<IVSFAB<EBStencil> > > BCOpStencil; //contains stencils in tagged Cells and BC cells for VTO calc

    virtual ~BCSetUp(){}

    BCSetUp(){};

    

    BCSetUp(const Vector< RefCountedPtr<LevelData<IVSFAB<Real>> > >& a_taggedFaces,
                const RefCountedPtr<LevelData<IVSFAB<Real> > >& a_taggedFluxCells,
                const RefCountedPtr<LevelData<NodeFArrayBox> >& a_psiNodes,
                const RefCountedPtr<LevelData<IVSFAB<Real> > >& a_moments,
                const RefCountedPtr< LevelData<IVSFAB<Real> > >& a_BCData,
                int a_Q,
                Real a_dx,
                Real a_weight)
    {
        define(a_taggedFaces, a_taggedFluxCells, a_psiNodes, a_moments, a_BCData, a_Q, a_dx, a_weight);
    }

    void define(const Vector< RefCountedPtr<LevelData<IVSFAB<Real>> > >& a_taggedFaces,
                const RefCountedPtr<LevelData<IVSFAB<Real> > >& a_taggedFluxCells,
                const RefCountedPtr<LevelData<NodeFArrayBox> >& a_psiNodes,
                const RefCountedPtr<LevelData<IVSFAB<Real> > >& a_moments,
                const RefCountedPtr< LevelData<IVSFAB<Real> > >& a_BCData,
                int a_Q,
                Real a_dx,
                Real a_weight)
    {
        m_taggedFaces.push_back(a_taggedFaces[0]);
        m_taggedFaces.push_back(a_taggedFaces[1]);
        m_taggedFluxCells = a_taggedFluxCells;
        m_psiNodes = a_psiNodes;
        m_moments = a_moments;
        m_BCData = a_BCData;
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
        LayoutData<IntVectSet> BCFact(dbl);
        IVSFABFactory<Real> ivf(BCFact);
        IVSFABFactory<EBStencil> ivfBC(BCFact);
        
        //not sure why i have to do this
        RefCountedPtr<BoxLayoutData<IVSFAB<EBStencil> > > BCOpStencilPtr( new BoxLayoutData<IVSFAB<EBStencil>>(dbl,2,ivfBC) );
        BCOpStencil = BCOpStencilPtr;
    }

    //need to make a function to fill in RHS data for calving front BC...but we can do that later


    
    void getFluxAndFriction(const LevelData<FArrayBox>& cellBeta, const LevelData<FArrayBox>& cellEta,
    const LevelData<IVSFAB<Real>>& phiCut, const LevelData<FArrayBox>& phi); 

    void getCoefficients(const IntVect center, const FArrayBox& etaFab, const FArrayBox& betaFab,
                        LAPACKMatrix& betaCoef, LAPACKMatrix& etaCoef,
                        const DataIndex& dit);
    void getCoefficientsFull(const IntVect center, const FArrayBox& etaFab, const FArrayBox& betaFab,
                        LAPACKMatrix& betaCoef, LAPACKMatrix& etaCoef,
                        const DataIndex& dit);
    //use faceIntersections that I pass you
    void getFluxFaceCell(const LAPACKMatrix& WMpW, const FArrayBox& cellNodes, 
                         const FluxBox& faceIntValsFab, const LAPACKMatrix& etaCoef, 
                         Vector<Real>& fluxFaceCell);

    void getFluxEBFaceCell(const LAPACKMatrix& WMpW, const IntVect iv, const DataIndex& dit,
                            const LAPACKMatrix& etaCoef, Vector<Real>& fluxEBFaceCell);

    void getEBFluxMoments(const Vector<Real>& nxMoments, const Vector<Real>& nyMoments,
                            const LAPACKMatrix& etaCoef, LAPACKMatrix& MFxGrd,
                            LAPACKMatrix& MFyGrd, LAPACKMatrix& MFxFlt, LAPACKMatrix& MFyFlt );
    
    void getFrictionEBCell(const LAPACKMatrix& WMpW, IntVect iv, const DataIndex& dit,
                            const LAPACKMatrix& betaCoef, Vector<Real>& frictionEBCell);

    //get the coefficients for one stencil
    void getPinvMSten(const IntVect center, const LAPACKMatrix& betaCoef,
                      const LAPACKMatrix& etaCoef, const DataIndex& dit,
                        LAPACKMatrix& WMpW);

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
    void updateBetaAndEta(LevelData<FArrayBox>& beta, 
                    LevelData<FArrayBox>& eta,
                    const Real AFace,
                    const LevelData<FArrayBox>& Ccoef,
                    const LevelData<FArrayBox>& hCoef,
                    const LevelData<IVSFAB<Real>>& phiCut,
                    const LevelData<FArrayBox>& phi,
                    const LevelData<FArrayBox>& a_H,
                    bool constFriction, bool constMu, Real muCGrd, Real muCFlt, Real betaScale);
    

};

#include "NamespaceFooter.H"

#endif