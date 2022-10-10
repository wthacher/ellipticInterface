#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#include "EBDriver.H"
#include "LevelDataOps.H"
#include "BiCGStabSolver.H"
#include "GMRESSolver.H"
#include "RelaxSolver.H"


#include "NamespaceHeader.H"


//need info about relationship between boxes and pts in space
void 
EBDriver::reconstructEB(LevelData<FArrayBox>& H,
                        ImplicitFunction& psi,
                        Real xl, Real yl, Real dx)
{
    CH_TIME("EBDriver::reconstructEB");
    m_IFPtr = &psi;
    //first set up the different objects that we need
    m_dx = dx;
    //see if we need extra order moments in cut cells
    int n = ((m_recOrder+1)*(m_recOrder+2))/2;
    int n1 = ((m_recOrder+2)*(m_recOrder+3))/2;
    int recOrder1 = m_recOrder+1;
    int velOrder1 = m_velOrder+1;
    if(!HOCut)
    {
        n1 = n;
        recOrder1 = m_recOrder;
        velOrder1 = m_velOrder;
    }

    DisjointBoxLayout dbl(H.disjointBoxLayout());
    DataIterator dit = dbl.dataIterator();
    IntVect ghost = H.ghostVect();
    
    LevelData<NodeFArrayBox>& psiNodes = *m_psiNodes;

    psiNodes.define(dbl, 1, ghost);

    LayoutData<IntVectSet> EBCells(dbl);

    Vector<LayoutData<IntVectSet>* > faceIntersections(2,NULL);
    for(int i =0;i<2;i++)
    {
        faceIntersections[i] = new LayoutData<IntVectSet> (dbl);
    }

    RefCountedPtr< LevelData<FluxBox> > faceIntValsPtr( new LevelData<FluxBox>(dbl,1,ghost));
    m_faceIntVals = faceIntValsPtr;

    if(m_verbose){pout() << "Reconstructing EB: " <<endl;}
    
    //this does the first step get in/out values and tag EBCells
    initializePsiNodes(psi, xl,yl,dx, psiNodes, EBCells, faceIntersections ) ;

    
    if(noEB)
    {
        for(dit.begin();dit.ok();++dit)
        {
            EBCells[dit()].define(); //clear this guy
            for(int d=0;d<2;d++)
            {
                LayoutData<IntVectSet>& faceLD = *faceIntersections[d];
                IntVectSet& taggedFaces = faceLD[dit()]; 
                taggedFaces.define();
            }
            
        }
    }
    
    //second step - using nodal info to find roots of b field
    
    findIntersections(psi, xl,yl,dx, psiNodes, EBCells,faceIntersections,*faceIntValsPtr);
    
    
    //set up outputs
    IVSFABFactory<Real> ivf(EBCells);
    RefCountedPtr<LevelData<IVSFAB<Real> > > momentsPtr(new LevelData<IVSFAB<Real>> (dbl,5*n1, ghost,ivf));
    m_moments = momentsPtr;
    LevelData<IVSFAB<Real>>& moments = *momentsPtr;
    RefCountedPtr< LevelData<IVSFAB<Real>> > EBCellsPtr(new LevelData<IVSFAB<Real>>(dbl, 1, ghost,ivf));
    
    //use roots to get moments using Greens thm stuff
    getMoments(psi, xl,yl,dx, psiNodes, *faceIntValsPtr,moments,recOrder1);

    
    //RefCountedPtr<LevelData<IVSFAB<Real> > > taggedFluxCells ( new LevelData<IVSFAB<Real>>(dbl, 1, ghost,ivf) );
    RefCountedPtr<LayoutData<IntVectSet>> taggedFluxCells(new LayoutData<IntVectSet>(dbl));
    
    // Vector< RefCountedPtr<LevelData<IVSFAB<Real > > > > taggedFaces(2);
    // RefCountedPtr<LevelData<IVSFAB<Real > > > taggedFacesX (new LevelData<IVSFAB<Real>>(dbl, 1, ghost,ivf));
    // RefCountedPtr<LevelData<IVSFAB<Real > > > taggedFacesY (new LevelData<IVSFAB<Real>>(dbl, 1, ghost,ivf));

    Vector< RefCountedPtr<LayoutData<IntVectSet>> > taggedFaces(2);
    RefCountedPtr<LayoutData<IntVectSet>> taggedFacesX (new LayoutData<IntVectSet>(dbl));
    RefCountedPtr<LayoutData<IntVectSet>> taggedFacesY (new LayoutData<IntVectSet>(dbl));

    taggedFaces[0] = taggedFacesX;
    taggedFaces[1] = taggedFacesY;
    
    if(m_verbose){pout() << "Setting up RHS: " <<endl;}
    tagCellsAndFaces(psiNodes, EBCells,  taggedFluxCells, taggedFaces, m_recOrder,m_velOrder);
    

    EBEllipticSetUp& EBEllipticSU = *m_EBEllipticSU;
    if(m_verbose){pout() << "Setting up EB Elliptic stencils: " << endl;}
    EBEllipticSU.m_faceIntVals = faceIntValsPtr;
    EBEllipticSU.useGradJump = useGradJump;
    EBEllipticSU.factorEta = factorEta;
    EBEllipticSU.conserve = conserve;


    EBEllipticSU.define(taggedFaces, taggedFluxCells, m_psiNodes, momentsPtr, velOrder1, dx, m_weight);  

    if(m_verbose){pout() << "Done. " << endl;} 
    
}


int
EBDriver::solver(const RefCountedPtr<LevelData<FArrayBox> >& a_cellBeta,
              const RefCountedPtr<LevelData<FArrayBox> >& a_cellEta,
              const RefCountedPtr<LevelData<IVSFAB<Real> > > uJumps,
              LevelData<FArrayBox>& a_phi,
              LevelData<IVSFAB<Real> >& a_phiCut,
              LevelData<FArrayBox>& a_rhs,
              LevelData<IVSFAB<Real> >& a_rhsCut)
{
    CH_TIME("EBDriver::solver");
    PetscErrorCode ierr;

    int n = ((m_recOrder+1)*(m_recOrder+2))/2;
    int n1 = ((m_recOrder+2)*(m_recOrder+3))/2;
    int recOrder1 = m_recOrder+1;
    int velOrder1 = m_velOrder+1;
    if(!HOCut)
    {
        n1 = n;
        recOrder1 = m_recOrder;
        velOrder1 = m_velOrder;
    }
    
    m_EBEllipticSU->m_jumps = uJumps;

    m_EBEllipticSU->getOperator(*a_cellBeta, *a_cellEta);

    ierr = labelCells(a_rhs, a_rhsCut);

    writeLevel(&a_rhs);
   
    #ifdef CH_USE_PETSC
    ierr = assembleSystem(a_phi, a_phiCut, a_rhs, a_rhsCut, *a_cellBeta, *a_cellEta);
    ierr = solveSystem(a_phi, a_phiCut);
    #endif

    #ifdef CH_USE_PETSC 
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    #endif

    return ierr;
}



#ifdef CH_USE_PETSC

//map from cell to matrix location using scan
#undef __FUNCT__
#define __FUNCT__ "labelCells"
PetscErrorCode
EBDriver::labelCells(const LevelData<FArrayBox>& a_rhs,
                    const LevelData<IVSFAB<Real>>& a_rhsCut)
{
    CH_TIME("EBDriver::labelCells");
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    

    DisjointBoxLayout dbl = a_rhs.disjointBoxLayout();
    cellIDs.define(dbl, 2, a_rhs.ghostVect());
    
    LevelData<IVSFAB<Real>>& moments = *m_moments;
    LevelData<NodeFArrayBox>& psiNodes = *m_psiNodes;

    PetscInt my0s[2] = {0,0}; //[starting id#, starting patch number]

    //get cell counts
    for(DataIterator dit = dbl.dataIterator();dit.ok();++dit )
    {
        Box region = dbl[dit];
        const IntVectSet& EBCellsFab = (moments[dit()]).getIVS();
        for(BoxIterator bit(region);bit.ok();++bit)
        {
            my0s[0] +=1; //1 DOF in each full cell
            if(EBCellsFab.contains(bit())) //1 more for cut cells
            {
                my0s[0] += 1;
            }
        }
        my0s[1] ++; //total patches
    }

    //do prefix scan to get start of this block
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Comm wcomm = Chombo_MPI::comm;
    PetscInt result[2];
    MPI_Datatype mtype;
    PetscDataTypeToMPIDataType(PETSC_INT,&mtype);
    MPI_Scan(my0s, result, 2, mtype, MPI_SUM, Chombo_MPI::comm);
    id0 = result[0] - my0s[0];
    patchid0 = result[1] - my0s[1]; 
    pout() << "rank " << world_rank << " id0 " << id0 << " patchid0 " << patchid0 << " \n";


    int id=0;
    for(DataIterator dit = dbl.dataIterator();dit.ok();++dit )
    {
        FArrayBox& cellIDFab = cellIDs[dit()];
        cellIDFab.setVal(-1); //petsc ignores -1s
        Box valid(cellIDFab.box());
        valid.grow(-1*cellIDs.ghostVect());
        IntVect iv;

        const IntVectSet& EBCellsFab = (moments[dit()]).getIVS();
        const IntVectSet& taggedCellsFab = (( *(m_EBEllipticSU->m_taggedFluxCells ) )[dit()]);

        IntVectSet validCutTagged = EBCellsFab; validCutTagged |= taggedCellsFab;
        int j;

        const FArrayBox& psiNodesFab = psiNodes[dit()].getFab();

        for(BoxIterator bit(valid);bit.ok();++bit)
        {
            iv = bit();
            if(EBCellsFab.contains(iv))
            {
                j = 2;
            }
            else
            {
                j=1;
            }

            for(int i =0;i<j;i++)
            {
                cellIDFab(iv,i) = id0 + id;
                id++;
            }
        } 
        //j+=10;
    }

    //writeLevel(&cellIDs);

    const ProblemDomain& pd = dbl.physDomain();
    
    cellIDs.exchange();

    //now set up the matrix and vectors
    ierr =MatCreate(wcomm,&A); CHKERRQ(ierr);
    ierr =MatSetSizes(A,my0s[0],my0s[0],PETSC_DECIDE,PETSC_DECIDE); CHKERRQ(ierr);
    
    ierr =MatSetFromOptions(A); CHKERRQ(ierr);//what options do I need to set other than aij

    //allocate memory - deal with this later, this is slightly wasteful but not bad - can do better on this if necessary
    PetscInt r_nnz;
    r_nnz =  2*( (2*m_velOrder+3) ) * ( (2*m_velOrder+3)); //max nnzs - use in cut cell rows
    if(HOCut)
    {
       r_nnz =  2*( (2*(m_velOrder+1) +3 ) ) * ( (2*(m_velOrder+1)+3)); //max nnzs - use in cut cell rows 
    }

    //r_nnz = 2*(m_velOrder+1)*(m_velOrder+1);
    ierr = MatMPIAIJSetPreallocation(A,r_nnz,NULL, r_nnz,NULL); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A,r_nnz,NULL); CHKERRQ(ierr);

    //preallocate enough memory for cut cell rows 
    PetscInt *d_nnz, *o_nnz;

    //set up vectors correctly - this seems to work
    pout() << "vec size: " << my0s[0] << "\n";
    ierr = VecCreate(wcomm, &x); CHKERRQ(ierr);
    ierr = VecSetSizes(x, my0s[0], PETSC_DECIDE); CHKERRQ(ierr);
    ierr = VecSetFromOptions(x); CHKERRQ(ierr);

    ierr = VecCreate(wcomm, &b); CHKERRQ(ierr);
    ierr = VecSetSizes(b, my0s[0], PETSC_DECIDE); CHKERRQ(ierr);
    ierr = VecSetFromOptions(b); CHKERRQ(ierr);

    //ierr = MatGetVecs(A,&x,&b); //CHKERRQ(ierr); //I think this only works if you have assembled mat
    
    //we can just form b right now actually

    PetscInt *ix, counter;
    PetscScalar *vals;

    for(DataIterator dit = dbl.dataIterator();dit.ok();++dit )
    {
        FArrayBox& cellIDFab = cellIDs[dit()];
        Box valid(cellIDFab.box());
        valid.grow(-1*cellIDs.ghostVect());
        IntVect iv;

        const FArrayBox& rhsFab = a_rhs[dit()];
        const IVSFAB<Real>& rhsCutFab = a_rhsCut[dit()];

        const IntVectSet& EBCellsFab = (moments[dit()]).getIVS();

        const IVSFAB<Real>& momentsFab = moments[dit()];
        const FArrayBox& psiNodesFab = psiNodes[dit()].getFab();

        Vector<Real> kappa(2,0);

        int j; int currID;
        // PetscInt ix [id];
        // PetscScalar y[id];
        counter = 0;

        //count number of points in this region so we can allocate memory
        for(BoxIterator bit(valid);bit.ok();++bit) 
        {
            iv  = bit();
            counter+=1;
            if(EBCellsFab.contains(iv))
            {
                counter +=1;
            }

        }
        
        //allocate memory
        PetscMalloc1(counter, &ix);
        PetscMalloc1(counter, &vals);

        counter =0; //reset

        for(BoxIterator bit(valid);bit.ok();++bit) 
        {
            iv  = bit();

            if(EBCellsFab.contains(iv))
            {
                j=2;
                if(psiNodesFab(iv) > 0) //lo side is grd
                {
                    kappa[0] = momentsFab(iv,0);
                    kappa[1] = 1.0 - kappa[0];
                }
                else //hi side is grd
                {
                    kappa[1] = momentsFab(iv,0);
                    kappa[0] = 1.0 - kappa[1];
                }
                for(int i =0;i<j;i++)
                {
                    currID = int( cellIDFab(iv,i) );
                    ix[counter] = currID;
                    if(1)
                    {
                        vals[counter] = kappa[i]*rhsCutFab(iv,i);
                    }
                    else
                    {
                        vals[counter] = rhsCutFab(iv,i);
                    }
                    
                    counter ++;
                }
            }
            else
            {
                j=1;
                for(int i =0;i<j;i++)
                {
                    currID = int( cellIDFab(iv,i) );
                    ix[counter] = currID;
                    vals[counter] = rhsFab(iv,i);
                    counter++;
                }
            }

        }


        //pout() << "counter: " << counter << "\n";
        ierr = VecSetValues(b, counter, ix, vals, INSERT_VALUES);

        PetscFree(vals);
        PetscFree(ix);

    }

    ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b); CHKERRQ(ierr);

    //debug check
    PetscInt low,high;

    ierr = VecGetOwnershipRange(b, &low,&high); CHKERRQ(ierr);
    pout() << "low and high: " << low << " " << high << "\n"; 

    PetscFunctionReturn(0);
}

//use weighted FD  - assembly step
#undef __FUNCT__
#define __FUNCT__ "assembleSystem"
PetscErrorCode 
EBDriver::assembleSystem(const LevelData<FArrayBox>& a_phi,
                        const LevelData<IVSFAB<Real>>& a_phiCut,
                        const LevelData<FArrayBox>& a_rhs,
                        const LevelData<IVSFAB<Real>>& a_rhsCut,
                        const LevelData<FArrayBox>& a_cellBeta,
                        const LevelData<FArrayBox>& a_cellEta)
{
    CH_TIME("EBDriver::assembleSystem");
    PetscFunctionBeginUser;

    #ifdef CH_MPI
    MPI_Comm wcomm = Chombo_MPI::comm;
    #else
    MPI_Comm wcomm = PETSC_COMM_SELF;
    #endif

    PetscErrorCode ierr;
    PetscViewer    viewer;
    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
    
    LevelData<IVSFAB<Real>>& moments = *m_moments;
    LevelData<NodeFArrayBox>& psiNodes = *m_psiNodes;

    PetscInt AStart, AEnd;

    ierr = MatGetOwnershipRange(A, &AStart, &AEnd); CHKERRQ(ierr);

    pout() << "Mat range: " << AStart << " " << AEnd << "\n";
    

    for(DataIterator dit = a_phi.dataIterator();dit.ok();++dit)
    {
        const FArrayBox& phiFab = a_phi[dit()]; //cell averaged phi
        const FArrayBox& rhsFab = a_rhs[dit()]; //cell averaged rhs
        const IVSFAB<Real>& phiCutFab = a_phiCut[dit()]; //4 comps, grd then floating
        const IVSFAB<Real>& rhsCutFab = a_rhsCut[dit()]; //4 comps, grd then floating

        // const FArrayBox& acoefFab = a_acoef[dit()]; //friction coefficient
     
        const IVSFAB<Real>& momentsFab = moments[dit()];
        const IntVectSet& EBCellsFab = momentsFab.getIVS();

        
        const FArrayBox& psiNodesFab = psiNodes[dit()].getFab();

        IVSFAB<EBStencil>& EBOpStencilFab = (*( (*m_EBEllipticSU).EBOpStencil) )[dit()];
        IVSFAB<EBStencil>& EBOpStencilCutFab = (*( (*m_EBEllipticSU).EBOpStencilCut) )[dit()];
        //IVSFAB<EBStencil>& BCStencilFab = (*( (*m_EBEllipticSU)BCStencil) )[dit()];

        const IntVectSet& EBOpCells = EBOpStencilFab.getIVS();

        Box valid = phiFab.box(); valid.grow(-1*a_phi.ghostVect());

        FArrayBox& cellIDFab = cellIDs[dit()];
        //Real beta, etax0, etax1, etay0, etay1;
        int counter;
        IntVect up(0,1); IntVect right(1,0); IntVect iv, sIv;

        Vector<Real> kappa(2,0);
        //Vector<Real> kappaScale(2,1);

        PetscScalar rJump [1];

        const FArrayBox& etaFab = (a_cellEta)[dit()];
        const FArrayBox& betaFab = (a_cellBeta)[dit()];

        int Q = m_velOrder;
        //int n = ((Q+1)*(Q+2) )/2;
        int r = (Q+1)/2;
        int row=0;
        int comp;
        int nCells = (2*r+1)*(2*r+1);

        int centerIdx = (nCells-1)/2;

        PetscScalar v[1000]; //holds data to be inserted in matrix
        PetscInt m = 2;
        PetscInt n;
        PetscInt idxm[2];
        PetscInt idxn[1000];


        LAPACKMatrix etaM(1,nCells);
        LAPACKMatrix betaM(1,nCells);

        IntVect position,sten;

        LAPACKMatrix etaTerm, betaTerm;

        Real etax0, etax1, etay0, etay1, beta;

        Real tol = 1e-12;
        PetscScalar rhsVal[1];
        PetscInt rhsIx[1];

        PetscScalar uVal[1];
        int phase;
        
        //try putting in exact solution for x to test residual
        for(BoxIterator bit(valid); bit.ok();++bit)
        {
            iv = bit();
            // if(iv[0]==26 && iv[1]==14)
            // {
            //     printf("Asdfasdf \n");
            // }

            if(EBCellsFab.contains(iv)) //cut cell
            {
                if(psiNodesFab(iv) > 0) //lo side is grd
                {
                    kappa[0] = momentsFab(iv,0);
                    kappa[1] = 1.0 - kappa[0];
                }
                else //hi side is grd
                {
                    kappa[1] = momentsFab(iv,0);
                    kappa[0] = 1.0 - kappa[1];
                }

                for(int p=0;p<2;p++) //for two phases
                {
                    for(int d=0;d<1;d++) //only one comp for poisson
                    {
                        counter = 0;
                        m = 1;
                        idxm[0] = cellIDFab(iv,p + d); //row in matrix
                        EBStencil& currEBStencil = EBOpStencilCutFab(iv,p + d);
                        
                        for(IVSIterator stenIt(currEBStencil.stencilIVSF.getIVS()); stenIt.ok();++stenIt)
                        {
                            if(EBCellsFab.contains(stenIt())) //its a cut cell with 2 dof
                            {
                                for(int i=0;i<2;i++)
                                {
                                    idxn[counter] = cellIDFab(stenIt(),i); //for phase we are grabbing from
                                    v[counter] = kappa[p] * currEBStencil(stenIt(),i)*(1.0/(kappa[p]*m_dx*m_dx));
                                    counter ++;
                                }
                            }
                            else //regular cell with 1 dof
                            {
                                for(int i=0;i<1;i++)
                                {
                                    idxn[counter] = cellIDFab(stenIt(),i); //for phase we are grabbing from
                                    v[counter] = kappa[p] * currEBStencil(stenIt(),i)*(1.0/(kappa[p]*m_dx*m_dx));
                                    counter ++;
                                }
                            }
                        }
                        rhsIx[0] = idxm[0];
                        rhsVal[0] =kappa[p] * ( (rhsCutFab(iv,p + d) ) -(1.0/(kappa[p]*m_dx*m_dx))*currEBStencil.jumpConstant );
                        uVal[0] = phiCutFab(iv,p);
                        //add contributions to RHS
                        ierr = VecSetValues(b, 1, rhsIx, rhsVal, INSERT_VALUES); CHKERRQ(ierr);
                        ierr = VecSetValues(x, 1, rhsIx, uVal, INSERT_VALUES); CHKERRQ(ierr);
                        
                        
                        ierr = MatSetValues(A, m, idxm, counter, idxn, v, INSERT_VALUES); CHKERRQ(ierr);
                        
                    }

                }
            
            
            }
            else if(EBOpStencilFab.getIVS().contains(iv)) //irregular cell
            {
                // if(iv[0]==13 && iv[1]==11)
                // {
                //     printf("qw4rt");
                // }
                
                for(int d=0;d<1;d++) //
                {
                    counter = 0;
                    m = 1;
                    idxm[0] = cellIDFab(iv,d);
                    EBStencil& currEBStencil = EBOpStencilFab(iv,d);
                    for(IVSIterator stenIt(currEBStencil.stencilIVSF.getIVS()); stenIt.ok();++stenIt)
                    {
                        if(EBCellsFab.contains(stenIt())) //its a cut cell with 2 dof
                        {
                            for(int i=0;i<2;i++)
                            {
                                idxn[counter] = cellIDFab(stenIt(),i); //for phase we are grabbing from
                                v[counter] = currEBStencil(stenIt(),i)*(1.0/(m_dx*m_dx));
                                counter ++;
                            }
                        }
                        else //regular cell with 1 dof
                        {
                            for(int i=0;i<1;i++)
                            {
                                idxn[counter] = cellIDFab(stenIt(),i); //for phase we are grabbing from
                                v[counter] = currEBStencil(stenIt(),i)*(1.0/(m_dx*m_dx));
                                counter ++;
                            }
                        }
                    }

                    ierr = MatSetValues(A, m, idxm, counter, idxn, v, INSERT_VALUES);  CHKERRQ(ierr);

                    rhsIx[0] = idxm[0];
                    rhsVal[0] = rhsFab(iv,d) - (1.0/(m_dx*m_dx))*currEBStencil.jumpConstant ;
                    uVal[0] = phiFab(iv,0);
                    //add contributions to RHS
                    ierr = VecSetValues(b, 1, rhsIx, rhsVal, INSERT_VALUES); CHKERRQ(ierr);
                    ierr = VecSetValues(x, 1, rhsIx, uVal, INSERT_VALUES); CHKERRQ(ierr);
                    
                                       
                }

            
                    
            }
            // else if(BCStencilFab.getIVS().contains(iv)) //BC cell
            // {

            // }
            else //regular cell
            {
                m = 1;

                idxm[0] = cellIDFab(iv,0);

                n = nCells;
                counter =0;

                phase=1;
                if(psiNodesFab(iv) >0 ){phase=0;}

                for(int d2 = 0;d2<1;d2++) //one components of u field
                {
                    row=0;
                    for(int i =-r;i<r+1;i++)
                    {
                        for(int j=-r;j<r+1;j++)
                        {
                            sIv[1] = iv[1] + i; 
                            sIv[0] = iv[0] + j;

                            idxn[counter] = cellIDFab(sIv, d2);
                            v[counter] =  0; //set to 0 
                            counter ++; 
                            
                            position[0] = j; position[1]= i;
                            sten = iv + position;
                            etaM(0,row) = etaFab(sten,phase);
                            betaM(0,row) = betaFab(sten,phase);
                            row+=1;
                        }
                    }
                }

                multiply(betaTerm, betaM, betaStencil);
                
                for(int i = 0;i<1;i++) //one comp of L
                {
                    multiply(etaTerm, etaM, etaStencil[0]); //u
                    for(int j=0;j<nCells;j++)
                    {
                        v[j] = etaTerm(0,j) * (1.0/(m_dx*m_dx));
                        
                    }

                    for(int j=0;j<nCells;j++)
                    {
                        v[j] += -1*betaTerm(0,j);
                    }
                }  
                rhsIx[0] = idxm[0];
                rhsVal[0] = rhsFab(iv,0); 
                uVal[0] = phiFab(iv,0);
                //add contributions to RHS
                ierr = VecSetValues(b, 1, rhsIx, rhsVal, INSERT_VALUES); CHKERRQ(ierr);
                ierr = VecSetValues(x, 1, rhsIx, uVal, INSERT_VALUES); CHKERRQ(ierr);
                
                //insert this stencil into ya matrix
                ierr = MatSetValues(A, m, idxm, n, idxn, v, INSERT_VALUES); CHKERRQ(ierr);
                
   
            }

        }

        
        //CHKERRQ(ierr);  
    }

    {
        CH_TIME("assembleMatrix::MAT ASSEMBLY");
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

        ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b); CHKERRQ(ierr);

        ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(x); CHKERRQ(ierr);
    }

    //test residual here with exact guess.
    //these lines are for debugging purposes, and only make sense for n <= 64

    // Vec res, LPhi;
    // ierr = MatCreateVecs(A, &LPhi, &res);  CHKERRQ(ierr);

    // ierr = MatMult(A, x, LPhi);  CHKERRQ(ierr);

    // ierr = VecCopy(LPhi, res); CHKERRQ(ierr);

    // ierr = VecAXPY(res, -1.0, b );  CHKERRQ(ierr);
    // ierr = VecScale(res, -1.0);   CHKERRQ(ierr);

    // PetscScalar resNorm;

    // ierr = VecNorm(res, NORM_INFINITY, &resNorm);    CHKERRQ(ierr);

    // LevelData<FArrayBox> residualLD(cellIDs.disjointBoxLayout(), 2, cellIDs.ghostVect() );
    // LevelData<FArrayBox> LPhiLD(cellIDs.disjointBoxLayout(), 2, cellIDs.ghostVect() );

    // writePETSC(res, cellIDs, residualLD);
    // writePETSC(LPhi, cellIDs, LPhiLD);  

    // writeLevel(&LPhiLD);

    PetscFunctionReturn(0);

    
}

#undef __FUNCT__
#define __FUNCT__ "solveSystem"
PetscErrorCode
EBDriver::solveSystem(LevelData<FArrayBox>& a_phi,
                        LevelData<IVSFAB<Real>>& a_phiCut)
{
    CH_TIME("EBDriver::solveSystem");
    PetscFunctionBeginUser;
    int N = ((m_velOrder+1)*(m_velOrder+2) )/2;
    KSP ksp;
    PetscErrorCode ierr;
    #ifdef CH_MPI
    MPI_Comm wcomm = Chombo_MPI::comm;
    #else
    MPI_Comm wcomm = PETSC_COMM_SELF;
    #endif

    //PetscFunctionBeginUser;

    PetscBool ism = PETSC_FALSE;

    ierr = KSPCreate(wcomm, &ksp); CHKERRQ(ierr);
    #if PETSC_VERSION_LT(3,5,0)
    ierr = KSPSetOperators(ksp, A, A, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
    #else
    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
    #endif
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    #if PETSC_VERSION_GE(3,7,0)
    PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-ksp_monitor",&ism,PETSC_NULL);
    #else
    PetscOptionsGetBool(PETSC_NULL,"-ksp_monitor",&ism,PETSC_NULL);
    #endif


    PC pc; 
    PetscInt gid,sz,bs,n,m,nGrids=1;
    // const int my0eq=CH_SPACEDIM*m_petscCompMat.m_gid0;
    #if PETSC_VERSION_LT(3,4,0) & PETSC_VERSION_RELEASE
    const PCType type;
    #else
    PCType type;
    #endif	

    ierr = KSPGetPC( ksp, &pc );     CHKERRQ(ierr);
    ierr = PCGetType( pc, &type );    CHKERRQ(ierr);
    ierr = MatGetBlockSize( A, &bs );               CHKERRQ( ierr );
    
    PetscViewer    viewer;
    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
    
    //KSPSetComputeEigenvalues(ksp,PETSC_TRUE);
    //ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE ); CHKERRQ(ierr);//resuse current x

    
    //CH_TIME("EBDriver::actual solve");
    ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr); //most of the time is preconditioner
    

    //this is to see what the hell is going on
    // Vec res, LPhi;
    // ierr = MatCreateVecs(A, &LPhi, &res);  CHKERRQ(ierr);

    // ierr = MatMult(A, x, LPhi);  CHKERRQ(ierr);

    // ierr = VecCopy(LPhi, res); CHKERRQ(ierr);

    // ierr = VecAXPY(res, -1.0, b );  CHKERRQ(ierr);
    // ierr = VecScale(res, -1.0);   CHKERRQ(ierr);

    // PetscReal realpart [50];
    // PetscReal complexpart [50];
    // PetscInt neig = 50;

    // KSPComputeEigenvalues(ksp,50,realpart,complexpart,&neig);

    // pout() << "eig = [ \n";
    // for(int i=0;i<neig;i++)
    // {
    //     pout() << realpart[i] << " + " << complexpart[i] << "i \n";
    // }
    // pout() << "]; \n";

    //need to find offset. the vector will be zero indexed 
    //but we need to subtract the beginning of this block
    PetscInt offset = id0;
    //pout() << "id0 from solve: " << id0 <<" \n";

    const PetscScalar *avec;
    // const PetscScalar *aLPhivec;
    // const PetscScalar *aRHSvec;
    // const PetscScalar *aResvec;
    //printf("Reading Vec \n");
    ierr = VecGetArrayRead(x,&avec); CHKERRQ(ierr);
    // ierr = VecGetArrayRead(LPhi,&aLPhivec); CHKERRQ(ierr);
    // ierr = VecGetArrayRead(b,&aRHSvec); CHKERRQ(ierr);
    // ierr = VecGetArrayRead(res,&aResvec); CHKERRQ(ierr);
    //printf("Read vec \n");
    LevelData<IVSFAB<Real>>& moments = *m_moments;
    LevelData<NodeFArrayBox>& psiNodes = *m_psiNodes;

    int currID; IntVect iv;
    Real volGrd,volFlt;

    //write solution info back to the LevelDatas
    for(DataIterator dit(a_phi.dataIterator());dit.ok();++dit)
    {
        FArrayBox& cellIDFab = cellIDs[dit()];
        FArrayBox& phiFab = a_phi[dit()];
        IVSFAB<Real>& phiCutFab = a_phiCut[dit()];
        const FArrayBox& psiNodesFab = psiNodes[dit()].getFab();

        IntVectSet EBCellsFab = phiCutFab.getIVS();

        const IVSFAB<Real>& momentsFab = moments[dit()];

        // FArrayBox lhsFab(phiFab.box(), 4); lhsFab.setVal(0);
        // FArrayBox rhsFab(phiFab.box(), 4); rhsFab.setVal(0);
        // FArrayBox resFab(phiFab.box(), 4); resFab.setVal(0);

        Box region = (a_phi.disjointBoxLayout() )[dit];

        for(BoxIterator bit(region);bit.ok();++bit)
        {
            iv = bit();
            if(EBCellsFab.contains(iv))
            {
                for(int i=0;i<2;i++)
                {
                    currID = cellIDFab(iv,i);
                    phiCutFab(iv,i) = avec[currID - offset];

                    // lhsFab(iv,i) = aLPhivec[currID- offset];
                    // rhsFab(iv,i) = aRHSvec[currID- offset];
                    // resFab(iv,i) = aResvec[currID- offset];

                }
                if(psiNodesFab(iv) > 0)
                {
                    volGrd = momentsFab(iv,0);
                    volFlt = 1.0 - volGrd;
                }
                else
                {
                    volFlt = momentsFab(iv,0);
                    volGrd = 1.0 - volFlt;
                }
                for(int i=0;i<1;i++)
                {
                    phiFab(iv,i) = volGrd*phiCutFab(iv,i) + volFlt*phiCutFab(iv,i +1);
                    //lhsFab(iv,i) = lhsFab(iv,i) + lhsFab(iv,i + 2);
                }

            }
            else
            {
                for(int i=0;i<1;i++)
                {
                    currID = cellIDFab(iv,i);
                    phiFab(iv,i) = avec[currID - offset];

                    // lhsFab(iv,i) = aLPhivec[currID- offset];
                    // rhsFab(iv,i) = aRHSvec[currID- offset];
                    // resFab(iv,i) = aResvec[currID- offset];
                }
            }
        }
        int asdg;
        asdg=0;
    }

    ierr = VecRestoreArrayRead(x,&avec);CHKERRQ(ierr);

    a_phi.exchange();
    a_phiCut.exchange();

    //printf("writ solution \n");

    //ierr = VecDestroy(&LPhi); 
    {
        CH_TIME("KSPSolve::KSPDestroy");
        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    }
    //PetscFunctionReturn(0);

    //printf("deestroyed ksp \n");

    PetscFunctionReturn(0);

    

}
#endif


#include "NamespaceFooter.H"