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
#include "SolutionSetUp.H"
#include "DebugDump.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "LevelDataOps.H"

//set up rhs, exact solution, and coefficient functions
void setUp(Real x0, Real y0, int order, Real dx,
           RealFunc& solFuncIn, RealFunc& solFuncOut, RealFunc& rhsFuncIn, RealFunc& rhsFuncOut,
           RealFunc& betaFuncIn, RealFunc& betaFuncOut, RealFunc& etaFuncIn, RealFunc& etaFuncOut,
           LevelData<FArrayBox>& cellU, LevelData<FArrayBox>& cellRhs, 
           LevelData<FArrayBox>& cellBeta, LevelData<FArrayBox>& cellEta,
           LevelData<IVSFAB<Real> >& cellUCut, LevelData<IVSFAB<Real> >& cellRhsCut,
           LevelData<IVSFAB<Real> >& cellBetaCut, LevelData<IVSFAB<Real> >& cellEtaCut,
           LevelData<IVSFAB<Real> >& uJumps,
           LevelData<NodeFArrayBox>& psiNodes, LevelData<IVSFAB<Real> >& moments)
{

        SolutionSetUp SSU(x0, y0, dx, order);
        //coefs

        SSU.setUpFullCellField(cellBeta, psiNodes, betaFuncIn, betaFuncOut, false, dx);
        SSU.setUpFullCellField(cellEta, psiNodes, etaFuncIn, etaFuncOut, false, dx);
        SSU.setUpCutCellField(cellBetaCut, moments, psiNodes, betaFuncIn, betaFuncOut, false, dx);
        SSU.setUpCutCellField(cellEtaCut, moments, psiNodes, etaFuncIn, etaFuncOut, false, dx);

        //transfer to levelDatas
        for(DataIterator dit = cellBeta.dataIterator(); dit.ok(); ++dit)
        {
            FArrayBox& cellBetaFab = cellBeta[dit()];
            FArrayBox& cellEtaFab = cellEta[dit()];
            IVSFAB<Real>& cellBetaCutFab = cellBetaCut[dit()];
            IVSFAB<Real>& cellEtaCutFab = cellEtaCut[dit()];

            for(IVSIterator ivIt(cellBetaCutFab.getIVS()); ivIt.ok(); ++ ivIt)
            {
                for(int i=0;i<2;i++)
                {
                    cellBetaFab(ivIt(), i) = cellBetaCutFab(ivIt(), i);
                    cellEtaFab(ivIt(), i) = cellEtaCutFab(ivIt(), i);
                }
            }
        }
        //RHS

        SSU.setUpFullCellField(cellRhs, psiNodes, rhsFuncIn, rhsFuncOut, true, dx);
        SSU.setUpCutCellField(cellRhsCut, moments, psiNodes, rhsFuncIn, rhsFuncOut, true, dx);
        
        //u
        SSU.setUpFullCellField(cellU, psiNodes, solFuncIn, solFuncOut, true, dx);
        SSU.setUpCutCellField(cellUCut, moments, psiNodes, solFuncIn, solFuncOut, true, dx);

        SSU.setUpJumps(uJumps, moments, psiNodes, solFuncIn, solFuncOut, etaFuncIn, etaFuncOut, dx);
}

int main(int argc, char* argv[])
{
    #ifdef CH_USE_PETSC
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc, &argv,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
    // dump loaded PETSc options into pout
    char *copts[1];
    #if PETSC_VERSION_GE(3,7,0)  
    PetscOptionsGetAll(PETSC_NULL,copts);
    #else
    PetscOptionsGetAll(copts);
    #endif
    
    pout() << copts[0] << "\n";
    #else
    #ifdef CH_MPI
    MPI_Init(&argc, &argv);
    #endif 
    #endif // end petsc conditional



    char* in_file = argv[1];
    ParmParse  pp(argc-2,argv+2,NULL,in_file);
    bool convergence, conserve,  factorEta, HOCut, noEB, noSolve;
    Real weight, betaInScale, betaOutScale, etaInScale, etaOutScale;
    int maxIter, recOrder, order, nCells, nCellsFine, geom;
    int solFuncInType, solFuncOutType;
    int betaFuncInType, betaFuncOutType, etaFuncInType, etaFuncOutType;
    int solN;

    Vector<Real> solCoefIn;
    Vector<Real> solCoefOut;
    Vector<Real> betaCoefOut;
    Vector<Real> betaCoefIn;
    Vector<Real> etaCoefOut;
    Vector<Real> etaCoefIn;

    std::string prefix;
    {
    ParmParse ppTest("test");
    ppTest.get("nCells",nCells);
    ppTest.get("nCellsFine",nCellsFine);
    ppTest.get("convergence",convergence);
    ppTest.get("noSolve", noSolve);
    ppTest.get("prefix",prefix);
    ppTest.get("geom", geom);

    ppTest.get("solFuncInType", solFuncInType );
    ppTest.get("solFuncOutType", solFuncOutType );

    ppTest.get("betaFuncInType", betaFuncInType );
    ppTest.get("betaFuncOutType", betaFuncOutType );

    ppTest.get("etaFuncInType", etaFuncInType );
    ppTest.get("etaFuncOutType", etaFuncOutType );

    ppTest.get("solN", solN);

    ppTest.getarr("solCoefIn", solCoefIn, 0, solN);
    ppTest.getarr("solCoefOut", solCoefOut, 0, solN);

    ppTest.getarr("betaCoefIn", betaCoefIn, 0, solN);
    ppTest.getarr("betaCoefOut", betaCoefOut, 0, solN);

    ppTest.getarr("etaCoefIn", etaCoefIn, 0, solN);
    ppTest.getarr("etaCoefOut", etaCoefOut, 0, solN);

    ppTest.get("betaInScale", betaInScale);
    ppTest.get("betaOutScale", betaOutScale);

    ppTest.get("etaInScale", etaInScale);
    ppTest.get("etaOutScale", etaOutScale);

    for(int i=0;i<solN;i++)
    {
        betaCoefIn[i]*=betaInScale; betaCoefOut[i]*=betaOutScale;
        etaCoefIn[i]*=etaInScale; etaCoefOut[i]*=etaOutScale;
    }
    ParmParse ppRec("rec"); 
    ppRec.get("order", order);
    ppRec.get("weight", weight);
    ppRec.get("factorEta", factorEta);
    ppRec.get("recOrder", recOrder);
    ppRec.get("HOCut", HOCut);
    ppRec.get("conserve",conserve);
    ppRec.get("noEB", noEB);
}
    
    if(!convergence)
    {  
        Real dx = 2.0 / (nCells);

        Real x0=-1; Real y0=-1; //lower left corner of domain [-1:1]^2

        int n = ((order+1)*(order+2) )/2;
        int r = order+1; 
        if(HOCut){r+=1;}

        Box b(IntVect::Zero,( (2.0/dx)-1)*IntVect::Unit);
        Vector<Box> boxes(1,b); Vector<int> proc_id(1,0);
        //divide up domain
        domainSplit(b, boxes, 64, 4); 
        LoadBalance( proc_id, boxes );

        bool per[2] ={true,true};

        ProblemDomain pd(b,per);
        DisjointBoxLayout dbl(boxes,proc_id,pd);

        Vector<Box> boxesS(1,b);
        Vector<int> proc_idS(1,0);
        DisjointBoxLayout dblSerial(boxesS, proc_idS, pd);

        //set options for discretization
        EBDriver EBRD;
        EBRD.m_velOrder = order;
        EBRD.m_recOrder = recOrder;
        EBRD.m_weight = weight;
        EBRD.factorEta = factorEta;
        EBRD.conserve = conserve;
        EBRD.HOCut = HOCut;
        EBRD.noEB = noEB;

        EBRD.define();

        //create data holders
        LevelData<FArrayBox> cellH(dbl,1,r*IntVect::Unit);

        Real radx = M_PI / 5.0;
        Real rady = M_PI / 5.0;

        EllipsoidIF psiE( 1.0/(radx*radx), 1.0/(rady*rady), 0, 0);

        annulusIF psiA;

        rhodoneaIF psiR;

        constIF psiC;
        
        //reconstruct EB
        if(geom==0) //ellipse
        {
            EBRD.reconstructEB(cellH, psiE, x0, y0, dx);
        }
        else if(geom==1) //annulus
        {
            EBRD.reconstructEB(cellH, psiA, x0, y0, dx);
        }
        else if(geom==2) //rhodonea
        {
            EBRD.reconstructEB(cellH, psiR, x0, y0, dx);
        }
        else //no geometry at all
        {
            EBRD.reconstructEB(cellH, psiC, x0, y0, dx);
        }
        
        //set up mpre data holders
        LevelData<IVSFAB<Real> >& moments = *(EBRD.m_moments);

        LevelData<FArrayBox> cellU(dbl,1,r*IntVect::Unit);
        LevelData<FArrayBox> cellRhs(dbl,1,r*IntVect::Unit);
        RefCountedPtr< LevelData<FArrayBox>> cellEta(new LevelData<FArrayBox>(dbl,2,r*IntVect::Unit));//friction coefficient
        RefCountedPtr< LevelData<FArrayBox>> cellBeta(new LevelData<FArrayBox>(dbl,2,r*IntVect::Unit));//friction coefficient

        LayoutData<IntVectSet> EBCells(dbl);
        IVSFABFactory<Real> ivf(EBCells);
        LevelData<IVSFAB<Real> >cellUCut(dbl,2,r*IntVect::Unit,ivf);
        LevelData<IVSFAB<Real> >cellRhsCut(dbl,2,r*IntVect::Unit,ivf);
        LevelData<IVSFAB<Real> >cellBetaCut(dbl,2,r*IntVect::Unit,ivf);
        LevelData<IVSFAB<Real> >cellEtaCut(dbl,2,r*IntVect::Unit,ivf);

        RefCountedPtr<LevelData<IVSFAB<Real> > > uJumps (new LevelData<IVSFAB<Real> >(dbl,2,r*IntVect::Unit, ivf) );
        LevelData<NodeFArrayBox>& psiNodes = *(EBRD.m_psiNodes);

        LevelDataOps<FArrayBox> LDO;

        int pQ, N;
        pQ = sqrt(solN) -1;
        N = (sqrt(solN) - 1)/2;
        Real w = M_PI;
        //pointers for function set up. the idea here is we will keep coming with new options for functions
        RealFunc *solFuncIn, *solFuncOut, *rhsFuncIn, *rhsFuncOut, *betaFuncIn,  *betaFuncOut, *etaFuncIn, *etaFuncOut;
        
        polyFunc betaPolyIn(pQ, betaCoefIn);
        periodicFunc betaPerIn(betaCoefIn, N,w);

        //beta

        if(betaFuncInType == 0){betaFuncIn = &betaPolyIn;}
        else if(betaFuncInType == 1){betaFuncIn = &betaPerIn;}
        else{CH_assert(false);}
        
        polyFunc betaPolyOut(pQ, betaCoefOut);
        periodicFunc betaPerOut(betaCoefOut, N,w);

        if(betaFuncOutType == 0){betaFuncOut = &betaPolyOut;}
        else if(betaFuncOutType == 1){betaFuncOut = &betaPerOut;}
        else{CH_assert(false);}

        //eta

        polyFunc etaPolyIn(pQ, etaCoefIn);
        periodicFunc etaPerIn(etaCoefIn, N,w);

        if(etaFuncInType == 0){etaFuncIn = &etaPolyIn;}
        else if(etaFuncInType == 1){etaFuncIn = &etaPerIn;}
        else{CH_assert(false);}
        
        polyFunc etaPolyOut(pQ, etaCoefOut);
        periodicFunc etaPerOut(etaCoefOut, N,w);

        if(etaFuncOutType == 0){etaFuncOut = &etaPolyOut;}
        else if(etaFuncOutType == 1){etaFuncOut = &etaPerOut;}
        else{CH_assert(false);}

        //solution

        polyFunc solPolyIn(pQ, solCoefIn);
        periodicFunc solPerIn(solCoefIn, N,w);

        if(solFuncInType == 0){solFuncIn = &solPolyIn;}
        else if(solFuncInType == 1){solFuncIn = &solPerIn;}
        else{CH_assert(false);}
        
        polyFunc solPolyOut(pQ, solCoefOut);
        periodicFunc solPerOut(solCoefOut, N,w);

        if(solFuncOutType == 0){solFuncOut = &solPolyOut;}
        else if(solFuncOutType == 1){solFuncOut = &solPerOut;}
        else{CH_assert(false);}

        RhsFunc rFuncIn(betaFuncIn, etaFuncIn, solFuncIn); 
        rhsFuncIn = &rFuncIn;
        RhsFunc rFuncOut(betaFuncOut, etaFuncOut, solFuncOut); 
        rhsFuncOut = &rFuncOut;

        //set up solution
        SolutionSetUp SSU(x0, y0, dx, order);
        //coefs

        SSU.setUpFullCellField(*cellBeta, psiNodes, *betaFuncIn, *betaFuncOut, false, dx);
        SSU.setUpFullCellField(*cellEta, psiNodes, *etaFuncIn, *etaFuncOut, false, dx);
        SSU.setUpCutCellField(cellBetaCut, moments, psiNodes, *betaFuncIn, *betaFuncOut, false, dx);
        SSU.setUpCutCellField(cellEtaCut, moments, psiNodes, *etaFuncIn, *etaFuncOut, false, dx);

        //transfer to levelDatas
        for(DataIterator dit = cellBeta->dataIterator(); dit.ok(); ++dit)
        {
            FArrayBox& cellBetaFab = (*cellBeta)[dit()];
            FArrayBox& cellEtaFab = (*cellEta)[dit()];
            IVSFAB<Real>& cellBetaCutFab = cellBetaCut[dit()];
            IVSFAB<Real>& cellEtaCutFab = cellEtaCut[dit()];

            for(IVSIterator ivIt(cellBetaCutFab.getIVS()); ivIt.ok(); ++ ivIt)
            {
                for(int i=0;i<2;i++)
                {
                    cellBetaFab(ivIt(), i) = cellBetaCutFab(ivIt(), i);
                    cellEtaFab(ivIt(), i) = cellEtaCutFab(ivIt(), i);
                }
            }
        }
        //RHS

        SSU.setUpFullCellField(cellRhs, psiNodes, *rhsFuncIn, *rhsFuncOut, true, dx);
        SSU.setUpCutCellField(cellRhsCut, moments, psiNodes, *rhsFuncIn, *rhsFuncOut, true, dx);
        
        //u
        SSU.setUpFullCellField(cellU, psiNodes, *solFuncIn, *solFuncOut, true, dx);
        SSU.setUpCutCellField(cellUCut, moments, psiNodes, *solFuncIn, *solFuncOut, true, dx);

        SSU.setUpJumps(*uJumps, moments, psiNodes, *solFuncIn, *solFuncOut, *etaFuncIn, *etaFuncOut, dx);

        // setUp(x0, y0, order, dx, *solFuncIn, *solFuncOut, *rhsFuncIn, *rhsFuncOut,
        // *betaFuncIn, *betaFuncOut, *etaFuncIn, *etaFuncOut,
        // cellU, cellRhs, *cellBeta, *cellEta,
        // cellUCut, cellRhsCut, cellBetaCut, cellEtaCut, *uJumps,
        // psiNodes, moments );

        //JCs and BCs need to get in here somehow ??
        if(!noSolve)
        {
            EBRD.solver(cellBeta, cellEta, uJumps, cellU, cellUCut, cellRhs, cellRhsCut);
        }

        //do convergence test and actually look at cut cell error
        //[phi, phiCut, rhs, rhsCut, beta, betaOut, eta, etaOut, psiNodes, moments]
        //just use the convergence test in richardson to calculate this pretty much
        //could consider using the EBCellStuff here to make nicer plots

        //[uGrd, uFlt, rhsGrd, rhsFlt, betaGrd, betaFlt, etaGrd, etaFlt, psiNodes, grdMoment, ]

        LevelData<FArrayBox> output(dbl,10,r*IntVect::Unit);
        LDO.setVal(output, -1);
        LevelData<FArrayBox> outputSerial(dblSerial,10,r*IntVect::Unit);

        for(DataIterator dit = output.dataIterator();dit.ok();++dit)
        {
            IVSFAB<Real>& momentsFab = moments[dit()];
            IntVectSet EBCellsFab = momentsFab.getIVS();
            FArrayBox&  psiNodesFab = ( (*(EBRD.m_psiNodes) ) [dit()]).getFab();
            FArrayBox& uFab = cellU[dit()];
            FArrayBox& rhsFab = cellRhs[dit()];
            FArrayBox& betaFab = (*cellBeta)[dit()];
            FArrayBox& etaFab = (*cellEta)[dit()];
            IVSFAB<Real>& uCutFab = cellUCut[dit()];
            IVSFAB<Real>& rhsCutFab = cellRhsCut[dit()];

            IntVect iv;

            int c;
            FArrayBox& outputFab = output[dit()];
            outputFab.setVal(-1, 9); //set all moments to neg 1
            for(BoxIterator bit( outputFab.box() );bit.ok();++bit)
            {
                iv = bit();

                if(EBCellsFab.contains(iv) )
                {
                    if(psiNodesFab(iv) > 0){outputFab(iv,9) = momentsFab(iv,0);}
                    else{outputFab(iv,9) = (1.0 - momentsFab(iv,0)); }

                    outputFab(iv,0) = uCutFab(iv,0);
                    outputFab(iv,1) = uCutFab(iv,1);
                    outputFab(iv,2) = rhsCutFab(iv,0);
                    outputFab(iv,3) = rhsCutFab(iv,1);
                    outputFab(iv,4) = betaFab(iv,0);
                    outputFab(iv,5) = betaFab(iv,1);
                    outputFab(iv,6) = etaFab(iv,0);
                    outputFab(iv,7) = etaFab(iv,1);
                    outputFab(iv,8) = psiNodesFab(iv);

                }
                else
                {
                    if(psiNodesFab(iv,0) > 0){c=0;}
                    else{c=1;}
                    for(c=0;c<2;c++)
                    {
                        outputFab(iv,c) = uFab(iv,0);
                        outputFab(iv,2+c) = rhsFab(iv,0);
                        outputFab(iv,4+c) = betaFab(iv,c);
                        outputFab(iv,6+c) = etaFab(iv,c);
                        outputFab(iv,8) = psiNodesFab(iv);
                    }
                    if(psiNodesFab(iv,0) > 0)
                    {
                        outputFab(iv,4) = betaFab(iv,0);
                        outputFab(iv,6) = etaFab(iv,0);
                    }
                    else
                    {
                        outputFab(iv,4) = betaFab(iv,1);
                        outputFab(iv,6) = etaFab(iv,1);
                    }
                    

                }
            }

        }

        output.exchange();
        output.copyTo(outputSerial);
        
        std::string phiS = "output";
        phiS += std::to_string(dx);
        phiS += ".hdf5";
        int nA = phiS.length();
        char phiA[nA+1];
        strcpy(phiA,phiS.c_str());
        writeLevelname(&outputSerial, phiA);    
    }
    else //run convergence test
    {
        Real dx = 2.0/ nCells;
        Real dxf = 2.0 / nCellsFine;
        std::string file_namec = prefix;
        file_namec += std::to_string(dx); 
        file_namec += ".hdf5";
        Vector<DisjointBoxLayout> vectGridsc;
        Vector<LevelData<FArrayBox>* > vectDatac;
        Vector<string> vectNamesc;
        Box domainc;
        Vector<int> refRatioc;
        int numLevelsc;
        Real dxc, dtc,timec;
        ReadAMRHierarchyHDF5(file_namec,vectGridsc,vectDatac,vectNamesc,domainc,dxc,dtc,timec,refRatioc,numLevelsc);

        LevelData<FArrayBox>& outputc = *(vectDatac[0]);

        std::string file_namef = prefix;
        file_namef += std::to_string(dxf); 
        file_namef += ".hdf5";
        Vector<DisjointBoxLayout> vectGridsf;
        Vector<LevelData<FArrayBox>* > vectDataf;
        Vector<string> vectNamesf;
        Box domainf;
        Vector<int> refRatiof;
        int numLevelsf;
        Real dxF, dtf,timef;
        ReadAMRHierarchyHDF5(file_namef,vectGridsf,vectDataf,vectNamesf,domainc,dxF,dtf,timef,refRatiof,numLevelsf);

        LevelData<FArrayBox>& outputf = *(vectDataf[0]);

        LAPACKMatrix err;

        convergenceTest(outputc, outputf, dx, dxf, order, err);

        print_matrix(err, 16);
    }


    CH_TIMER_REPORT(); 
    #ifdef CH_USE_PETSC
    ierr = PetscFinalize(); CHKERRQ(ierr);
    #else
    #ifdef CH_MPI
    MPI_Finalize();
    #endif // mpi conditional
    #endif
    int s =0;
    return s; 
}