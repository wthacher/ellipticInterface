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
#include "DebugDump.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"


static char help[] = "Appends to an ASCII file.\n\n";

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
    Real rad, coefIn, coefOut, omegaF, omegaG;
    bool piecewiseL,center,exactGeom,useArea,convergence;
    Real weight;
    Real tol;
    int maxIter;
    int Q=2;//order of polynomial. 
    bool conserve;
    bool useGradJump;
    Real lambdaStart, penalty, scale, betaFloat;
    int aQ, recOrder;
    int solverType, relaxIter, maxOuterIt, nCells, nCellsFine;
    Real eps, relaxScale, hang, picardTol, muInitFloat;
    Real rho,rhoIce, g, rateFactor, muInit, beta;
    bool verbose; int maxIterp;
    bool useTagged;int order;bool momC,truncationTest,flatIce,useEBFlux,useMF, useConstrained, usePolyFits;
    bool useJumpConstraint,useCentroid,useUJump, constFriction, factorMu; int boundaryCond;
    bool KS, GRS;
    bool HOCut,useL, noEB;
    std::string prefix;
    {
    ParmParse ppTest("test");
    ppTest.get("nCells",nCells);
    ppTest.get("nCellsFine",nCellsFine);
    ppTest.get("convergence",convergence);
    ppTest.get("prefix",prefix);
    ppTest.get("bc",boundaryCond);
  

    ParmParse ppRec("rec"); 
    ppRec.get("order", order);
    ppRec.get("weight", weight);
    ppRec.get("useGradJump",useGradJump);
    ppRec.get("factorMu", factorMu);
    ppRec.get("recOrder", recOrder);
    ppRec.get("HOCut", HOCut);
    ppRec.get("conserve",conserve);
    ppRec.get("useL",useL);
    ppRec.get("noEB", noEB);

    ParmParse ppPhys("phys");
    ppPhys.get("muInit",muInit);
    ppPhys.get("beta",beta);
    ppPhys.get("betaFloat",betaFloat);
    ppPhys.get("muInitFloat", muInitFloat);

    ParmParse ppSolver("solver");
    
    ppSolver.get("maxOuterIt",maxOuterIt);
    ppSolver.get("picardTol",picardTol);
}
    
    Real radx = M_PI / 10.0;
    Real rady = M_PI / 20.0;

    //EllipsoidIF psi( 1.0/(radx*radx), 1.0/(rady*rady), .5, .5);

    annulusIF psi;

    //rhodoneaIF psi;

    Real dxc = 2.0 / (nCells);
    Real dxf = 2.0 / (nCellsFine);
    //printf("dx: %1.10e \n", dx);

    int n = ((order+1)*(order+2) )/2;
    int r = order+1; //to make me lfye easier
    if(HOCut){r+=1;}

    Box bc(IntVect::Zero,( (2.0/dxc)-1)*IntVect::Unit);
    Vector<Box> boxesc(1,bc); Vector<int> proc_idc(1,0);

    bool per[2] ={true,true};

    ProblemDomain pdc(bc,per);
    DisjointBoxLayout dblc(boxesc,proc_idc,pdc);

    //set options
    EBDriver EBRDc;
    EBRDc.m_velOrder = order;
    EBRDc.m_recOrder = recOrder;

    EBRDc.define();

    //create data holders
    LevelData<FArrayBox> cellHc(dblc,1,r*IntVect::Unit);
    
    //reconstruct EB, set up data holders for solve
    EBRDc.reconstructEB(cellHc, psi, -1, -1, dxc);

    LevelData<IVSFAB<Real> >& momentsc = *(EBRDc.m_moments);




    ////fine level
    Box bf(IntVect::Zero,((2.0/dxf)-1)*IntVect::Unit);
    Vector<Box> boxesf(1,bf); Vector<int> proc_idf(1,0);


    ProblemDomain pdf(bf,per);
    DisjointBoxLayout dblf(boxesf,proc_idf,pdf);

    //set options
    EBDriver EBRDf;
    EBRDf.m_velOrder = order;
    EBRDf.m_recOrder = recOrder;

    EBRDf.define();

    //create data holders
    LevelData<FArrayBox> cellHf(dblf,1,r*IntVect::Unit);

    //radius .5 centered at .5,.5
    
    //reconstruct EB, set up data holders for solve
    EBRDf.reconstructEB(cellHf, psi, -1, -1, dxf);

    LevelData<IVSFAB<Real> >& momentsf = *(EBRDf.m_moments);

    //now examine errors


    IntVect right(1,0); IntVect up(0,1);
    Real fineSol, coarseSol,fineSolArea, coarseSolArea;



    //now compare the coarse solution with the fine
    Vector<Real> errors(6,0); Real nReg=0; int numGradFacesGrd, numGradFacesFlt; 
    
    Vector<Real> momErrors(6,0); //vol*3, area*3
    
    Vector<Real> cutCellErrors(6,0);
    Vector<Real> taggedCellErrors(6,0);

    Real totalGrdArea=0;
    Real totalEBArea = 0;
    Real totalEBAreaCoarse = 0;
    Real nTagged = 0; Real nCut = 0;
    Real nRhsReg = 0; Real nRhsCut = 0; Real nRhsTagged = 0;
    IntVect iv;
    int velOrder = order;
    int nVel = n;
    for(DataIterator dit = dblc.dataIterator();dit.ok();++dit)
    {
        FArrayBox& psiNodesc = (*(EBRDc.m_psiNodes))[dit()].getFab();
        IVSFAB<Real>& momFabc = momentsc[dit()];
        FluxBox& faceIntValscFab = (*(EBRDc.m_faceIntVals))[dit()];

        Box validCoarse = dblc[dit()];
        
        for(DataIterator dit = dblf.dataIterator();dit.ok();++dit)
        {
            FArrayBox& psiNodesf = (*(EBRDf.m_psiNodes))[dit()].getFab();
            IVSFAB<Real>& momFabf = momentsf[dit()];
            FluxBox& faceIntValsfFab = (*(EBRDf.m_faceIntVals))[dit()];

            
            //validCoarse.grow(-(r+2*(1.0/dxc))*IntVect::Unit); //so we are covering same domain
            FArrayBox dif(validCoarse,4);
            dif.setVal(0);
            FArrayBox momDif(validCoarse,nVel);
            momDif.setVal(0);
            FArrayBox rhsDif(validCoarse,2);
            rhsDif.setVal(0);


            Vector<Real> momC(nVel,0); Vector<Real> momCa(nVel,0); Vector<Real> momCnx(nVel,0);  Vector<Real> momCny(nVel,0);  
            Vector<Real> momF(nVel,0); Vector<Real> momFa(nVel,0); Vector<Real> momFnx(nVel,0);  Vector<Real> momFny(nVel,0);  
            Vector<Real> momFTot(nVel,0); Vector<Real> momFTota(nVel,0); Vector<Real> momFTotnx(nVel,0); Vector<Real> momFTotny(nVel,0);
            Vector<Real> momErrorsAll(nVel,0); Vector<Real> momErrorsAlla(nVel,0); Vector<Real> momErrorsAllnx(nVel,0); Vector<Real> momErrorsAllny(nVel,0);
            Real sx,sy,x,x2,y,y2,momFull;
            int comp;
            int nCut = 0;
            bool llGrdC = false;
            int nSign = 0;
            if(true) //look at convergence of volume/EB length moments for positive or grd parts of cell
            {
                IntVectSet validEB  = momFabc.getIVS();
                validEB &= validCoarse;
                Box fineBox;
                for(IVSIterator ivIt(validEB);ivIt.ok();++ivIt)
                {
                    nCut += 1;
                    iv = ivIt();
                    momFTot.assign(0); momFTota.assign(0); momFTotnx.assign(0); momFTotny.assign(0);

                    llGrdC = psiNodesc(iv ) > 0;
                    
                    if(psiNodesc(iv) > 0)
                    {
                        for(int px =0; px<=velOrder; px++)
                        {
                            for(int py = 0; py<=velOrder-px;py++)
                            {
                                comp=  px*velOrder + (3*px - px*px)/2 + py;
                                momC[comp] = momFabc(iv,comp)*pow(dxc,2+px+py);
                                momCa[comp] = momFabc(iv,comp + 2*nVel)*pow(dxc,1+px+py);
                                momCnx[comp] = ( momFabc(iv, comp + 3*nVel)*pow(dxc,1+px+py) );
                                momCny[comp] = ( momFabc(iv, comp + 4*nVel)*pow(dxc,1+px+py) );
                            }
                        }
                            
                        
                        coarseSol = momFabc(iv,0)*dxc*dxc; //grd area on coarse level  
                    }
                    else
                    {
                        for(int px =0; px<=velOrder; px++)
                        {
                            for(int py = 0; py<=velOrder-px;py++)
                            {
                                comp=  px*velOrder + (3*px - px*px)/2 + py;
                                momC[comp] = momFabc(iv,comp+nVel)*pow(dxc,2+px+py);
                                momCa[comp] = momFabc(iv,comp + 2*nVel)*pow(dxc,1+px+py);
                                momCnx[comp] = ( momFabc(iv, comp + 3*nVel)*pow(dxc,1+px+py) );
                                momCny[comp] = ( momFabc(iv, comp + 4*nVel)*pow(dxc,1+px+py) );
                            }
                        }
                        coarseSol = momFabc(iv,nVel)*dxc*dxc;
                    }
                    coarseSolArea = momFabc(iv,2*nVel)*dxc;

                    fineBox.define(2*iv,2*iv +IntVect::Unit);
                    fineSol = 0;
                    fineSolArea = 0;
                    for(BoxIterator bit(fineBox);bit.ok();++bit)
                    {
                        sx = dxf* ( -2*iv[0] +  bit()[0] - .5 );
                        sy = dxf* ( -2*iv[1] +  bit()[1] - .5 ) ;
                        if(momFabf.getIVS().contains(bit()))
                        {
                            fineSolArea += momFabf(bit(),2*nVel)*dxf;
                            nSign =  (psiNodesf(bit()) > 0) == llGrdC ? 1 : -1;

                            if(psiNodesf(bit()) > 0)
                            {
                                fineSol += momFabf(bit(),0)*dxf*dxf; //grd area on coarse level  
                                for(int i = 0;i<nVel;i++)
                                {
                                    momF[i] = momFabf(bit(),i);
                                    momFa[i] = momFabf(bit(),i + 2*nVel);
                                    momFnx[i] = nSign*momFabf(bit(), i + 3*nVel) ;
                                    momFny[i] = nSign*momFabf(bit(), i + 4*nVel) ;
                                }
                            }
                            else
                            {
                                fineSol += momFabf(bit(),nVel)*dxf*dxf;
                                for(int i = 0;i<nVel;i++)
                                {
                                    momF[i] = momFabf(bit(),i+nVel);
                                    momFa[i] = momFabf(bit(),i + 2*nVel);
                                    momFnx[i] = nSign*momFabf(bit(), i + 3*nVel) ;
                                    momFny[i] = nSign*momFabf(bit(), i + 4*nVel) ;
                                }
                            }
                            
                            for(int px =0; px<=velOrder; px++)
                            {
                                for(int py = 0; py<=velOrder-px;py++)
                                {
                                    comp=  px*velOrder + (3*px - px*px)/2 + py;
                                    momF[comp] = momF[comp]*(pow(dxf,2+px+py));
                                    momFa[comp] = momFa[comp]*(pow(dxf,1+px+py));
                                    momFnx[comp] = momFnx[comp]*(pow(dxf,1+px+py));
                                    momFny[comp] = momFny[comp]*(pow(dxf,1+px+py));

                                }
                            }
                            shiftMoments(momF,velOrder,sx,sy); //convert units and shift
                            shiftMoments(momFa,velOrder,sx,sy); //convert units and shift
                            shiftMoments(momFnx,velOrder,sx,sy); //convert units and shift
                            shiftMoments(momFny,velOrder,sx,sy); //convert units and shift

                            for(int i = 0;i<nVel;i++)
                            {
                                momFTot[i] += momF[i];
                                momFTota[i] += momFa[i];
                                momFTotnx[i] += momFnx[i];
                                momFTotny[i] += momFny[i];
                            }
                        }
                        else
                        {
                            if(psiNodesf(bit()) > 0) //add full grd cells
                            {
                                fineSol += dxf*dxf;

                                x = sx - .5*dxf; x2 = sx + .5*dxf;
                                y = sy-.5*dxf; y2 = sy + .5*dxf;
                                for(int px =0; px<=velOrder; px++)
                                {
                                    for(int py = 0; py<=velOrder-px;py++)
                                    {
                                        comp=  px*velOrder + (3*px - px*px)/2 + py;
                                        momFull = (pow(x2,px+1) - pow(x,px+1) ) / (px+1);
                                        momFull*= (pow(y2,py+1) - pow(y,py+1) ) / (py+1);
                                        momF[comp] = momFull;
                                    }
                                }
                                for(int i = 0;i<nVel;i++)
                                {
                                    momFTot[i] += momF[i];
                                }

                            }
                        }
                        
                        
                    }
                    //momDif(iv,0) = abs(fineSol - coarseSol);
                    //momDif(iv,1) = abs(fineSolArea - coarseSolArea);
                    momErrors[0] += abs(fineSol-coarseSol)*fineSol;
                    momErrors[1] += pow(fineSol-coarseSol,2)*fineSol;
                    momErrors[2] = std::max(momErrors[2],abs(fineSol-coarseSol));

                    momErrors[3] += abs(fineSolArea-coarseSolArea)*fineSolArea;
                    momErrors[4] += pow(fineSolArea-coarseSolArea,2)*fineSolArea;
                    momErrors[5] = std::max(momErrors[5],abs(fineSolArea-coarseSolArea));

                    totalGrdArea += fineSol;
                    totalEBArea += fineSolArea;
                    totalEBAreaCoarse += coarseSolArea;

                    for(int i = 0;i<nVel;i++)
                    {
                        momErrorsAll[i] += abs(momFTot[i] - momC[i]);
                        momErrorsAlla[i] += abs(momFTota[i] - momCa[i]);
                        momErrorsAllnx[i] += abs(momFTotnx[i] - momCnx[i]);
                        momErrorsAllny[i] += abs(momFTotny[i] - momCny[i]);
                        momDif(iv,i) = abs(momFTot[i] - momC[i]);
                    }

                    
                }
            }

            int asd;
            asd+=1;
            printf("volMoms%d = [", int(1.0/dxc) );
            for(int i =0;i<nVel;i++)
            {
                printf("%1.6e, " , momErrorsAll[i]/nCut);
            }
            printf("] \n");

            printf("volMomsa%d = [", int(1.0/dxc) );
            for(int i =0;i<nVel;i++)
            {
                printf("%1.6e, " , momErrorsAlla[i]/nCut);
            }
            printf("] \n");

            printf("volMomsnx%d = [", int(1.0/dxc) );
            for(int i =0;i<nVel;i++)
            {
                printf("%1.6e, " , momErrorsAllnx[i]/nCut);
            }
            printf("] \n");

            printf("volMomsny%d = [", int(1.0/dxc) );
            for(int i =0;i<nVel;i++)
            {
                printf("%1.6e, " , momErrorsAllny[i]/nCut);
            }
            printf("] \n");

            //get errors in face interections
            Real faceIntError = 0;
            int nFace =  0;
            Real faceIntCoarse, faceIntFine;
            validCoarse.grow(-1*IntVect::Unit);
            Real nodeError = 0;
            int nNodes = 0;
            for(BoxIterator bit(validCoarse);bit.ok();++bit)
            {
                iv = bit();
                nodeError += abs( psiNodesc(iv) - psiNodesf(2*iv) );
                nNodes += 1;
                if(psiNodesc(iv)*psiNodesc(iv+right) < 0)
                {
                    nFace +=1;
                    faceIntCoarse = dxc * faceIntValscFab[1](iv);//  ( psiNodesc(iv) / (psiNodesc(iv) - psiNodesc(iv+right) ) );
                    if(psiNodesf(2*iv)*psiNodesf(2*iv+right) < 0)
                    {
                        faceIntFine =  dxf * faceIntValsfFab[1](2*iv); //( psiNodesf(2*iv) / (psiNodesf(2*iv) - psiNodesf(2*iv+right) ) );
                    }
                    else if(psiNodesf(2*iv + right)*psiNodesf(2*iv+2*right) < 0)
                    {
                        faceIntFine =  dxf + dxf * faceIntValsfFab[1](2*iv + right); //( psiNodesf(2*iv + right) / (psiNodesf(2*iv+right) - psiNodesf(2*iv+2*right) ) );
                    }
                    else
                    {
                        pout() << "face whiff at " <<iv<<endl;
                        if(psiNodesf(2*iv + 2*right)*psiNodesf(2*iv+3*right) < 0)
                        {
                            faceIntFine = 2*dxf +  dxf * ( psiNodesf(2*iv+2*right) / (psiNodesf(2*iv+2*right) - psiNodesf(2*iv+3*right) ) );
                        }
                        else if(psiNodesf(2*iv)*psiNodesf(2*iv - right) < 0)
                        {
                            faceIntFine =  -1*dxf + dxf * ( psiNodesf(2*iv - right) / (psiNodesf(2*iv-right) - psiNodesf(2*iv) ) );
                        }
                        else
                        {
                            CH_assert(false);
                        }

                    }
                    faceIntError += abs( abs(faceIntCoarse) - abs(faceIntFine) );
                    
                    
                }
                if(psiNodesc(iv)*psiNodesc(iv+up)<0)
                {
                    nFace += 1;
                    faceIntCoarse = dxc * faceIntValscFab[0](iv); //( psiNodesc(iv) / (psiNodesc(iv) - psiNodesc(iv+up) ) );
                    if(psiNodesf(2*iv)*psiNodesf(2*iv+up) < 0)
                    {
                        faceIntFine =  dxf * faceIntValsfFab[0](2*iv); //( psiNodesf(2*iv) / (psiNodesf(2*iv) - psiNodesf(2*iv+up) ) );
                    }
                    else if(psiNodesf(2*iv + up)*psiNodesf(2*iv+2*up) < 0)
                    {
                        faceIntFine =  dxf + dxf * faceIntValsfFab[0](2*iv + up); //( psiNodesf(2*iv + up) / (psiNodesf(2*iv+up) - psiNodesf(2*iv+2*up) ) );
                    }
                    else
                    {
                        pout() << "face whiff at " <<iv<<endl;
                        if(psiNodesf(2*iv + 2*up)*psiNodesf(2*iv+3*up) < 0)
                        {
                            faceIntFine = 2*dxf +  dxf * ( psiNodesf(2*iv+2*up) / (psiNodesf(2*iv+2*up) - psiNodesf(2*iv+3*up) ) );
                        }
                        else if(psiNodesf(2*iv)*psiNodesf(2*iv - up) < 0)
                        {
                            faceIntFine =  -1*dxf + dxf * ( psiNodesf(2*iv -up) / (psiNodesf(2*iv-up) - psiNodesf(2*iv) ) );
                        }
                        else
                        {
                            CH_assert(false);
                        }
                    }
                    faceIntError += abs (abs(faceIntCoarse) - abs(faceIntFine) );
                    
                    
                }
            }
            printf("faceIntError: %1.6e \n ",faceIntError/nFace);
            printf("nodeError: %1.6e \n ", nodeError/nNodes);
            
        }
    

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