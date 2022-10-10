#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBEllipticSetUp.H"
#include "NamespaceHeader.H"

/* 
create the operator in cut and irregular cells
*/


void
EBEllipticSetUp::getOperator(const LevelData<FArrayBox>& cellBeta, const LevelData<FArrayBox>& cellEta)
{
    CH_TIME("EBEllipticSetUp::getOperator");
    //this will hold duStencils for ghost cells on the edges that we can exchange to complete stencils on edge
    DisjointBoxLayout dbl ( (*m_psiNodes).disjointBoxLayout() );
    
    for(DataIterator dit = (*m_psiNodes).dataIterator(); dit.ok(); ++dit)
    {
        const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit()]).getFab();
        const IVSFAB<Real>& momentsFab = (*m_moments)[dit()];


        const IntVectSet& EBCellsFab = momentsFab.getIVS();
        
        const FArrayBox& etaFab = cellEta[dit()];
        const FArrayBox& betaFab = cellBeta[dit()];

        const IntVectSet& taggedCellsFab = ((*m_taggedFluxCells)[dit()]);

        Box valid = dbl[dit];
        valid.grow(IntVect::Unit); //so we can do calculations in the first ghost row

        IntVectSet validEB = EBCellsFab;
        validEB&=valid;

        IntVectSet validTagged = taggedCellsFab;
        validTagged&=valid;

        IntVectSet validReOp = validEB;
        
        validReOp |= validTagged; 
        
        

        IVSFAB<EBStencil>& EBOpStencilFab = (*EBOpStencil)[dit()];
        IVSFAB<EBStencil>& EBOpStencilCutFab = (*EBOpStencilCut)[dit()];
        //EBOpStencilFab.define(validReOp,2);
        //if we haven't already defined on correct IVS, define it
        bool initial = false;
        if( !(EBOpStencilFab.getIVS() == validReOp))
        {
           EBOpStencilFab.define(validReOp,2); 
           EBOpStencilCutFab.define(validEB,4);
           initial = true;
        }
        

        FArrayBox cellNodes, addStencil;
        Box stenBox, nodeStenBox, cellBox;

        IntVectSet EBCellsSten;
        LAPACKMatrix WMpW, GBeta, betaS, A;
        //GBeta.define(2,2*n); //this extracts friction integral. multiplies 2*ncoef
        

        IntVect iv, normShift, tangShift, nIv, neighIv,position;
        Interval range(0,5*n-1);
        int otherd,counter;
        Box stencilBox;


        const FluxBox& faceIntValsFab = (*m_faceIntVals)[dit()];
        FluxBox faceIntValsBox;
        IntVect right(1,0);
        IntVect up(0,1);
        Real x,y;

        
        int cx = 1*m_Q + (3 - 1)/2 + 0;
        int cy = 1;

        int l,k,row,comp; bool loGrd;

        LAPACKMatrix betaCoef, etaCoef;
        // LAPACKMatrix phiRhs(2*nStenCells,1); phiRhs.setVal(0);
        // LAPACKMatrix phiCoef(2*n,1);

        LAPACKMatrix CRhs, CJump;

        Vector<Real> fluxFaceCell, fluxEBFaceCell, fluxVect, frictionCell, frictionEBCell, frictionVect, regFaceFlux;
        Vector<Real> fluxFaceCellCRhs, fluxFaceCellCJump,fluxEBFaceCellCRhs, fluxEBFaceCellCJump,frictionEBCellCRhs, frictionEBCellCJump;
        fluxFaceCell.resize(6*2*nStenCells);
        fluxEBFaceCell.resize(2*2*nStenCells);
        frictionCell.resize(2*nStenCells);
        frictionEBCell.resize(2*2*nStenCells);

        fluxVect.resize(2*nStenCells);
        frictionVect.resize(2*nStenCells);

        regFaceFlux.resize(2*nStenCells);

        Real faceW = .5;
        if(!conserve)
        {
            faceW = 1;
        }

        //initialize and clear
        Box initBox;
        for(IVSIterator ivIt(validEB);ivIt.ok();++ivIt)
        {
            iv = ivIt();
            if(initial)
            {
                for(int c = 0;c<2;c++) //loop through the 2 phases, 
                {
                    EBStencil& currEBStencil = EBOpStencilCutFab(iv,c);
                    initBox.define(iv - (m_Q+1)*IntVect::Unit, iv + (m_Q+1)*IntVect::Unit);
                    currEBStencil.defineBox(2, initBox);
                }
            }

            for(int c = 0;c<2;c++) //loop through the 2 phases
            {
                EBStencil& currEBStencil = EBOpStencilCutFab(iv,c);
                currEBStencil.clear();
                currEBStencil.jumpConstant = 0;
                currEBStencil.rhsConstant = 0;
            }
        }

        for(IVSIterator ivIt(validTagged); ivIt.ok(); ++ivIt)
        {
            iv = ivIt();

            if(initial)
            {
                for(int c = 0;c<1;c++) //just one flux
                {
                    EBStencil& currEBStencil = EBOpStencilFab(iv,c);
                    initBox.define(iv - (m_Q+1)*IntVect::Unit, iv + (m_Q+1)*IntVect::Unit);
                    currEBStencil.defineBox(2, initBox);
                }
            }
            
            for(int c = 0;c<1;c++) //dont need this anymore
            {
                EBStencil& currEBStencil = EBOpStencilFab(iv,c);
                currEBStencil.clear();
                currEBStencil.jumpConstant = 0;
                currEBStencil.rhsConstant = 0;
            }
        }

        Real volGrd;

        Real centerX, centerY;
        LAPACKMatrix centerVal, centerM;

        for(IVSIterator ivIt(validEB);ivIt.ok();++ivIt)
        {
            iv = ivIt();

            // if(iv[0]==13 && iv[1]==40)
            // {
            //     printf("Asdfasdf \n");
            // }

            //get box from which we take neighbors
            stenBox.define(iv-r*IntVect::Unit, iv+r*IntVect::Unit);

            getCoefficients(iv, etaFab, betaFab, betaCoef, etaCoef, dit());
            //WMc = Wphi. c=[c_u^g c_u^f] phi = [u u_f]
            
            getPinvMSten(iv, betaCoef, etaCoef, dit(), WMpW, CJump); //get WM^+W for this stencil
            //WMpW only has stencil info. CRhs and CJump are vectors that we will use later to push into stencils

            counter=0;

            cellBox.define(iv,iv);
            cellBox.surroundingNodes();  cellBox.growHi(0); cellBox.growHi(1);
            cellNodes.define(cellBox,1);
            cellNodes.copy(psiNodesFab);

            getFluxFaceCell(WMpW, cellNodes, faceIntValsFab, etaCoef, fluxFaceCell,CJump,fluxFaceCellCJump ); //6 vectors of length 2*nStenCells: L_x, L_y

            getFluxEBFaceCell(WMpW, iv, dit(), etaCoef, fluxEBFaceCell,CJump, fluxEBFaceCellCJump); //2*2 vectors of length 2*nStenCells

            getFrictionEBCell(WMpW, iv, dit(), betaCoef, frictionEBCell,CJump, frictionEBCellCJump); //2 vectors of length 2*nStenCells

            if(psiNodesFab(iv) > 0){volGrd = momentsFab(iv,0);}
            else{volGrd = momentsFab(iv,hiVol);}
            
            counter = 0;
            for(int i=0;i<2;i++) //for phases
            {
                EBStencil& currEBStencil = EBOpStencilCutFab(iv, i); //phase 
                
                for(int j=0;j<2*nStenCells;j++)
                {
                    fluxVect[j] = fluxEBFaceCell[counter + j];
                }
                
                currEBStencil.addVect(stenBox,fluxVect,1,initial);
                currEBStencil.jumpConstant += fluxEBFaceCellCJump[i];

                for(int j=0;j<2*nStenCells;j++)
                {
                    frictionVect[j] = frictionEBCell[counter + j];
                }

                counter += 2*nStenCells;

                currEBStencil.addVect(stenBox, frictionVect, -1*m_dx*m_dx, initial);

                currEBStencil.jumpConstant += -1*m_dx*m_dx*frictionEBCellCJump[i];
                            
            }
            assert(counter == 2*2*nStenCells);

            //loop through the cell faces and build the stencil the flux stencil FOR THIS CELL AND NEIGHBOR
            
            counter = 0; //vectors
            for(int d = 0;d<2;d++) //loop through dims
            {
                otherd = (d+1)%2;
                normShift = IntVect::Zero;
                tangShift = IntVect::Unit;
                normShift[d] = 1;
                tangShift[d] = 0;
                for(int s = 0;s<2;s++) //loop through sides
                {
                    nIv = iv + s*normShift;
                    neighIv = nIv - ((s+1)%2)*normShift; //neighboring normal center

                    if(psiNodesFab(nIv) * psiNodesFab(nIv+tangShift) < 0) 
                    {
                        l=2;
                    }
                    else
                    {
                        l=1;
                    }

                    loGrd = false;
                    if(psiNodesFab(nIv) > 0) //lo side of face is grd
                    {
                        loGrd = true;
                    }


                    for(int i =0;i<l;i++) //loop through lo and hi side of FACE if cut
                    {
                        if(loGrd){k = i;}
                        else{k = (i+1)%2;}

                        
                        EBStencil& currEBStencil = EBOpStencilCutFab(iv,k ); //grd grd flt flt

                        for(int j=0;j<2*nStenCells;j++)
                        {
                            fluxVect[j] = fluxFaceCell[counter + j];
                        }

                        currEBStencil.jumpConstant += faceW*fluxFaceCellCJump[counter/(2*nStenCells)];
                        
                        currEBStencil.addVect(stenBox, fluxVect, faceW, initial);

                        if(conserve)
                        {
                            if(validEB.contains(neighIv) )//neighbor is fellow cut cell
                            {
                                EBStencil& currNeighEBStencil = EBOpStencilCutFab(neighIv,k); //grd flt
                                currNeighEBStencil.addVect(stenBox, fluxVect, -.5, initial); //normal is going the other way
                                
                                currNeighEBStencil.jumpConstant += -.5*fluxFaceCellCJump[counter/(2*nStenCells)];

                                counter += 2*nStenCells;
                            }
                            else if(validTagged.contains(neighIv))//neighbor is a tagged cell
                            {
                                EBStencil& currNeighEBStencil = EBOpStencilFab(neighIv,0); //grd grd flt flt
                                currNeighEBStencil.addVect(stenBox, fluxVect, -.5, initial); //normal is going the other way

                                currNeighEBStencil.jumpConstant += -.5*fluxFaceCellCJump[counter/(2*nStenCells)];

                                counter += 2*nStenCells;
                            }
                            else
                            { 
                                counter += 2*nStenCells;
                                continue;
                            }
                        }

                        else{counter+=2*nStenCells;}
                    }

                }
            }
            
            CH_assert(counter == (6*nStenCells*2));
  
        }

        Vector<Real> frictionCellCRhs;
        //Vector<Real> fluxFaceCellCRhs;
        for(IVSIterator ivIt(validTagged); ivIt.ok(); ++ivIt)
        {
            iv = ivIt();
            // if(iv[0]==33 && iv[1]==43)
            // {
            //     printf("asdfa \n");
            // }

            getCoefficientsFull(iv, etaFab, betaFab, betaCoef, etaCoef, dit());

            getPinvMStenFull(iv, betaCoef, etaCoef, dit(), WMpW); //get WM^+W for this stencil

            cellBox.define(iv,iv);
            cellBox.surroundingNodes();  cellBox.growHi(0); cellBox.growHi(1);
            cellNodes.define(cellBox,1);
            cellNodes.copy(psiNodesFab);

            stenBox.define(iv-r*IntVect::Unit, iv + r*IntVect::Unit);

            //get du stencil on all faces of this cell
            getFluxFaceCellFull(WMpW, etaCoef, cellNodes, fluxFaceCell); //2*4 vectors of length 2*nStenCells

            
            getFrictionCell(WMpW, iv, dit(), betaCoef, frictionCell);
            
            EBStencil& currEBStencil = EBOpStencilFab(iv,0); //grd grd flt flt

            for(int j=0;j<2*nStenCells;j++)
            {
                frictionVect[j] = frictionCell[j];
            }

            currEBStencil.addVect(stenBox, frictionVect, -1*m_dx*m_dx, initial);
            
            //loop through the cell faces and build the stencil for flux
            
            counter =0; 
            for(int d = 0;d<2;d++)
            {
                normShift = IntVect::Zero;
                tangShift = IntVect::Unit;
                normShift[d] = 1;
                tangShift[d] = 0;
                for(int s = 0;s<2;s++)
                {
                    nIv = iv + s*normShift;
                    neighIv = nIv - ((s+1)%2)*normShift; //neighboring normal center
                    if( ( (*(m_taggedFaces[d]) )[dit()] ).contains(nIv)  )
                    {
                        
                        EBStencil& currEBStencil = EBOpStencilFab(iv,0); //grd grd flt flt

                        for(int j=0;j<2*nStenCells;j++)
                        {
                            fluxVect[j] = fluxFaceCell[counter + j];
                            // if(isnan(fluxVect[j]))
                            // {
                            //     printf("die \n");
                            // }
                        }
                        
                        currEBStencil.addVect(stenBox, fluxVect, faceW, initial);
                        if(conserve)
                        {
                            if(validEB.contains(neighIv) )//neighbor is a cut cell
                            {
                                if(psiNodesFab(nIv) > 0) 
                                {
                                    EBStencil& currNeighEBStencil = EBOpStencilCutFab(neighIv,0); //grd grd flt flt
                                    currNeighEBStencil.addVect(stenBox, fluxVect, -.5, initial); //normal is going the other way 
                                   
                                }
                                else
                                {
                                    EBStencil& currNeighEBStencil = EBOpStencilCutFab(neighIv,1); //grd grd flt flt
                                    currNeighEBStencil.addVect(stenBox, fluxVect,-.5, initial); //normal is going the other way  
                    
                                }
                                
                            }
                            else if(validTagged.contains(neighIv))
                            {
                                EBStencil& currNeighEBStencil = EBOpStencilFab(neighIv,0); //grd grd flt flt
                                currNeighEBStencil.addVect(stenBox, fluxVect, -.5, initial); //normal is going the other way
                                
                            }
                            else
                            {
                                counter+=2*nStenCells;
                                continue;
                            }
                        }
                        counter += 2*nStenCells;

                        
                        

                    }
                    else //not tagged face -- need to add regular face to this boyo
                    {
                        if(conserve)
                        {
                            getRegFaceStencil(s,d,etaFab,iv,dit(),regFaceFlux);
                        }
                        
                        
                        if(conserve)
                        {
                            for(int j=0;j<2*nStenCells;j++)
                            {
                                fluxVect[j] = regFaceFlux[j];
                                //fluxVect[j] = fluxFaceCell[counter + j];
                            }
                        }
                        else
                        {
                            for(int j=0;j<2*nStenCells;j++)
                            {
                                fluxVect[j] = fluxFaceCell[counter + j];
                            }
                        }
                        

                        EBStencil& currEBStencil = EBOpStencilFab(iv,0);

                        currEBStencil.addVect(stenBox, fluxVect, 1, initial);

                        counter += 2*nStenCells;

                    }
    
                }
            }
            CH_assert(counter == (4*2*nStenCells));

        }
              
    }


    
}

void 
EBEllipticSetUp::getCoefficients(const IntVect center, const FArrayBox& etaFab, const FArrayBox& betaFab,
LAPACKMatrix& betaCoef, LAPACKMatrix& etaCoef,
                        const DataIndex& dit)
{
    //interpolate local (point) values of beta and eta on each side of the EB to get coefficients
    const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit]).getFab();
    const IVSFAB<Real>& momentsFab = (*m_moments)[dit];

    const IntVectSet& EBCellsFab = momentsFab.getIVS();
    IntVectSet EBCellsSten;
    EBCellsSten.define(EBCellsFab);
    Box stenBox(center - r*IntVect::Unit, center + r*IntVect::Unit);
    EBCellsSten &= stenBox;
    int nCutCells = EBCellsSten.numPts();

    int m = nCutCells + nStenCells; //bc we have one val in full cells, 2 in cut cells
    

    LAPACKMatrix M(m, 2*n); M.setVal(0); //to interoplate points - grd then floating
    LAPACKMatrix phiEta(m,1); phiEta.setVal(0);
    LAPACKMatrix phiBeta(m,1); phiBeta.setVal(0); //full, full, lo, hi  etc

    LAPACKMatrix W(m,m); W.setVal(0);

    IntVect iv,position;

    int cx = 1*m_Q + (3 - 1)/2;
    int cy = 1;
    int row=0;
    Real loCtdX, loCtdY, hiCtdX, hiCtdY;
    Real d,w; int comp;

    Vector<Real> volMomentsLo(n,0);
    Vector<Real> volMomentsHi(n,0);

    for(int l = -r;l<r+1;l++)
    {
        for(int c = -r;c<r+1;c++) //c++ hahhaa
        {
            position[0] = c; position[1]= l;
            iv = center + position;
            d = sqrt(position[0]*position[0] + position[1]*position[1]) + 1;
            w = pow(d, -1*m_weight);
            

            //find centroid of this cell and calculate point moments
            if(EBCellsFab.contains(iv)) //its a cut cell
            {
                for(int i =0;i<n;i++)
                {
                    volMomentsLo[i] = momentsFab(iv,i);
                    volMomentsHi[i] = momentsFab(iv,i+hiVol);
                }

                shiftMoments(volMomentsLo, m_Q, position[0], position[1]); 
                shiftMoments(volMomentsHi, m_Q, position[0], position[1]);

                loCtdX = volMomentsLo[cx] / volMomentsLo[0];
                loCtdY = volMomentsLo[cy] / volMomentsLo[0];

                hiCtdX = volMomentsHi[cx] / volMomentsHi[0];
                hiCtdY = volMomentsHi[cy] / volMomentsHi[0];

                if(psiNodesFab(iv) > 0)
                {
                    for(int px =0; px<=m_Q; px++)
                    {
                        for(int py = 0; py<=m_Q-px;py++)
                        {
                            comp =  px*m_Q + (3*px - px*px)/2 + py;
                            M(row,comp) = w*pow(loCtdX,px) * pow(loCtdY, py);
                            M(row + 1,comp + n) = w*pow(hiCtdX,px) * pow(hiCtdY, py);
                            
                        }
                    }
                }
                else
                {
                    for(int px =0; px<=m_Q; px++)
                    {
                        for(int py = 0; py<=m_Q-px;py++)
                        {
                            comp =  px*m_Q + (3*px - px*px)/2 + py;
                            M(row,comp) = w*pow(hiCtdX,px) * pow(hiCtdY, py);
                            M(row + 1,comp+n) = w*pow(loCtdX,px) * pow(loCtdY, py);
                            
                        }
                    }
                }

                phiBeta(row,0) = w*betaFab(iv,0);
                phiBeta(row+1,0) = w*betaFab(iv,1);

                phiEta(row,0) = w*etaFab(iv,0);
                phiEta(row+1,0) = w*etaFab(iv,1);

                row +=2;

                
            }
            else //full cell
            {
                loCtdX = position[0];
                loCtdY = position[1];
                if(psiNodesFab(iv) > 0)
                {
                    for(int px =0; px<=m_Q; px++)
                    {
                        for(int py = 0; py<=m_Q-px;py++)
                        {
                            comp =  px*m_Q + (3*px - px*px)/2 + py;
                            M(row,comp) = w*pow(loCtdX,px) * pow(loCtdY, py);
                        }
                    }
                    phiBeta(row,0) = w*betaFab(iv,0);

                    phiEta(row,0) = w*etaFab(iv,0);
                }
                else
                {
                    for(int px =0; px<=m_Q; px++)
                    {
                        for(int py = 0; py<=m_Q-px;py++)
                        {
                            comp =  px*m_Q + (3*px - px*px)/2 + py;
                            M(row,comp + n) = w*pow(loCtdX,px) * pow(loCtdY, py);
                            
                        }
                    }

                    phiBeta(row,0) = w*betaFab(iv,1);

                    phiEta(row,0) = w*etaFab(iv,1);
                }

                

                row +=1;
            }
            

        }
    }

    pseudoInvert(M);
    multiply(betaCoef, M, phiBeta);
    multiply(etaCoef, M, phiEta);

    Real tol = 1e-13;
    for(int i=0;i<2*n;i++)
    {
        if( abs(betaCoef(i,0) ) < tol)
        {
            betaCoef(i,0) = 0;
        }
        if( abs(etaCoef(i,0) ) < tol)
        {
            etaCoef(i,0) = 0;
        }
    }
}

void 
EBEllipticSetUp::getCoefficientsFull(const IntVect center, const FArrayBox& etaFab, const FArrayBox& betaFab,
LAPACKMatrix& betaCoef, LAPACKMatrix& etaCoef,
                        const DataIndex& dit)
{
    //interpolate local (point) values of beta and eta on each side of the EB to get coefficients
    const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit]).getFab();
    const IVSFAB<Real>& momentsFab = (*m_moments)[dit];

    const IntVectSet& EBCellsFab = momentsFab.getIVS();
    IntVectSet EBCellsSten;
    EBCellsSten.define(EBCellsFab);
    Box stenBox(center - r*IntVect::Unit, center + r*IntVect::Unit);
    EBCellsSten &= stenBox;
    int nCutCells = EBCellsSten.numPts();

    int m = nCutCells + nStenCells; //bc we have one val in full cells, 2 in cut cells
    

    LAPACKMatrix M(m, n); M.setVal(0); //to interoplate points - grd then floating
    LAPACKMatrix phiEta(m,1); phiEta.setVal(0);
    LAPACKMatrix phiBeta(m,1); phiBeta.setVal(0);  //full, full, lo, hi  etc

    LAPACKMatrix W(m,m); W.setVal(0);

    IntVect iv,position;

    int cx = 1*m_Q + (3 - 1)/2;
    int cy = 1;
    int row=0;
    Real ctdX, ctdY;
    Real d,w; int comp;

    Vector<Real> volMoments(n,0);

    bool grd = (psiNodesFab(center) > 0);

    for(int l = -r;l<r+1;l++)
    {
        for(int c = -r;c<r+1;c++) //c++ hahhaa
        {
            position[0] = c; position[1]= l;
            iv = center + position;
            d = sqrt(position[0]*position[0] + position[1]*position[1]) + 1;
            w = pow(d, -1*m_weight);
            

            //find centroid of this cell and calculate point moments
            if(EBCellsFab.contains(iv)) //its a cut cell, use the correct side
            {
                if( (psiNodesFab(iv)>0) ==grd ) //use lo side
                {
                    for(int i =0;i<n;i++)
                    {
                        volMoments[i] = momentsFab(iv,i);
                    }

                }
                else
                {
                    for(int i =0;i<n;i++)
                    {
                        volMoments[i] = momentsFab(iv,i + hiVol);
                    }
                }
                

                shiftMoments(volMoments, m_Q, position[0], position[1]); 

                ctdX = volMoments[cx] / volMoments[0];
                ctdY = volMoments[cy] / volMoments[0];


                
                for(int px =0; px<=m_Q; px++)
                {
                    for(int py = 0; py<=m_Q-px;py++)
                    {
                        comp =  px*m_Q + (3*px - px*px)/2 + py;
                        M(row,comp) = w*pow(ctdX,px) * pow(ctdY, py);
                    }
                }
         

                
                if(grd)
                {
                    phiEta(row,0) = w*etaFab(iv,0);
                    phiBeta(row,0) = w*betaFab(iv,0);
                }
                else
                {
                    phiEta(row,0) = w*etaFab(iv,1);
                    phiBeta(row,0) = w*betaFab(iv,1);
                }
                
                row +=1;

            }
            else //full cell
            {
                ctdX = position[0];
                ctdY = position[1];
                if( (psiNodesFab(iv) > 0) == grd) //if same type of cell
                {
                    for(int px =0; px<=m_Q; px++)
                    {
                        for(int py = 0; py<=m_Q-px;py++)
                        {
                            comp =  px*m_Q + (3*px - px*px)/2 + py;
                            M(row,comp) = w*pow(ctdX,px) * pow(ctdY, py);
                        }
                    }
                    if(grd)
                    {
                        phiEta(row,0) = w*etaFab(iv,0);
                        phiBeta(row,0) = w*betaFab(iv,0);
                    }
                    else
                    {
                        phiEta(row,0) = w*etaFab(iv,1);
                        phiBeta(row,0) = w*betaFab(iv,1);
                    }
                }
                row +=1;
            }
            

        }
    }

    pseudoInvert(M);
    multiply(betaCoef, M, phiBeta);
    multiply(etaCoef, M, phiEta);
    Real tol = 1e-13;
    for(int i=0;i<n;i++)
    {
        if( abs(betaCoef(i,0) ) < tol)
        {
            betaCoef(i,0) = 0;
        }
        if( abs(etaCoef(i,0) ) < tol)
        {
            etaCoef(i,0) = 0;
        }
    }
}


void 
EBEllipticSetUp::getFluxFaceCell(const LAPACKMatrix& WMpW, const FArrayBox& cellNodes, 
                         const FluxBox& faceIntValsFab, const LAPACKMatrix& etaCoef, 
                         Vector<Real>& fluxFaceCell, const LAPACKMatrix& CJump,
                         Vector<Real>& fluxFaceCellCJump)
{
    fluxFaceCell.assign(0);
    int counter=0;//where you are in the vector
    Real loPart, hiPart, loCtd, hiCtd;
    IntVect normShift, tangShift, iv;
    IntVect ll = (cellNodes.box()).smallEnd();
    Real sgn,x1,x2,x3,y1,y2,y3,momXLo,momYLo,momXHi,momYHi,nx,ny; int comp;
    Real loMid, hiMid;
    int row;
    int Q = m_Q;

    LAPACKMatrix MFULo, MFUHi, MF;
    LAPACKMatrix MFUGrd;
    LAPACKMatrix MFUFlt, SF;
    MF.define(4,2*n);

    LAPACKMatrix etaCoefGrd(n,1);
    LAPACKMatrix etaCoefFlt(n,1);  
    for(int i=0;i<n;i++)
    {
        etaCoefGrd(i,0) = etaCoef(i,0);
        etaCoefFlt(i,0) = etaCoef(i+n,0);
    }
    etaCoefGrd.transpose(); etaCoefFlt.transpose();

    Vector<Real> bds(2,0);
    Vector<Real> bdsHi(2,0);
    bool loGrd;

    fluxFaceCellCJump.resize(6);

    int cV = 0;

    LAPACKMatrix SFCRhs, SFCJump;

    for(int d = 0;d<2;d++) //directions
    {
        normShift = IntVect::Zero;
        tangShift = IntVect::Unit;
        normShift[d] = 1;
        tangShift[d] = 0;
        for(int s = 0;s<2;s++) //sides of cell
        {
            sgn = (s==0) ? -1 : 1;
            nx = ((d+1)%2)*sgn;//if d is 0, nx is 1
            ny = d*sgn; //if d is 0, ny is 0

            MFULo.define(n, n); MFULo.setVal(0);

            MFUHi.define(n, n); MFUHi.setVal(0);

            iv = ll + s*normShift;

            int l;

            if(cellNodes(iv)*cellNodes(iv+tangShift) < 0) //this is a cut face
            {
                loPart = faceIntValsFab[d](iv);
                l=2;
            }
            else
            {
                l=1;
            }
            loGrd = false;
            if(cellNodes(iv) > 0){loGrd = true;}

            if(cellNodes(iv)*cellNodes(iv+tangShift) < 0) //this is a cut face
            {
                MF.define(2,2*n); MF.setVal(0);
                bds[0]=-.5; bds[1]=loPart-.5;
                bdsHi[0]=loPart-.5;  bdsHi[1]=.5;
                row=0;
                for(int qx = 0; qx<=Q; qx++) //i know youre not supposed to do loops like this but whatever these numbers are small ¯\_(ツ)_/¯ 
                {
                    for(int qy = 0;qy <= Q - qx; qy++)
                    {
                        for(int px = 0; px<=Q; px++)
                        {
                            for(int py =0; py<= Q -px; py++)
                            {
                                comp =  (px)*Q + (3*(px) - (px)*(px)) /2 + (py);

                                MFULo(row,comp) = qx * faceMom(px+qx-1,py+qy,d,s,bds) * nx + qy * faceMom(px+qx,py+qy-1,d,s,bds) * ny;

                                MFUHi(row,comp) = qx * faceMom(px+qx-1,py+qy,d,s,bdsHi) * nx + qy * faceMom(px+qx,py+qy-1,d,s,bdsHi) * ny;

                                
                            }
                        }
                        row+=1;
                    }
                }
                
                MFULo.transpose();

                MFUHi.transpose();

                if(loGrd)
                {
                    multiply(MFUGrd, etaCoefGrd, MFULo);
                    multiply(MFUFlt, etaCoefFlt, MFUHi);

                    for(int i = 0;i<n;i++)
                    {
                        MF(0,i) = MFUGrd(0,i); 

                        MF(1,n+i) = MFUFlt(0,i);

                    }

                }
                else
                {
                    multiply(MFUGrd, etaCoefGrd, MFUHi);
                    multiply(MFUFlt, etaCoefFlt, MFULo);

                    for(int i = 0;i<n;i++)
                    {
                        MF(0, n + i) = MFUFlt(0,i); //lo side use floating coef

                        MF(1, i) = MFUGrd(0,i);

                    }
                }

                multiply(SF, MF, WMpW);

               
                multiply(SFCJump, MF, CJump);

                for(int i=0;i<2;i++)
                {
                    fluxFaceCellCJump[cV+i] = SFCJump(i,0);
                }

                cV+=2;

                for(int i=0;i<2*nStenCells;i++)
                {
                    fluxFaceCell[counter + i] = SF(0, i); //L(u) lo
                    fluxFaceCell[counter + i + 2*nStenCells] = SF(1, i); //L(u) Hi

                }
                counter += 2*2*nStenCells;


            }

            else
            {
                row=0;
                MF.define(1,2*n); MF.setVal(0);
                bds[0]=-.5; bds[1]=.5;
                for(int qx = 0; qx<=Q; qx++) //i know youre not supposed to do loops like this but whatever these numbers are small ¯\_(ツ)_/¯ 
                {
                    for(int qy = 0;qy <= Q - qx; qy++)
                    {
                        for(int px = 0; px<=Q; px++)
                        {
                            for(int py =0; py<= Q -px; py++)
                            {
                                comp =  (px)*Q + (3*(px) - (px)*(px)) /2 + (py);

                                MFULo(row,comp) = qx * faceMom(px+qx-1,py+qy,d,s,bds) * nx + qy * faceMom(px+qx,py+qy-1,d,s,bds) * ny;
                            }
                        }
                        row+=1;
                    }
                }
                
                MFULo.transpose();

                if(loGrd)
                {
                    multiply(MFUGrd, etaCoefGrd, MFULo);

                    for(int i = 0;i<n;i++)
                    {
                        MF(0,i) = MFUGrd(0,i); 
                    }

                }
                else
                {
                    multiply(MFUFlt, etaCoefFlt, MFULo);

                    for(int i = 0;i<n;i++)
                    {
                        MF(0,n+i) = MFUFlt(0,i); 
                    }
                }

                multiply(SF, MF, WMpW);


                multiply(SFCJump, MF, CJump);

                for(int i=0;i<1;i++)
                {
                    fluxFaceCellCJump[cV+i] = SFCJump(i,0);
                }

                cV+=1;

                for(int i=0;i<2*nStenCells;i++)
                {
                    fluxFaceCell[counter + i] = SF(0, i);
                }
                counter += 2*nStenCells;
            }
        }
    }
    CH_assert(cV == 6);
}

void 
EBEllipticSetUp::getFluxEBFaceCell(const LAPACKMatrix& WMpW, const IntVect iv, const DataIndex& dit,
                        const LAPACKMatrix& etaCoef, Vector<Real>& fluxEBFaceCell,
                        const LAPACKMatrix& CJump,Vector<Real>& fluxEBFaceCellCJump)
{
    //need to use EB moments, but this is similar to getFluxFaceCell
    //should probably make another function to create the stencil for EB flux as we will 
    //also use it in getPinMSten :) 

    //returns Fx grd, Fy grd, Fx flt, Fy flt
    const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit]).getFab();
    const IVSFAB<Real>& momentsFab = (*m_moments)[dit];
    bool loGrd = (psiNodesFab(iv) > 0);
    int grdSgn=-1; int fltSgn=1;
    if(loGrd){grdSgn = 1; fltSgn = -1;}

    fluxEBFaceCell.assign(0);

    Vector<Real> nxMoments(n,0);
    Vector<Real> nyMoments(n,0);

    for(int i=0;i<n;i++)
    {
        nxMoments[i] = momentsFab(iv,nX + i);
        nyMoments[i] = momentsFab(iv,nY + i);
    }

    //etaCoef has grd then floating coefs

    LAPACKMatrix MFGrd, MFFlt,  FGrdSten, FFltSten ;

    LAPACKMatrix  FGrdCJump, FFltCJump;

    getEBFluxMoments(nxMoments, nyMoments, etaCoef, MFGrd, MFFlt );

    multiply(FGrdSten, MFGrd, WMpW);

    multiply(FFltSten, MFFlt, WMpW);


    fluxEBFaceCellCJump.resize(2);
    

    multiply(FGrdCJump, MFGrd, CJump); 
    multiply(FFltCJump, MFFlt, CJump); 

    fluxEBFaceCellCJump[0] = grdSgn*FGrdCJump(0,0); fluxEBFaceCellCJump[1] = fltSgn*FFltCJump(0,0);


    for(int i=0;i<2*nStenCells;i++)
    {
        fluxEBFaceCell[i] = grdSgn*FGrdSten(0,i);
        fluxEBFaceCell[i+ 2*nStenCells] = fltSgn*FFltSten(0,i);

    }

}

 void 
 EBEllipticSetUp::getEBFluxMoments(const Vector<Real>& nxMoments, const Vector<Real>& nyMoments,
                            const LAPACKMatrix& etaCoef, LAPACKMatrix& MFGrd,LAPACKMatrix& MFFlt)
{
    LAPACKMatrix MF;
    LAPACKMatrix MFUGrd; 
    LAPACKMatrix MFUFlt;
    
    MF.define(n, n); MF.setVal(0);
    MFUGrd.define(n, n); MFUGrd.setVal(0);
    MFUFlt.define(n, n); MFUFlt.setVal(0);

    int row=0;
    int comp, compy, compx;
    int Q = m_Q;

    Real nXx, nYx, nXy, nYy;

    LAPACKMatrix etaCoefGrd(n,1);
    LAPACKMatrix etaCoefFlt(n,1);  
    for(int i=0;i<n;i++)
    {
        etaCoefGrd(i,0) = etaCoef(i,0);
        etaCoefFlt(i,0) = etaCoef(i+n,0);
    }
    etaCoefGrd.transpose(); etaCoefFlt.transpose();


    
    for(int qx = 0; qx<=Q; qx++) 
    {
        for(int qy = 0;qy <= Q - qx; qy++)
        {
            for(int px = 0; px<=Q; px++)
            {
                for(int py =0; py<= Q -px; py++)
                {
                    comp =  (px)*Q + (3*(px) - (px)*(px)) /2 + (py);
                    compx = (px+qx-1)*Q + (3*(px+qx-1) - (px+qx-1)*(px+qx-1)) /2 + (py+qy);
                    compy = (px+qx)*Q + (3*(px+qx) - (px+qx)*(px+qx)) /2 + (py+qy-1);

                    if(px + qx + py + qy - 1 > Q)
                    {
                        continue;
                    } //we dont need this high order junk

                    nXx = (qx>=1) ? nxMoments[compx]: 0;

                    nYy = (qy>=1) ? nyMoments[compy]: 0;

                    MF(row,comp) = qx * nXx + qy * nYy;
                }
            }
            row+=1;
        }
    }
    
    MF.transpose();

    multiply(MFUGrd, etaCoefGrd, MF);
    multiply(MFUFlt, etaCoefFlt, MF);

    MFGrd.define(1,2*n); MFGrd.setVal(0);
    MFFlt.define(1,2*n); MFFlt.setVal(0);
    

    for(int i = 0;i<n;i++)
    {
        MFGrd(0,i) = MFUGrd(0,i); 

        MFFlt(0,n+i) = MFUFlt(0,i); 

    }


}

void 
EBEllipticSetUp::getFrictionEBCell(const LAPACKMatrix& WMpW, IntVect iv, const DataIndex& dit,
                        const LAPACKMatrix& betaCoef, Vector<Real>& frictionEBCell,
                        const LAPACKMatrix& CJump, Vector<Real>& frictionEBCellCJump)
{
    //this should be pretty simple, just create Mbeta
    const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit]).getFab();
    const IVSFAB<Real>& momentsFab = (*m_moments)[dit];
    bool loGrd = psiNodesFab(iv) > 0;
    LAPACKMatrix betaCoefGrd(n,1);
    LAPACKMatrix betaCoefFlt(n,1);
    int Q = m_Q;
    LAPACKMatrix MBeta2Grd, MBeta2Flt;
    for(int i=0;i<n;i++)
    {
        betaCoefGrd(i,0) = betaCoef(i,0);
        betaCoefFlt(i,0) = betaCoef(i+n,0);
    }

    Vector<Real> momentsGrd(n,0);
    Vector<Real> momentsFlt(n,0);

    if(loGrd)
    {
        for(int i=0;i<n;i++)
        {
           momentsGrd[i] = momentsFab(iv,i);
           momentsFlt[i] = momentsFab(iv,i + hiVol);
        }
    }
    else
    {
        for(int i=0;i<n;i++)
        {
           momentsFlt[i] = momentsFab(iv,i);
           momentsGrd[i] = momentsFab(iv,i + hiVol);
        }
    }

    LAPACKMatrix MBetaGrd(n,n); MBetaGrd.setVal(0);
    LAPACKMatrix MBetaFlt(n,n); MBetaFlt.setVal(0);
    int row=0;
    int pxx,pyy;
    int comp, compP;

    for(int qx = 0; qx<=Q; qx++) //i know youre not supposed to do loops like this but whatever these numbers are small ¯\_(ツ)_/¯ 
    {
        for(int qy = 0;qy <= Q - qx; qy++)
        {
            for(int px = 0; px<=Q; px++)
            {
                for(int py =0; py<= Q -px; py++)
                {
                    pxx = px+qx;
                    pyy = py+qy;
                    if(pxx + pyy > Q){continue;}
                    comp=  (px)*Q + (3*(px) - (px)*(px)) /2 + (py);
                    compP = (pxx)*Q + (3*(pxx) - (pxx)*(pxx)) /2 + (pyy);
                    MBetaGrd(row,comp) = momentsGrd[compP]; //column corresponds to coefficient of beta, whcih has p
                    MBetaFlt(row,comp) = momentsFlt[compP];
                }
            }
            row+=1;
        }
    }

    
    LAPACKMatrix MBetaU;
    MBetaU.define(2, 2*n); MBetaU.setVal(0);

    multiply(MBeta2Grd, MBetaGrd, betaCoefGrd );
    multiply(MBeta2Flt, MBetaFlt, betaCoefFlt );

    MBeta2Grd.transpose();
    MBeta2Flt.transpose();

    for(int i=0;i<n;i++)
    {
        MBetaU(0,i) = MBeta2Grd(0,i);
        MBetaU(1,i+n) = MBeta2Flt(0,i);
    }

    LAPACKMatrix betaStenU, betaStenV;
    multiply(betaStenU, MBetaU, WMpW);

    frictionEBCellCJump.resize(2);

    LAPACKMatrix betaCJumpU;

    multiply(betaCJumpU, MBetaU, CJump);

    frictionEBCellCJump[0] = betaCJumpU(0,0);
    frictionEBCellCJump[1] = betaCJumpU(1,0);


    for(int i=0;i<2*nStenCells;i++)
    {
        frictionEBCell[i] = betaStenU(0,i);
        frictionEBCell[i + 2*nStenCells] = betaStenU(1,i);
    }


}

//check
void 
EBEllipticSetUp::getFluxFaceCellFull(const LAPACKMatrix& WMpW, const LAPACKMatrix& etaCoef, 
                        const FArrayBox& cellNodes, Vector<Real>& fluxFaceCell)
{
    fluxFaceCell.assign(0);
    int counter=0;//where you are in the vector
    Real loPart, hiPart, loCtd, hiCtd;
    IntVect normShift, tangShift, iv;
    IntVect ll = (cellNodes.box()).smallEnd();
    Real sgn,x1,x2,x3,y1,y2,y3,momXLo,momYLo,momXHi,momYHi,nx,ny; int comp;
    Real loMid, hiMid;
    int row;
    int Q = m_Q;

    LAPACKMatrix MF;
    LAPACKMatrix MFU;
    LAPACKMatrix SF;
    MF.define(4,2*n);

    LAPACKMatrix etaCoefV(n,1);
    for(int i=0;i<n;i++)
    {
        etaCoefV(i,0) = etaCoef(i,0);
    }
    etaCoefV.transpose(); 

    Vector<Real> bds(2,0);
    Vector<Real> bdsHi(2,0);
    bool loGrd;


    for(int d = 0;d<2;d++) //directions
    {
        normShift = IntVect::Zero;
        tangShift = IntVect::Unit;
        normShift[d] = 1;
        tangShift[d] = 0;
        for(int s = 0;s<2;s++) //sides of cell
        {
            sgn = (s==0) ? -1 : 1;
            nx = ((d+1)%2)*sgn;//if d is 0, nx is 1
            ny = d*sgn; //if d is 0, ny is 0

            MF.define(n, n); MF.setVal(0);

            iv = ll + s*normShift;

            loGrd = false;
            if(cellNodes(iv) > 0){loGrd = true;}

            row=0;
            //MF.define(2,2*n); MF.setVal(0);
            bds[0]=-.5; bds[1]=.5;
            for(int qx = 0; qx<=Q; qx++) //i know youre not supposed to do loops like this but whatever these numbers are small ¯\_(ツ)_/¯ 
            {
                for(int qy = 0;qy <= Q - qx; qy++)
                {
                    for(int px = 0; px<=Q; px++)
                    {
                        for(int py =0; py<= Q -px; py++)
                        {
                            comp =  (px)*Q + (3*(px) - (px)*(px)) /2 + (py);

                            MF(row,comp) = qx * faceMom(px+qx-1,py+qy,d,s,bds) * nx + qy * faceMom(px+qx,py+qy-1,d,s,bds) * ny;
                        }
                    }
                    row+=1;
                }
            }
            
            MF.transpose();

            multiply(MFU, etaCoefV, MF);

            multiply(SF, MFU, WMpW);

            for(int i=0;i<2*nStenCells;i++)
            {
                fluxFaceCell[counter + i] = SF(0, i);
            }
            counter += 2*nStenCells;
            
        }
    }
}
//check
void 
EBEllipticSetUp::getFrictionCell(const LAPACKMatrix& WMpW, IntVect iv, const DataIndex& dit,
                        const LAPACKMatrix& betaCoef, Vector<Real>& frictionCell)

{
    int row=0;
    int comp;
    Real mom;
    LAPACKMatrix MBeta(n,n); MBeta.setVal(0);
    LAPACKMatrix MBeta2;
    int pxx,pyy;
    int Q = m_Q;
    for(int qx = 0; qx<=Q; qx++) //i know youre not supposed to do loops like this but whatever these numbers are small ¯\_(ツ)_/¯ 
    {
        for(int qy = 0;qy <= Q - qx; qy++)
        {
            for(int px = 0; px<=Q; px++)
            {
                for(int py =0; py<= Q -px; py++)
                {
                    pxx = px+qx;
                    pyy = py+qy;
                    comp=  (px)*Q + (3*(px) - (px)*(px)) /2 + (py);
                    mom = (pow(0+.5,pxx+1) - pow(0-.5,pxx+1) ) / (pxx+1);
                    mom*= (pow(0+.5,pyy+1) - pow(0-.5,pyy+1) ) / (pyy+1);
                    MBeta(row,comp) = mom; //column corresponds to coefficient of beta, whcih has p
                }
            }
            row+=1;
        }
    }

    LAPACKMatrix MBetaU, MBetaV;
    MBetaU.define(1, n); MBetaU.setVal(0);

    multiply(MBeta2, MBeta, betaCoef );

    MBeta2.transpose();

    for(int i=0;i<n;i++)
    {
        MBetaU(0,i) = MBeta2(0,i);
    }

    LAPACKMatrix betaStenU, betaStenV, betaCRhsU, betaCRhsV;
    multiply(betaStenU, MBetaU, WMpW);


    for(int i=0;i<2*nStenCells;i++)
    {
        frictionCell[i] = betaStenU(0,i);
    }
}

void 
EBEllipticSetUp::getRegFaceStencil(int s, int d, const FArrayBox& etaFab, 
                            const IntVect iv, const DataIndex& dit,Vector<Real>& regFaceFlux )
{
    regFaceFlux.assign(0);
    Real sgn = s==0 ? -1:1; //for hi and lo side of CELL
    Real nx = ((d+1)%2)*sgn;//if d is 0, nx is 1
    Real ny = d*sgn; //if d is 0, ny is 0

    int Q = m_Q;
    int nCells = nStenCells;

    const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit]).getFab();

    int phase = 0;
    if(psiNodesFab(iv) < 0){phase=1;}

    LAPACKMatrix MF;
    LAPACKMatrix MpEta, MvU,  STemp;

    MF.define(n, n); MF.setVal(0);

    MvU.define(nCells,n);   MvU.setVal(0);

    MpEta.define(nCells,n); MpEta.setVal(0);

    LAPACKMatrix eta(1,nCells);

    int row = 0;int comp;
    Real mom;
    IntVect position;

    int rad = (Q+1)/2;

    for(int l = -r;l<r+1;l++)
    {
        for(int c = -r;c<r+1;c++) 
        {
            position[0] = c; position[1]= l;
            //skip these rows bc they are not in the stencil for this face
            if(abs(l) > rad || abs(c) > rad){row+=1; continue;}
            if(d==0 && s==0 && c==rad){row+=1; continue;}
            if(d==0 && s==1 && c==-rad){row+=1; continue;}
            if(d==1 && s==0 && l==rad){row+=1; continue;}
            if(d==1 && s==1 && l==-rad){row+=1; continue;}

            eta(0,row) = etaFab(iv+position,phase);


            for(int px =0; px<=Q; px++)
            {
                for(int py = 0; py<=Q-px;py++)
                {
                    if(d == 0 && px == Q){continue;} //dont need this columns bc we dont have the rank for them
                    if(d == 1 && py == Q){continue;}

                    comp=  px*Q + (3*px - px*px)/2 + py;
                    mom = (pow(c+.5,px+1) - pow(c-.5,px+1) ) / (px+1);
                    mom*= (pow(l+.5,py+1) - pow(l-.5,py+1) ) / (py+1);
                    MvU(row,comp) = mom;

                    MpEta(row,comp) = pow(c,px) * pow(l,py);
                }
            }
            row+=1;
        }
    }

    //now build up the MFx etc matrices
    row = 0; //corresponds to the q coefficient of u

    for(int qx = 0; qx<=Q; qx++) //i know youre not supposed to do loops like this but whatever these numbers are small ¯\_(ツ)_/¯ 
    {
        for(int qy = 0;qy <= Q - qx; qy++)
        {
            for(int px = 0; px<=Q; px++)
            {
                for(int py =0; py<= Q -px; py++)
                {
                    comp =  (px)*Q + (3*(px) - (px)*(px)) /2 + (py);

                    MF(row,comp) = qx * faceMom(px+qx-1,py+qy,d,s) * nx + qy * faceMom(px+qx,py+qy-1,d,s) * ny;
                    
                }
            }
            row+=1;
        }
    }


    //SEta += pinv(MpEta)^T MFx^T pinv(Mv) 
    MvU.pseudoInvertUsingSVD(10,1e-10);

    LAPACKMatrix pM(n, 2*nStenCells); 
    pM.setVal(0);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<nStenCells;j++)
        {
            pM(i,j) = MvU(i,j);
        }
    }

    MpEta.pseudoInvertUsingSVD(10,1e-10);
    MpEta.transpose();
    LAPACKMatrix etaCoef, etaCoefU, etaCoefV;

    multiply(etaCoef, eta, MpEta);

    MF.transpose();

    LAPACKMatrix SFxU, SFxV, SFyU, SFyV,SF;
    SF.define(1,n); SF.setVal(0);

    multiply(SF, etaCoef, MF);

    multiply(STemp, SF, pM);
    for(int i=0;i<2*nStenCells;i++)
    {
        regFaceFlux[i] = STemp(0,i);
    }

}

//C vects have the correction terms
void 
EBEllipticSetUp::getPinvMSten(const IntVect center, const LAPACKMatrix& betaCoef,
                const LAPACKMatrix& etaCoef, const DataIndex& dit,
                LAPACKMatrix& WMpW, LAPACKMatrix& CJump)
{
    CH_TIME("EBEllipticSetUp::getPinvMSten");
    const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit]).getFab();
    const IVSFAB<Real>& momentsFab = (*m_moments)[dit];
    const IVSFAB<Real>& jumpFab = (*m_jumps)[dit];

    const IntVectSet& EBCellsFab = momentsFab.getIVS();
    IntVectSet EBCellsSten;
    EBCellsSten.define(EBCellsFab);
    Box stenBox(center - r*IntVect::Unit, center + r*IntVect::Unit);
    EBCellsSten &= stenBox;
    int nCutCells = EBCellsSten.numPts();

    Real wv,wa,d;//this is our weight for volume integrals AND derivatives - based on d^-m_weight
    int comp;
    Real mom;
    Real sx, sy;
    Vector<Real> moments(n,0);

    Vector<Real> volMomentsLo(n,0);
    Vector<Real> volMomentsHi(n,0); 
    Vector<Real> areaMoments(n,0);
    Vector<Real> nXMoments(n,0);
    Vector<Real> nYMoments(n,0);


    int nCells = nStenCells;
    int nFullCells = nCells - nCutCells;
    int nJ,col;
    int nJCells = nCutCells;
    nJ = 4; //vol, vol,  u, eta f

    int m = 2*nStenCells + nCutCells * 2; //nStenCells + nCutCells is volumes + 2*nCutCells
    
    
    LAPACKMatrix M, W, Ma, Wa, WvGrd, WvFlt;
    LAPACKMatrix Mu, Mv, Muf, Mvf;
    LAPACKMatrix PhiJump;// 
    PhiJump.define(nCutCells*2,1); PhiJump.setVal(0);
    
    Mu.define(nStenCells, 2*n); Mu.setVal(0); //all full cells and grd part of cut cells
    Muf.define(nStenCells, 2*n);  Muf.setVal(0); //floating cut cells
    
    M.define(m,2*n);
    M.setVal(0);
    W.define(m,m); W.setToIdentity(); //weighting matrix

    Ma.define(2*nCutCells, 2*n); Ma.setVal(0); //just contains jump conditions - one area and one derivative


    Wa.define(2*nCutCells,1); //area weight matrix diag
    WvGrd.define(nStenCells,1); //volume weight matrix diag
    WvFlt.define(nStenCells,1); //volume weight matrix diag

    WvGrd.setVal(1); WvFlt.setVal(1);

    //need to keep track of these seperately
    int rowA,rowV;
    rowA = rowV = 0;
    
    IntVect iv; IntVect right(1,0); IntVect up(0,1);
    bool centerPos = (psiNodesFab(center) >0); //is my center cell positive phase?
    IntVect position;
    int compx,compy;

    Vector<Real> uNRow(2*n,0);
    Vector<Real> fXRow(2*n,0);
    Vector<Real> fYRow(2*n,0);
    Vector<Real> uxRow(2*n,0);
    Vector<Real> uyRow(2*n,0);
    LAPACKMatrix LxRow, LyRow;

    int cx = 1*m_Q + (3 - 1)/2 + 0;
    int cy = 1;

    //iterate through from lower left, and across, etc
    //Mc = phi ; c=[c_u^g c_u^f c_v^g c_v^f]
    //star shaped stencil
    // Box centerBox(center - (m_Q-1)*IntVect::Unit, center + (m_Q-1)*IntVect::Unit );
    // IntVectSet nbSet(centerBox); nbSet |= center + m_Q*up; nbSet |= center - m_Q*up;
    // nbSet |= center + m_Q*right; nbSet |= center - m_Q*right; 

    int rc = 0;

    Vector<Real> gEBg(2*n,0);
    Real nx,ny;

    LAPACKMatrix MFGrd, MFFlt;

    LAPACKMatrix etaCtd, betaCtd; LAPACKMatrix ctdMoments(1,2*n);
    LAPACKMatrix ctdMomentsB(1,n);
    Real etaCtdAvg;
    Real betaCtdAvg;
    Real wf,wL;
    int rowL = 0;
    LAPACKMatrix LxRowFlt, LyRowFlt, LxRowGrd, LyRowGrd; 

    LAPACKMatrix etaLo, etaHi;
    LAPACKMatrix etaIntLo, etaIntHi;
    LAPACKMatrix momLo, momHi; momLo.define(n,1); momHi.define(n,1);
    LAPACKMatrix etaCoefGrd, etaCoefFlt;
    etaCoefGrd.define(1,n); etaCoefFlt.define(1,n); 
    for(int i=0;i<n;i++){etaCoefGrd(0,i) = etaCoef(i,0); etaCoefFlt(0,i)=etaCoef(i+n,0);}

    Real rhsScale = m_dx*m_dx;

    Real gradScale; int grdSgn, fltSgn;

    for(int l = -r;l<r+1;l++)
    {
        for(int c = -r;c<r+1;c++) //c++ hahhaa
        {
            position[0] = c; position[1]= l;
            iv = center + position;
            
            sx = position[0]; sy = position[1];

            d = sqrt(sx*sx + sy*sy) + 1; //we offset this to center node. now center cell will have weight 1
            wv= pow(d,-1*m_weight); //experiment with weighting. Hans says use Q+1
            wa = wv; //not sure how to weight the area integrals.
            wL = wv;

            if(EBCellsFab.contains(iv)) //its a cut cell
            {
                for(int i =0;i<n;i++)
                {
                    volMomentsLo[i] = momentsFab(iv,i);
                    volMomentsHi[i] = momentsFab(iv,i+hiVol);
                    areaMoments[i] = momentsFab(iv,i+area);
                    nXMoments[i] = momentsFab(iv,i+nX);
                    nYMoments[i] = momentsFab(iv,i+nY);
                }
                shiftMoments(volMomentsLo, m_Q, sx, sy); 
                shiftMoments(volMomentsHi, m_Q, sx, sy);

                
                //can optionally just add in jump conditions for neighboring cells
                if(1)//iv == center || iv == (center+right) || iv == (center-right) || iv == (center+up) || iv == (center-up) )
                {
                    shiftMoments(areaMoments, m_Q, sx,sy);
                    shiftMoments(nXMoments, m_Q, sx, sy);
                    shiftMoments(nYMoments, m_Q, sx, sy);

                    for(int i=0;i<n;i++)
                    {
                        ctdMoments(0,i) = areaMoments[i] / areaMoments[0];
                        ctdMoments(0,i+n) = areaMoments[i] / areaMoments[0];
                    }
                    multiply(etaCtd, ctdMoments, etaCoef);
                    multiply(betaCtd, ctdMoments, betaCoef);
                    etaCtdAvg = etaCtd(0,0) / 2;
                    betaCtdAvg = betaCtd(0,0);

                    
                    wf = wa / etaCtdAvg;// * (1.0 - exp( -10.0*etaCtdAvg/(betaCtdAvg * pow(m_dx,2) ) ) ); //dimensionless quantity
                    if(!factorEta){wf = wa;}

                    //get a stencil for flux so we can put in jump
                    getEBFluxMoments(nXMoments, nYMoments, etaCoef, MFGrd, MFFlt);
                    
                    for(int i =0;i<n;i++)
                    {
                        Ma(rowA,i) = wa*areaMoments[i] / areaMoments[0]; Ma(rowA,i+n) = -1.0*wa*areaMoments[i] / areaMoments[0];

                        Ma(rowA+1,i) = wf*MFGrd(0,i) / areaMoments[0]; Ma(rowA+1,i+n) = -1.0*wf*MFFlt(0,i+n) / areaMoments[0]; 
                    }
                    PhiJump(rowA,0) = jumpFab(iv,0); PhiJump(rowA+1,0) = jumpFab(iv,1);

                    Wa(rowA,0) = wa; //area
                    Wa(rowA+1,0) = wf; //area
                        
                    rowA+=2;
                }

                //fill in volume stuff
                if( psiNodesFab(iv) > 0) //lower part of this cell is in positive phase
                {

                    for(int i =0;i<n;i++)
                    {
                        Mu(rowV, i) = wv*volMomentsLo[i] / volMomentsLo[0]; //u grounded

                        Muf(rowV, i + n) = wv*volMomentsHi[i] / volMomentsHi[0]; //u floating

                    }
                }
                else
                {
                    for(int i =0;i<n;i++) //upper part of cell is is positive phase
                    {
                        Mu(rowV, i) = wv*volMomentsHi[i] / volMomentsHi[0]; //u grounded

                        Muf(rowV, i + n) = wv*volMomentsLo[i] / volMomentsLo[0]; //u floating
                    }
                }


                WvGrd(rowV,0) = wv;
                WvFlt(rowV,0) = wv;
            
                rowV+=1;
            }
            else //its a full cell
            {
                //wvGrd = wvFlt = 1;
                for(int px =0; px<=m_Q; px++)
                {
                    for(int py = 0; py<=m_Q-px;py++)
                    {
                        comp=  px*m_Q + (3*px - px*px)/2 + py;
                        mom = (pow(sx+.5,px+1) - pow(sx-.5,px+1) ) / (px+1);
                        mom*= (pow(sy+.5,py+1) - pow(sy-.5,py+1) ) / (py+1);
                        moments[comp] = mom;
                        momLo(comp,0) = mom;

                    }
                }
                if( psiNodesFab(iv) > 0 ) 
                {
                    for(int i = 0;i<n;i++)
                    {
                        Mu(rowV,i) = wv*moments[i];
                    }

                    WvGrd(rowV,0) = wv;
                

                }
                else 
                {
                    for(int i = 0;i<n;i++)
                    {
                        Mu(rowV,i + n) = wv*moments[i];
                    }

                    WvGrd(rowV,0) = wv;

                }

                rowV+=1;
            }
        }
    }

    for(int i =0;i<nStenCells;i++)
    {
        W(i,i) = WvGrd(i,0);

        W(i+nStenCells,i+nStenCells) = WvFlt(i,0);

        for(int j = 0;j<2*n;j++)
        {
            M(i,j) = Mu(i,j);
            M(i + nStenCells,j) = Muf(i,j);
        }
    }

    //printf("I took out the area rows \n");
    //fill bottom rows with area so we can throw em away later
    for(int i =2*nStenCells;i<m;i++)
    {
        W(i,i) = Wa(i-(2*nStenCells),0);
        for(int j = 0;j<2*n;j++)
        {
            M(i,j) = Ma(i-2*nStenCells,j);
        }
    }



    int rv;
    rv = pseudoInvert(M); //M = W*M already 
    LAPACKMatrix WMpW1, WMpWL, WMpWA;

    WMpWA.define(2*n, 2*nCutCells);
    WMpW.define(2*n, 2*nStenCells);

    multiply(WMpW1, M, W); //multiply by W on the right
    for(int j= 0;j<2*nStenCells;j++)
    {
        for(int i=0;i<2*n;i++)
        {
            WMpW(i,j) = WMpW1(i,j);
        }
    }

    for(int j= 2*nStenCells;j<m;j++)
    {
        for(int i=0;i<2*n;i++)
        {
            WMpWA(i,j-(2*nStenCells) ) = WMpW1(i,j);
        }
    }


    multiply(CJump, WMpWA,PhiJump);

}

void 
EBEllipticSetUp::getPinvMStenFull(const IntVect center, const LAPACKMatrix& betaCoef,
                    const LAPACKMatrix& etaCoef, const DataIndex& dit,
                    LAPACKMatrix& WMpW)

{
    CH_TIME("EBEllipticSetUp::getPinvMStenFull");
    const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit]).getFab();
    const IVSFAB<Real>& momentsFab = (*m_moments)[dit];

    const IntVectSet& EBCellsFab = momentsFab.getIVS();
    IntVectSet EBCellsSten;
    EBCellsSten.define(EBCellsFab);
    Box stenBox(center - r*IntVect::Unit, center + r*IntVect::Unit);
    EBCellsSten &= stenBox;
    int nCutCells = EBCellsSten.numPts();

    int comp;
    Real mom;
    Real sx,sy;

    Vector<Real> moments(n,0);
    Vector<Real> volMomentsLo(n,0);
    Vector<Real> volMomentsHi(n,0);

    int nCells = nStenCells;

    LAPACKMatrix Mg, Mf, M, Wg, Wf, W, ML, WL;

    Mg.define(nCells, n); Mg.setVal(0);
    Mf.define(nCells,n); Mf.setVal(0);
    Wg.define(nCells,1); Wg.setVal(0);
    Wf.define(nCells,1); Wf.setVal(0);

    IntVect iv; IntVect right(1,0); IntVect up(0,1);
    bool centerPos = (psiNodesFab(center) >0); //is my center cell positive phase?
    IntVect position;

    LAPACKMatrix LxRow, LyRow;
    // Box centerBox(center - (m_Q-1)*IntVect::Unit, center + (m_Q-1)*IntVect::Unit );
    // IntVectSet nbSet(centerBox); nbSet |= center + m_Q*up; nbSet |= center - m_Q*up;
    // nbSet |= center + m_Q*right; nbSet |= center - m_Q*right; 

    bool grd = false;

    if(psiNodesFab(center) > 0) //center cell is grd
    {
        grd = true;
    }

    int row = 0;
    Real d,wv;
    int rowL = 0;
    Real wL;
    //we need to scale the rows by the size of this coefficient for conditioning porpises
    LAPACKMatrix etaLo, etaHi;
    LAPACKMatrix etaIntLo, etaIntHi;
    LAPACKMatrix momLo, momHi; momLo.define(n,1); momHi.define(n,1);
    LAPACKMatrix etaCoefGrd, etaCoefFlt;
    etaCoefGrd.define(1,n); etaCoefFlt.define(1,n); 
    for(int i=0;i<n;i++){etaCoefGrd(0,i) = etaCoef(i,0); etaCoefFlt(0,i)=etaCoef(i,0);}
    Real rhsScale = m_dx*m_dx;

    for(int l = -r;l<r+1;l++)
    {
        for(int c = -r;c<r+1;c++) //c++ hahhaa
        {
            position[0] = c; position[1]= l;
            iv = center + position;
            
            sx = position[0]; sy = position[1];

            d = sqrt(sx*sx + sy*sy) + 1; //we offset this to center node. now center cell will have weight 1
            wv = pow(d,-1*m_weight); //experiment with weighting. Hans says use Q+1
            wL = wv;
            
            if(EBCellsSten.contains(iv)) //its a cut cell. use only the part that is the same phase
            {
                for(int i =0;i<n;i++)
                {
                    volMomentsLo[i] = momentsFab(iv,i);
                    volMomentsHi[i] = momentsFab(iv,i+hiVol);
                }
                shiftMoments(volMomentsLo, m_Q, sx, sy); 
                shiftMoments(volMomentsHi, m_Q, sx, sy);

                for(int i =0;i<n;i++)
                {
                    momLo(i,0) = volMomentsLo[i]; momHi(i,0) = volMomentsHi[i];
                }

                //fill in volume stuff
                if(grd) //center cell is grounded
                {
                    if( (psiNodesFab(iv) > 0) == grd) //lower part of cell is same phase
                    {
                        for(int i =0;i<n;i++)
                        {
                            Mg(row, i) = wv*volMomentsLo[i] / volMomentsLo[0]; //u grounded
                        }
                    }
                    else
                    {
                        for(int i =0;i<n;i++)
                        {
                            Mg(row, i) = wv*volMomentsHi[i] / volMomentsHi[0]; //u grounded
                            
                        }

                    }
                    Wg(row,0) = wv;
                    
                }
                else
                {
                    if( (psiNodesFab(iv) > 0) == grd) //lower part of cell is same phase
                    {
                        for(int i =0;i<n;i++)
                        {
                            Mf(row, i) = wv*volMomentsLo[i] / volMomentsLo[0]; //u floatin
                        }
                    }
                    else
                    {       
                        for(int i =0;i<n;i++)
                        {
                            Mf(row, i) = wv*volMomentsHi[i] / volMomentsHi[0]; //u floatin
                        }
                    }
                    Wf(row,0) = wv;
                    
                }

                row+=1;
                
            }
            else if((psiNodesFab(iv) > 0) ==  grd) //full cell of same phase
            {
                for(int px =0; px<=m_Q; px++)
                {
                    for(int py = 0; py<=m_Q-px;py++)
                    {
                        comp=  px*m_Q + (3*px - px*px)/2 + py;
                        mom = (pow(sx+.5,px+1) - pow(sx-.5,px+1) ) / (px+1);
                        mom*= (pow(sy+.5,py+1) - pow(sy-.5,py+1) ) / (py+1);
                        moments[comp] = mom;
                        momLo(comp,0) = mom;
                    }
                }

                for(int i = 0;i<n;i++)
                {
                    Mg(row,i) = wv*moments[i];
                }

                Wg(row,0) = wv;
                
            
                row+=1;
            }
            else
            {
                row +=1; //just have a row of zeros
            }
        }
    }
    
    M.define(2*nCells, n); M.setVal(0);
    W.define(2*nCells, 2*nCells); W.setVal(0);

    for(int i =0;i<nStenCells;i++)
    {
        W(i,i) = Wg(i,0);
        W(i+nStenCells,i+nStenCells) = Wf(i,0);

        for(int j = 0;j<n;j++)
        {
            M(i,j) = Mg(i,j);
            M(i + nStenCells,j) = Mf(i,j);
        }
    }
    

    pseudoInvert(M); //M = W*M already
    
    multiply(WMpW, M, W); //multiply by W on the right

}



#include "NamespaceFooter.H"