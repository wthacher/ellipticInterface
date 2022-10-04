#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BCSetUp.H"
#include "NamespaceHeader.H"

/* 
1) Given beta and eta fields: go solve for coefficients
    0a) locally solve for beta and eta coefficients
    1a) create stencils for flux across EB - use in LS system
    2a) once you have solved for pinv(WM)W, use it to build flux stencils on faces - just direclty put in BCOPStencilFAB
    3a) Create eval and gradient stencils - ACTUALLY dont think this is necessary - can just interpolate values of u from that side
     - this doesnt need to be coupled into jump conditions
2) go 'head and solve the linear systemo
3) using solution, eval new eta and beta fields at centroids - simple interpolation :)


*/

 
void
BCSetUp::getFluxAndFriction(const LevelData<FArrayBox>& cellBeta, const LevelData<FArrayBox>& cellEta,
const LevelData<IVSFAB<Real>>& phiCut, const LevelData<FArrayBox>& phi)
{
    /* CH_TIME("BCSetUp::getFluxAndFriction");
    //this will hold duStencils for ghost cells on the edges that we can exchange to complete stencils on edge
    DisjointBoxLayout dbl ( (*m_psiNodes).disjointBoxLayout() );
    
    for(DataIterator dit = (*m_psiNodes).dataIterator(); dit.ok(); ++dit)
    {
        const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit()]).getFab();
        const IVSFAB<Real>& momentsFab = (*m_moments)[dit()];

        const IVSFAB<Real>& phiCutFab = phiCut[dit()];
        const FArrayBox& phiFab = phi[dit()];

        const IntVectSet& BCCellsFab = momentsFab.getIVS();
        
        const FArrayBox& etaFab = cellEta[dit()];
        const FArrayBox& betaFab = cellBeta[dit()];

        const IntVectSet& taggedCellsFab = ((*m_taggedFluxCells)[dit()]).getIVS();

        Box valid = dbl[dit];
        valid.grow(IntVect::Unit); //so we can do calculations in the first ghost row
        // valid.enclosedCells();
        // valid.grow(-r*IntVect::Unit);

        IntVectSet validBC = BCCellsFab;
        validBC&=valid;

        IntVectSet validTagged = taggedCellsFab;
        validTagged&=valid;

        IntVectSet validReOp = validBC;

        IVSFAB<EBStencil>& BCOpStencilFab = (*BCOpStencil)[dit()];

        //if we haven't already defined on correct IVS, define it
        bool initial = false;
        if( !(BCOpStencilFab.m_ivs == validReOp))
        {
           BCOpStencilFab.define(validReOp,2); 
           //BCOpStencilCutFab.define(validBC,4);
           initial = true;
        }
        

        FArrayBox cellNodes, addStencil;
        Box stenBox, nodeStenBox, cellBox;

        IntVectSet BCCellsSten;
        LAPACKMatrix WMpW, GBeta, betaS, A;
        

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
        LAPACKMatrix phiRhs(4*nStenCells,1); phiRhs.setVal(0);
        LAPACKMatrix phiCoef(4*n,1);

        Vector<Real> fluxFaceCell, fluxEBFaceCell, fluxVect, frictionCell, frictionEBCell, frictionVect, regFaceFlux;
        fluxFaceCell.resize(6*2*4*nStenCells);
        fluxEBFaceCell.resize(2*2*4*nStenCells);
        frictionCell.resize(2*4*nStenCells);
        frictionEBCell.resize(2*4*nStenCells);

        fluxVect.resize(4*nStenCells);
        frictionVect.resize(4*nStenCells);

        regFaceFlux.resize(2*4*nStenCells);

        Real faceW = .5;
        if(!conserve)
        {
            faceW = 1;
        }

        //initialize and clear
        Box initBox;
        for(IVSIterator ivIt(validBC);ivIt.ok();++ivIt)
        {
            iv = ivIt();
            if(initial)
            {
                for(int c = 0;c<4;c++) //loop through the 2 phases, 2 flux components
                {
                    EBStencil& currEBStencil = BCOpStencilCutFab(iv,c);
                    initBox.define(iv - (m_Q+1)*IntVect::Unit, iv + (m_Q+1)*IntVect::Unit);
                    currEBStencil.defineBox(4, initBox);
                }
            }

            for(int c = 0;c<4;c++) //loop through the 2 phases, 2 flux components
            {
                EBStencil& currEBStencil = BCOpStencilCutFab(iv,c);
                currEBStencil.clear();
            }
        }

        for(IVSIterator ivIt(validTagged); ivIt.ok(); ++ivIt)
        {
            iv = ivIt();

            if(initial)
            {
                for(int c = 0;c<2;c++) //loop through the 2 flux components
                {
                    EBStencil& currEBStencil = BCOpStencilFab(iv,c);
                    initBox.define(iv - (m_Q+1)*IntVect::Unit, iv + (m_Q+1)*IntVect::Unit);
                    currEBStencil.defineBox(4, initBox);
                }
            }
            
            for(int c = 0;c<2;c++) //loop through the 2 flux components
            {
                EBStencil& currEBStencil = BCOpStencilFab(iv,c);
                currEBStencil.clear();
            }
        }

        Real volGrd;

        Real centerX, centerY;
        LAPACKMatrix centerVal, centerM;

        for(IVSIterator ivIt(validBC);ivIt.ok();++ivIt)
        {
            iv = ivIt();

            // if( iv[0] == 39 && iv[1] ==37 )
            // {
            //     printf("adsfadsfadf \n");
            // }

            //get box from which we take neighbors
            stenBox.define(iv-r*IntVect::Unit, iv+r*IntVect::Unit);

            getCoefficients(iv, etaFab, betaFab, betaCoef, etaCoef, dit());
            //WMc = Wphi. c=[c_u^g c_u^f c_v^g c_v^f] phi = [u v u_f v_f]
            
            getPinvMSten(iv, betaCoef, etaCoef, dit(), WMpW); //get WM^+W for this stencil

            counter=0;

            cellBox.define(iv,iv);
            cellBox.surroundingNodes();  cellBox.growHi(0); cellBox.growHi(1);
            cellNodes.define(cellBox,1);
            cellNodes.copy(psiNodesFab);

            getFluxFaceCell(WMpW, cellNodes, faceIntValsFab, etaCoef, fluxFaceCell); //2*6 vectors of length 4*nStenCells: L_x, L_y

            getFluxEBFaceCell(WMpW, iv, dit(), etaCoef, fluxEBFaceCell); //2*2 vectors of length 4*nStenCells

            getFrictionEBCell(WMpW, iv, dit(), betaCoef, frictionEBCell);//2 vectors of length 4*nStenCells

            if(psiNodesFab(iv) > 0){volGrd = momentsFab(iv,0);}
            else{volGrd = momentsFab(iv,hiVol);}
            
            counter = 0;
            for(int i=0;i<2;i++) //for grd and floating
            {
                for(int c=0;c<2;c++) //for each comp of flux
                {
                    EBStencil& currEBStencil = BCOpStencilCutFab(iv,2*i + c); //phase + component
                    
                    for(int j=0;j<4*nStenCells;j++)
                    {
                        fluxVect[j] = fluxEBFaceCell[counter + j];
                    }
                    counter += 4*nStenCells;

                    currEBStencil.addVect(stenBox,fluxVect,1,initial);

                    if(i == 0) //for grd boys only
                    {
                        for(int j=0;j<4*nStenCells;j++)
                        {
                            frictionVect[j] = frictionEBCell[c*4*nStenCells + j];
                        }
                        currEBStencil.addVect(stenBox, frictionVect, -1*m_dx*m_dx, initial);
                        //currEBStencil.addIV(iv, -1.0*betaFab(iv)*m_dx*m_dx,c,initial);
                    }
                    

                }
            }
            assert(counter == 2*2*4*nStenCells);

            //loop through the cell faces and build the stencil the flux stencil FOR THIS CELL AND NEIGHBOR
            
            counter = 0; //counter for duFaceVect
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

                        for(int c = 0;c<2;c++) //loop through x and y components of FLUX stencils, which both use u AND v
                        {
                            EBStencil& currEBStencil = BCOpStencilCutFab(iv,2*k + c); //grd grd flt flt

                            for(int j=0;j<4*nStenCells;j++)
                            {
                                fluxVect[j] = fluxFaceCell[counter + j];
                            }
                            counter += 4*nStenCells;

                            currEBStencil.addVect(stenBox, fluxVect, faceW, initial);

                            if(conserve)
                            {
                                if(validBC.contains(neighIv) )//neighbor is fellow cut cell
                                {
                                    EBStencil& currNeighEBStencil = BCOpStencilCutFab(neighIv,2*k + c); //grd grd flt flt
                                    currNeighEBStencil.addVect(stenBox, fluxVect, -.5, initial); //normal is going the other way

                                }
                                else if(validTagged.contains(neighIv))//neighbor is a tagged cell
                                {
                                    EBStencil& currNeighEBStencil = BCOpStencilFab(neighIv,c); //grd grd flt flt
                                    currNeighEBStencil.addVect(stenBox, fluxVect, -.5, initial); //normal is going the other way
                                }
                                else{continue;}
                            }


                            
                        }
                    }

                }
            }
            
            CH_assert(counter == (2*6*nStenCells*4));

            //now find the actual coefficients
            phiRhs.define(WMpW.dims().second,1);
            phiRhs.setVal(0);//uFull, vFull, uFlt, vFlt
            row=0;

            for(int l = -r;l<r+1;l++)
            {
                for(int c = -r;c<r+1;c++) //c++ hahhaa
                {
                    position[0] = c; position[1]= l;

                    if(BCCellsFab.contains(iv+position) ) //its a cut cell
                    {
                        phiRhs(row,0) = phiCutFab(iv+position,0);
                        phiRhs(row + nStenCells,0) = phiCutFab(iv+position,1);
                        phiRhs(row + 2*nStenCells,0) = phiCutFab(iv+position,2);
                        phiRhs(row + 3*nStenCells,0) = phiCutFab(iv+position,3);
                    }
                    else
                    {
                        phiRhs(row,0) = phiFab(iv+position,0);
                        phiRhs(row + nStenCells,0) = phiFab(iv+position,1);
                    }
                    row+=1;
                }
            }

            multiply(phiCoef, WMpW, phiRhs);
            //phiCoef.setVal(0);

            // centerX = -1*iv[0]; centerY = -1*iv[1];
            // centerM.define(1,4*n); centerM.setVal(0);
            // for(int px =0; px<=m_Q; px++)
            // {
            //     for(int py = 0; py<=m_Q-px;py++)
            //     {
            //         comp=  px*m_Q + (3*px - px*px)/2 + py;
            //         centerM(0,comp+3*n) = pow(centerX,px) * pow(centerY,py);
            //         //phiCoef(comp,0) = etaCoef(comp+n,0);
            //     }
            // }

            // multiply(centerVal, centerM, phiCoef);

            // printf("CenterVal: at iv %d %d: %1.10e \n", iv[0],iv[1],centerVal(0,0)); 
  
        }

            
        for(IVSIterator ivIt(validTagged); ivIt.ok(); ++ivIt)
        {
            iv = ivIt();

            getCoefficientsFull(iv, etaFab, betaFab, betaCoef, etaCoef, dit());

            getPinvMStenFull(iv, betaCoef, etaCoef, dit(), WMpW); //get WM^+W for this stencil

            cellBox.define(iv,iv);
            cellBox.surroundingNodes();  cellBox.growHi(0); cellBox.growHi(1);
            cellNodes.define(cellBox,1);
            cellNodes.copy(psiNodesFab);

            stenBox.define(iv-r*IntVect::Unit, iv + r*IntVect::Unit);

            //get du stencil on all faces of this cell
            getFluxFaceCellFull(WMpW, etaCoef, cellNodes, fluxFaceCell); //2*4 vectors of length 4*nStenCells
            if(psiNodesFab(iv) > 0)
            {
                getFrictionCell(WMpW, iv, dit(), betaCoef, frictionCell);
                for(int c=0;c<2;c++)
                {
                    EBStencil& currEBStencil = BCOpStencilFab(iv,c); //grd grd flt flt

                    for(int j=0;j<4*nStenCells;j++)
                    {
                        frictionVect[j] = frictionCell[c*4*nStenCells + j];
                    }

                    currEBStencil.addVect(stenBox, frictionVect, -1*m_dx*m_dx, initial);
                }
            }

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
                    if( (((*m_taggedFaces[d])[dit()] ).getIVS() ).contains(nIv)  )
                    {
                        for(int c = 0;c<2;c++) //loop through x and y components of FLUX stencils, which both use u AND v
                        {
                            EBStencil& currEBStencil = BCOpStencilFab(iv,c); //grd grd flt flt

                            for(int j=0;j<4*nStenCells;j++)
                            {
                                fluxVect[j] = fluxFaceCell[counter + j];
                                // if(isnan(fluxVect[j]))
                                // {
                                //     printf("die \n");
                                // }
                            }
                            counter += 4*nStenCells;

                            currEBStencil.addVect(stenBox, fluxVect, faceW, initial);
                            if(conserve)
                            {
                                if(validBC.contains(neighIv) )//neighbor is a cut cell
                                {
                                    if(psiNodesFab(nIv) > 0)
                                    {
                                        EBStencil& currNeighEBStencil = BCOpStencilCutFab(neighIv,c); //grd grd flt flt
                                        currNeighEBStencil.addVect(stenBox, fluxVect, -.5, initial); //normal is going the other way  
                                    }
                                    else
                                    {
                                        EBStencil& currNeighEBStencil = BCOpStencilCutFab(neighIv,2 + c); //grd grd flt flt
                                        currNeighEBStencil.addVect(stenBox, fluxVect,-.5, initial); //normal is going the other way  
                                    }
                                    
                                }
                                else if(validTagged.contains(neighIv))
                                {
                                    EBStencil& currNeighEBStencil = BCOpStencilFab(neighIv,c); //grd grd flt flt
                                    currNeighEBStencil.addVect(stenBox, fluxVect, -.5, initial); //normal is going the other way
                                }
                                else
                                {
                                    continue;
                                }
                            }
                        }
                        

                    }
                    else //not tagged face -- need to add regular face to this boyo
                    {
                        if(conserve)
                        {
                            getRegFaceStencil(s,d,etaFab,iv,dit(),regFaceFlux);
                        }
                        
                        for(int c = 0;c<2;c++) //loop through x and y components of FLUX stencils, which both use u AND v
                        {
                            if(conserve)
                            {
                                for(int j=0;j<4*nStenCells;j++)
                                {
                                    fluxVect[j] = regFaceFlux[c*4*nStenCells + j];
                                    //fluxVect[j] = fluxFaceCell[counter + j];
                                }
                            }
                            else
                            {
                                for(int j=0;j<4*nStenCells;j++)
                                {
                                    fluxVect[j] = fluxFaceCell[counter + j];
                                }
                            }
                            

                            EBStencil& currEBStencil = BCOpStencilFab(iv,c);

                            currEBStencil.addVect(stenBox, fluxVect, 1, initial);

                            counter += 4*nStenCells;


                        }


                    }
                    
                    
                    
                }
            }
            CH_assert(counter == (4*2*nStenCells*4));

        }
              
    } */


    
}

/*

void 
BCSetUp::getCoefficients(const IntVect center, const FArrayBox& etaFab, const FArrayBox& betaFab,
LAPACKMatrix& betaCoef, LAPACKMatrix& etaCoef,
                        const DataIndex& dit)
{
    //interpolate local (point) values of beta and eta on each side of the EB to get coefficients
    const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit]).getFab();
    const IVSFAB<Real>& momentsFab = (*m_moments)[dit];

    const IntVectSet& BCCellsFab = momentsFab.getIVS();
    IntVectSet BCCellsSten;
    BCCellsSten.define(BCCellsFab);
    Box stenBox(center - r*IntVect::Unit, center + r*IntVect::Unit);
    BCCellsSten &= stenBox;
    int nCutCells = BCCellsSten.numPts();

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
            if(BCCellsFab.contains(iv)) //its a cut cell
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
                }

                phiBeta(row,0) = w*betaFab(iv,0);

                phiEta(row,0) = w*etaFab(iv,0);

                row +=1;
            }
            

        }
    }

    M.pseudoInvertUsingSVD(10,1e-10);
    multiply(betaCoef, M, phiBeta);
    multiply(etaCoef, M, phiEta);
}

void 
BCSetUp::getCoefficientsFull(const IntVect center, const FArrayBox& etaFab, const FArrayBox& betaFab,
LAPACKMatrix& betaCoef, LAPACKMatrix& etaCoef,
                        const DataIndex& dit)
{
    //interpolate local (point) values of beta and eta on each side of the EB to get coefficients
    const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit]).getFab();
    const IVSFAB<Real>& momentsFab = (*m_moments)[dit];

    const IntVectSet& BCCellsFab = momentsFab.getIVS();
    IntVectSet BCCellsSten;
    BCCellsSten.define(BCCellsFab);
    Box stenBox(center - r*IntVect::Unit, center + r*IntVect::Unit);
    BCCellsSten &= stenBox;
    int nCutCells = BCCellsSten.numPts();

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
            if(BCCellsFab.contains(iv)) //its a cut cell, use the correct side
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
         

                phiBeta(row,0) = w*betaFab(iv,0);
                if(grd)
                {
                    phiEta(row,0) = w*etaFab(iv,0);
                }
                else
                {
                    phiEta(row,0) = w*etaFab(iv,1);
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
                    phiBeta(row,0) = w*betaFab(iv,0);

                    phiEta(row,0) = w*etaFab(iv,0);
                }
                row +=1;
            }
            

        }
    }

    M.pseudoInvertUsingSVD(10,1e-10);
    multiply(betaCoef, M, phiBeta);
    multiply(etaCoef, M, phiEta);
}


void 
BCSetUp::getFluxFaceCell(const LAPACKMatrix& WMpW, const FArrayBox& cellNodes, 
                         const FluxBox& faceIntValsFab, const LAPACKMatrix& etaCoef, 
                         Vector<Real>& fluxFaceCell)
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

    LAPACKMatrix MFxULo, MFxVLo, MFyULo, MFyVLo, MFxUHi, MFxVHi, MFyUHi, MFyVHi, MF;
    LAPACKMatrix MFxUGrd, MFxVGrd,  MFyUGrd,  MFyVGrd; 
    LAPACKMatrix MFxUFlt, MFxVFlt,  MFyUFlt,  MFyVFlt, SF;
    MF.define(4,4*n);

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

            MFxULo.define(n, n); MFxULo.setVal(0);
            MFxVLo.define(n, n); MFxVLo.setVal(0);
            MFyULo.define(n, n); MFyULo.setVal(0);
            MFyVLo.define(n, n); MFyVLo.setVal(0);

            MFxUHi.define(n, n); MFxUHi.setVal(0);
            MFxVHi.define(n, n); MFxVHi.setVal(0);
            MFyUHi.define(n, n); MFyUHi.setVal(0);
            MFyVHi.define(n, n); MFyVHi.setVal(0);

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
                MF.define(4,4*n); MF.setVal(0);
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

                                MFxULo(row,comp) = 4 * qx * faceMom(px+qx-1,py+qy,d,s,bds) * nx + qy * faceMom(px+qx,py+qy-1,d,s,bds) * ny;
                                MFxVLo(row,comp) = 2 * qy * faceMom(px+qx,py+qy-1,d,s,bds) * nx + qx * faceMom(px+qx-1,py+qy,d,s,bds) * ny;
                                MFyULo(row,comp) = qy * faceMom(px+qx,py+qy-1,d,s,bds) * nx     + 2 * qx * faceMom(px+qx-1,py+qy,d,s,bds) * ny;
                                MFyVLo(row,comp) = qx * faceMom(px+qx-1,py+qy,d,s,bds) * nx     + 4 * qy * faceMom(px+qx,py+qy-1,d,s,bds) * ny;

                                MFxUHi(row,comp) = 4 * qx * faceMom(px+qx-1,py+qy,d,s,bdsHi) * nx + qy * faceMom(px+qx,py+qy-1,d,s,bdsHi) * ny;
                                MFxVHi(row,comp) = 2 * qy * faceMom(px+qx,py+qy-1,d,s,bdsHi) * nx + qx * faceMom(px+qx-1,py+qy,d,s,bdsHi) * ny;
                                MFyUHi(row,comp) = qy * faceMom(px+qx,py+qy-1,d,s,bdsHi) * nx     + 2 * qx * faceMom(px+qx-1,py+qy,d,s,bdsHi) * ny;
                                MFyVHi(row,comp) = qx * faceMom(px+qx-1,py+qy,d,s,bdsHi) * nx     + 4 * qy * faceMom(px+qx,py+qy-1,d,s,bdsHi) * ny;
                                
                                
                            }
                        }
                        row+=1;
                    }
                }
                
                MFxULo.transpose();
                MFxVLo.transpose();
                MFyULo.transpose();
                MFyVLo.transpose();

                MFxUHi.transpose();
                MFxVHi.transpose();
                MFyUHi.transpose();
                MFyVHi.transpose();

                if(loGrd)
                {
                    multiply(MFxUGrd, etaCoefGrd, MFxULo);
                    multiply(MFxVGrd, etaCoefGrd, MFxVLo);
                    multiply(MFxUFlt, etaCoefFlt, MFxUHi);
                    multiply(MFxVFlt, etaCoefFlt, MFxVHi);

                    multiply(MFyUGrd, etaCoefGrd, MFyULo);
                    multiply(MFyVGrd, etaCoefGrd, MFyVLo);
                    multiply(MFyUFlt, etaCoefFlt, MFyUHi);
                    multiply(MFyVFlt, etaCoefFlt, MFyVHi);

                    for(int i = 0;i<n;i++)
                    {
                        MF(0,i) = MFxUGrd(0,i); MF(0,2*n + i) = MFxVGrd(0,i);
                        MF(1,i) = MFyUGrd(0,i); MF(1,2*n + i) = MFyVGrd(0,i);
                        MF(2,n+i) = MFxUFlt(0,i); MF(2,3*n + i) = MFxVFlt(0,i);
                        MF(3,n+i) = MFyUFlt(0,i); MF(3,3*n + i) = MFyVFlt(0,i);

                    }

                }
                else
                {
                    multiply(MFxUGrd, etaCoefGrd, MFxUHi);
                    multiply(MFxVGrd, etaCoefGrd, MFxVHi);
                    multiply(MFxUFlt, etaCoefFlt, MFxULo);
                    multiply(MFxVFlt, etaCoefFlt, MFxVLo);

                    multiply(MFyUGrd, etaCoefGrd, MFyUHi);
                    multiply(MFyVGrd, etaCoefGrd, MFyVHi);
                    multiply(MFyUFlt, etaCoefFlt, MFyULo);
                    multiply(MFyVFlt, etaCoefFlt, MFyVLo);

                    for(int i = 0;i<n;i++)
                    {
                        MF(2,i) = MFxUGrd(0,i); MF(2,2*n + i) = MFxVGrd(0,i);
                        MF(3,i) = MFyUGrd(0,i); MF(3,2*n + i) = MFyVGrd(0,i);
                        MF(0,n+i) = MFxUFlt(0,i); MF(0,3*n + i) = MFxVFlt(0,i);
                        MF(1,n+i) = MFyUFlt(0,i); MF(1,3*n + i) = MFyVFlt(0,i);

                    }
                }

                multiply(SF, MF, WMpW);

                for(int i=0;i<4*nStenCells;i++)
                {
                    fluxFaceCell[counter + i] = SF(0, i);
                    fluxFaceCell[counter + i + 4*nStenCells] = SF(1, i);
                    fluxFaceCell[counter + i + 2*4*nStenCells] = SF(2, i);
                    fluxFaceCell[counter + i + 3*4*nStenCells] = SF(3, i);

                }
                counter += 4*4*nStenCells;


            }

            else
            {
                row=0;
                MF.define(2,4*n); MF.setVal(0);
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

                                MFxULo(row,comp) = 4 * qx * faceMom(px+qx-1,py+qy,d,s,bds) * nx + qy * faceMom(px+qx,py+qy-1,d,s,bds) * ny;
                                MFxVLo(row,comp) = 2 * qy * faceMom(px+qx,py+qy-1,d,s,bds) * nx + qx * faceMom(px+qx-1,py+qy,d,s,bds) * ny;
                                MFyULo(row,comp) = qy * faceMom(px+qx,py+qy-1,d,s,bds) * nx     + 2 * qx * faceMom(px+qx-1,py+qy,d,s,bds) * ny;
                                MFyVLo(row,comp) = qx * faceMom(px+qx-1,py+qy,d,s,bds) * nx     + 4 * qy * faceMom(px+qx,py+qy-1,d,s,bds) * ny;
                            }
                        }
                        row+=1;
                    }
                }
                
                MFxULo.transpose();
                MFxVLo.transpose();
                MFyULo.transpose();
                MFyVLo.transpose();

                if(loGrd)
                {
                    multiply(MFxUGrd, etaCoefGrd, MFxULo);
                    multiply(MFxVGrd, etaCoefGrd, MFxVLo);

                    multiply(MFyUGrd, etaCoefGrd, MFyULo);
                    multiply(MFyVGrd, etaCoefGrd, MFyVLo);

                    for(int i = 0;i<n;i++)
                    {
                        MF(0,i) = MFxUGrd(0,i); MF(0,2*n + i) = MFxVGrd(0,i);
                        MF(1,i) = MFyUGrd(0,i); MF(1,2*n + i) = MFyVGrd(0,i);
                    }

                }
                else
                {
                    multiply(MFxUFlt, etaCoefFlt, MFxULo);
                    multiply(MFxVFlt, etaCoefFlt, MFxVLo);

                    multiply(MFyUFlt, etaCoefFlt, MFyULo);
                    multiply(MFyVFlt, etaCoefFlt, MFyVLo);

                    for(int i = 0;i<n;i++)
                    {
                        MF(0,n+i) = MFxUFlt(0,i); MF(0,3*n + i) = MFxVFlt(0,i);
                        MF(1,n+i) = MFyUFlt(0,i); MF(1,3*n + i) = MFyVFlt(0,i);
                    }
                }

                multiply(SF, MF, WMpW);

                for(int i=0;i<4*nStenCells;i++)
                {
                    fluxFaceCell[counter + i] = SF(0, i);
                    fluxFaceCell[counter + i + 4*nStenCells] = SF(1, i);
                }
                counter += 2*4*nStenCells;
            }
        }
    }
}

void 
BCSetUp::getFluxEBFaceCell(const LAPACKMatrix& WMpW, const IntVect iv, const DataIndex& dit,
                        const LAPACKMatrix& etaCoef, Vector<Real>& fluxEBFaceCell)
{
    //need to use EB moments, but this is similar to getFluxFaceCell
    //should probably make another function to create the stencil for EB flux as we will 
    //also use it in getPinMSten :) 

    //returns Fx grd, Fy grd, Fx flt, Fy flt
    const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit]).getFab();
    const IVSFAB<Real>& momentsFab = (*m_moments)[dit];
    bool loGrd = psiNodesFab(iv) > 0;
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

    LAPACKMatrix MFxGrd, MFyGrd, MFxFlt, MFyFlt, FxGrdSten, FyGrdSten, FxFltSten, FyFltSten;
    getEBFluxMoments(nxMoments, nyMoments, etaCoef, MFxGrd, MFyGrd, MFxFlt, MFyFlt);

    multiply(FxGrdSten, MFxGrd, WMpW);
    multiply(FyGrdSten, MFyGrd, WMpW);
    multiply(FxFltSten, MFxFlt, WMpW);
    multiply(FyFltSten, MFyFlt, WMpW);

    for(int i=0;i<4*nStenCells;i++)
    {
        fluxEBFaceCell[i] = grdSgn*FxGrdSten(0,i);
        fluxEBFaceCell[i+   4*nStenCells] = grdSgn*FyGrdSten(0,i);
        fluxEBFaceCell[i+ 2*4*nStenCells] = fltSgn*FxFltSten(0,i);
        fluxEBFaceCell[i+ 3*4*nStenCells] = fltSgn*FyFltSten(0,i);

    }

}

 void 
 BCSetUp::getEBFluxMoments(const Vector<Real>& nxMoments, const Vector<Real>& nyMoments,
                            const LAPACKMatrix& etaCoef, LAPACKMatrix& MFxGrd,
                            LAPACKMatrix& MFyGrd, LAPACKMatrix& MFxFlt, LAPACKMatrix& MFyFlt)
{
    LAPACKMatrix MFxU, MFxV, MFyU, MFyV;
    LAPACKMatrix MFxUGrd, MFxVGrd,  MFyUGrd,  MFyVGrd; 
    LAPACKMatrix MFxUFlt, MFxVFlt,  MFyUFlt,  MFyVFlt;
    
    MFxU.define(n, n); MFxU.setVal(0);
    MFxV.define(n, n); MFxV.setVal(0);
    MFyU.define(n, n); MFyU.setVal(0);
    MFyV.define(n, n); MFyV.setVal(0);

    MFxGrd.define(1,4*n); MFxGrd.setVal(0);
    MFyGrd.define(1,4*n); MFyGrd.setVal(0);
    MFxFlt.define(1,4*n); MFxFlt.setVal(0);
    MFyFlt.define(1,4*n); MFyFlt.setVal(0);


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
                    nXy = (qy>=1) ? nxMoments[compy]: 0;
                    nYx = (qx>=1) ? nyMoments[compx]: 0;
                    nYy = (qy>=1) ? nyMoments[compy]: 0;

                    MFxU(row,comp) = 4 * qx * nXx + qy * nYy;
                    MFxV(row,comp) = 2 * qy * nXy + qx * nYx;
                    MFyU(row,comp) = qy * nXy     + 2 * qx * nYx;
                    MFyV(row,comp) = qx * nXx     + 4 * qy * nYy;
                }
            }
            row+=1;
        }
    }
    
    MFxU.transpose();
    MFxV.transpose();
    MFyU.transpose();
    MFyV.transpose();

    multiply(MFxUGrd, etaCoefGrd, MFxU);
    multiply(MFxVGrd, etaCoefGrd, MFxV);
    multiply(MFxUFlt, etaCoefFlt, MFxU);
    multiply(MFxVFlt, etaCoefFlt, MFxV);

    multiply(MFyUGrd, etaCoefGrd, MFyU);
    multiply(MFyVGrd, etaCoefGrd, MFyV);
    multiply(MFyUFlt, etaCoefFlt, MFyU);
    multiply(MFyVFlt, etaCoefFlt, MFyV);

    for(int i = 0;i<n;i++)
    {
        MFxGrd(0,i) = MFxUGrd(0,i); MFxGrd(0,2*n + i) = MFxVGrd(0,i);
        MFyGrd(0,i) = MFyUGrd(0,i); MFyGrd(0,2*n + i) = MFyVGrd(0,i);
        MFxFlt(0,n+i) = MFxUFlt(0,i); MFxFlt(0,3*n + i) = MFxVFlt(0,i);
        MFyFlt(0,n+i) = MFyUFlt(0,i); MFyFlt(0,3*n + i) = MFyVFlt(0,i);

    }


}

void 
BCSetUp::getFrictionEBCell(const LAPACKMatrix& WMpW, IntVect iv, const DataIndex& dit,
                        const LAPACKMatrix& betaCoef, Vector<Real>& frictionEBCell)
{
    //this should be pretty simple, just create Mbeta
    const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit]).getFab();
    const IVSFAB<Real>& momentsFab = (*m_moments)[dit];
    bool loGrd = psiNodesFab(iv) > 0;
    LAPACKMatrix betaCoefFull(n,1);
    int Q = m_Q;
    LAPACKMatrix MBeta2;
    for(int i=0;i<n;i++)
    {
        betaCoefFull(i,0) = betaCoef(i,0);
    }

    Vector<Real> moments(n,0);
    if(loGrd)
    {
        for(int i=0;i<n;i++)
        {
           moments[i] = momentsFab(iv,i);
        }
    }
    else
    {
        for(int i=0;i<n;i++)
        {
           moments[i] = momentsFab(iv,i+ hiVol);
        }
    }

    LAPACKMatrix MBeta(n,n); MBeta.setVal(0);
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
                    MBeta(row,comp) = moments[compP]; //column corresponds to coefficient of beta, whcih has p
                }
            }
            row+=1;
        }
    }

    
    LAPACKMatrix MBetaU, MBetaV;
    MBetaU.define(1, 4*n); MBetaU.setVal(0);
    MBetaV.define(1, 4*n); MBetaV.setVal(0);

    multiply(MBeta2, MBeta, betaCoefFull );

    MBeta2.transpose();

    for(int i=0;i<n;i++)
    {
        MBetaU(0,i) = MBeta2(0,i);
        MBetaV(0,i + 2*n) = MBeta2(0,i);
    }

    LAPACKMatrix betaStenU, betaStenV;
    multiply(betaStenU, MBetaU, WMpW);
    multiply(betaStenV, MBetaV, WMpW);

    for(int i=0;i<4*nStenCells;i++)
    {
        frictionEBCell[i] = betaStenU(0,i);
        frictionEBCell[i + 4*nStenCells] = betaStenV(0,i);
    }


}


void 
BCSetUp::getPinvMSten(const IntVect center, const LAPACKMatrix& betaCoef,
                const LAPACKMatrix& etaCoef, const DataIndex& dit,
                    LAPACKMatrix& WMpW)
{
    CH_TIME("BCSetUp::getPinvMSten");
    const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit]).getFab();
    const IVSFAB<Real>& momentsFab = (*m_moments)[dit];

    const IntVectSet& BCCellsFab = momentsFab.getIVS();
    IntVectSet BCCellsSten;
    BCCellsSten.define(BCCellsFab);
    Box stenBox(center - r*IntVect::Unit, center + r*IntVect::Unit);
    BCCellsSten &= stenBox;
    int nCutCells = BCCellsSten.numPts();

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
    nJ = 8; //vol, vol, vol, vol, u,y mu fx, mu fy
    if(!useGradJump)
    {
        nJ = 6;
    }
    int m = nJCells*(nJ-4) + 4*nStenCells; //# of rows in M. for cut cells, have volume, area, normals (maybe)
    
    
    LAPACKMatrix M, W, Ma, Wa, WvGrd, WvFlt;
    LAPACKMatrix Mu, Mv, Muf, Mvf;

    Mu.define(nStenCells, 4*n); Mv.define(nStenCells, 4*n); Mu.setVal(0); Mv.setVal(0);
    Muf.define(nStenCells, 4*n); Mvf.define(nStenCells, 4*n); Muf.setVal(0); Mvf.setVal(0);
    
    M.define(m,4*n);
    M.setVal(0);
    W.define(m,m); W.setToIdentity(); //weighting matrix

    Ma.define(nJCells*(nJ-4), 4*n); Ma.setVal(0); //just contains jump conditions

    int nConstrained = 1;

    LAPACKMatrix MaConstrained(nConstrained*(nJ-4), 4*n); MaConstrained.setVal(0);
    //Mv.define(nStenCells,2*n); Mv.setVal(0); //just contains volume info
    
    Wa.define(nJCells*(nJ-4),1); //area weight matrix diag
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

    Vector<Real> uNRow(4*n,0);
    Vector<Real> fXRow(4*n,0);
    Vector<Real> fYRow(4*n,0);
    Vector<Real> uxRow(4*n,0);
    Vector<Real> uyRow(4*n,0);

    int cx = 1*m_Q + (3 - 1)/2 + 0;
    int cy = 1;

    //iterate through from lower left, and across, etc
    //Mc = phi ; c=[c_u^g c_u^f c_v^g c_v^f]
    //star shaped stencil
    // Box centerBox(center - (m_Q-1)*IntVect::Unit, center + (m_Q-1)*IntVect::Unit );
    // IntVectSet nbSet(centerBox); nbSet |= center + m_Q*up; nbSet |= center - m_Q*up;
    // nbSet |= center + m_Q*right; nbSet |= center - m_Q*right; 

    int rc = 0;

    Vector<Real> gEBg(4*n,0);
    Real nx,ny;

    LAPACKMatrix MFxGrd, MFyGrd, MFxFlt, MFyFlt;

    LAPACKMatrix etaCtd, betaCtd; LAPACKMatrix ctdMoments(1,2*n);
    LAPACKMatrix ctdMomentsB(1,n);
    Real etaCtdAvg;
    Real betaCtdAvg;
    Real wf;

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

            if(BCCellsFab.contains(iv)) //its a cut cell
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

                //fill in volume stuff
                if( psiNodesFab(iv) > 0) //lower part of this cell is in positive phase
                {
                    for(int i =0;i<n;i++)
                    {
                        Mu(rowV, i) = wv*volMomentsLo[i] / volMomentsLo[0]; //u grounded
                        Mv(rowV, i+2*n) = wv*volMomentsLo[i] / volMomentsLo[0]; //v grounded
                        Muf(rowV, i + n) = wv*volMomentsHi[i] / volMomentsHi[0]; //u floating
                        Mvf(rowV, i + 3*n) = wv*volMomentsHi[i] / volMomentsHi[0]; //v floating
                    }
                }
                else
                {
                    for(int i =0;i<n;i++) //upper part of cell is is positive phase
                    {
                        Mu(rowV, i) = wv*volMomentsHi[i] / volMomentsHi[0]; //u grounded
                        Mv(rowV, i+2*n) = wv*volMomentsHi[i] / volMomentsHi[0]; //v grounded
                        Muf(rowV, i + n) = wv*volMomentsLo[i] / volMomentsLo[0]; //u floating
                        Mvf(rowV, i + 3*n) = wv*volMomentsLo[i] / volMomentsLo[0]; //v floating
                    }
                }
                WvGrd(rowV,0) = wv;
                WvFlt(rowV,0) = wv;
                rowV+=1;
                //add in jump conditions for center cell and two neighbors?
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

                    wf = wa;// * (1.0 - exp( -10.0*etaCtdAvg/(betaCtdAvg * pow(m_dx,2) ) ) ); //dimensionless quantity
                    

                    //get a stencil for flux so we can put in jump
                    getEBFluxMoments(nXMoments, nYMoments, etaCoef, MFxGrd, MFyGrd, MFxFlt, MFyFlt);
                    //NOTE YOULL NEED TO NORMALIZE THIS ROW PROBABLY TO PREVENT WEIRD CONDITIONING


                    if(useConstrained && iv == center)// || iv == (center+right) || iv == (center-right) || iv == (center+up) || iv == (center-up) )
                    {
                        for(int i =0;i<n;i++)
                        {
                            MaConstrained(rc,i) = wa*areaMoments[i]; MaConstrained(rc,i+n) = -1.0*wa*areaMoments[i];
                            MaConstrained(rc+1,i+2*n) = wa*areaMoments[i]; MaConstrained(rc+1,i+3*n) = -1.0*wa*areaMoments[i];
                            if(useGradJump)
                            {
                                MaConstrained(rc+2,i) = wf*MFxGrd(0,i)/etaCtdAvg; MaConstrained(rc+2,i+n) = -1.0*wf*MFxFlt(0,i+n)/etaCtdAvg; 
                                MaConstrained(rc+2,i+2*n) = wf*MFxGrd(0,i+2*n)/etaCtdAvg; MaConstrained(rc+2,i+3*n) = -1.0*wf*MFxFlt(0,i+3*n)/etaCtdAvg; 

                                MaConstrained(rc+3,i) = wf*MFyGrd(0,i)/etaCtdAvg; MaConstrained(rc+3,i+n) = -1*wf*MFyFlt(0,i+n)/etaCtdAvg; 
                                MaConstrained(rc+3,i+2*n) = wf*MFyGrd(0,i+2*n)/etaCtdAvg; MaConstrained(rc+3,i+3*n) = -1.0*wf*MFyFlt(0,i+3*n)/etaCtdAvg; 

                            }
                        }
                        rc += nJ - 4;
                    }
                    else
                    {
                        for(int i =0;i<n;i++)
                        {
                            Ma(rowA,i) = wa*areaMoments[i]; Ma(rowA,i+n) = -1.0*wa*areaMoments[i];
                            Ma(rowA+1,i+2*n) = wa*areaMoments[i]; Ma(rowA+1,i+3*n) = -1.0*wa*areaMoments[i];
                            if(useGradJump)
                            {
                                Ma(rowA+2,i) = wf*MFxGrd(0,i)/etaCtdAvg; Ma(rowA+2,i+n) = -1.0*wf*MFxFlt(0,i+n)/etaCtdAvg; 
                                Ma(rowA+2,i+2*n) = wf*MFxGrd(0,i+2*n)/etaCtdAvg; Ma(rowA+2,i+3*n) = -1.0*wf*MFxFlt(0,i+3*n)/etaCtdAvg; 

                                Ma(rowA+3,i) = wf*MFyGrd(0,i)/etaCtdAvg; Ma(rowA+3,i+n) = -1*wf*MFyFlt(0,i+n)/etaCtdAvg; 
                                Ma(rowA+3,i+2*n) = wf*MFyGrd(0,i+2*n)/etaCtdAvg; Ma(rowA+3,i+3*n) = -1.0*wf*MFyFlt(0,i+3*n)/etaCtdAvg; 
                            }
                        }
                    }
                        
                    
                    Wa(rowA,0) = wa; //area
                    Wa(rowA+1,0) = wa; //area
                    if(useGradJump)
                    {
                        Wa(rowA+2,0) = wf; //flux
                        Wa(rowA+3,0) = wf;  //flux
                        // Wa(rowA+4,0) = wa; 
                        // Wa(rowA+5,0) = wa; 
                    }
                    
                    rowA+=nJ-4;
                }
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

                    }
                }
                if( psiNodesFab(iv) > 0 ) 
                {
                    for(int i = 0;i<n;i++)
                    {
                        Mu(rowV,i) = wv*moments[i];
                        Mv(rowV,i + 2*n) = wv*moments[i];
                    }
                    WvGrd(rowV,0) = wv;
                }
                else 
                {
                    for(int i = 0;i<n;i++)
                    {
                        Mu(rowV,i + n) = wv*moments[i];
                        Mv(rowV,i + 3*n) = wv*moments[i];
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
        W(i+nStenCells,i+nStenCells) = WvGrd(i,0);
        W(i+2*nStenCells,i+2*nStenCells) = WvFlt(i,0);
        W(i+3*nStenCells,i+3*nStenCells) = WvFlt(i,0);

        for(int j = 0;j<4*n;j++)
        {
            M(i,j) = Mu(i,j);
            M(i + nStenCells,j) = Mv(i,j);
            M(i + 2*nStenCells,j) = Muf(i,j);
            M(i + 3*nStenCells,j) = Mvf(i,j);
        }
    }
    //fill bottom rows with area so we can throw em away later
    for(int i =4*nStenCells;i<m;i++)
    {
        W(i,i) = Wa(i-4*nStenCells,0);
        for(int j = 0;j<4*n;j++)
        {
            M(i,j) = Ma(i-4*nStenCells,j);
        }
    }


    if(useConstrained)
    {
        int p  = nJ - 4; // # of constraints

        LAPACKMatrix K(4*n + p,4*n+p); K.setVal(0);
        LAPACKMatrix MT = M; //LS matrix
        MT.transpose();
        LAPACKMatrix MTM;
        multiply(MTM, MT, M);

        LAPACKMatrix MaConstrainedT = MaConstrained;
        MaConstrainedT.transpose();

        //fill in MTM
        for(int i =0;i<4*n;i++)
        {
            for(int j = 0;j<4*n;j++)
            {
                K(i,j) = MTM(i,j);
            }
        }
        //fill in D (constraint)
        for(int i =4*n;i<4*n+p;i++)
        {
            for(int j = 0;j<4*n;j++)
            {
                K(i,j) = MaConstrained(i - 4*n,j);
            }
        }
        //fill in DT
        for(int i =0;i<4*n;i++)
        {
            for(int j = 4*n;j<4*n+p;j++)
            {
                K(i,j) = MaConstrainedT(i,j - 4*n);
            }
        }

        Real cNum = getInverseOfConditionNumber(K);
        if(1.0/cNum > 1.0e5)
        {
            printf("Condition Number of LS Matrix at iv: %d %d %1.10e \n", center[0],center[1],1.0/cNum);
        }
        
        K.invert();
        

        LAPACKMatrix KSquare(4*n, 4*n); //we only need this much bc everything else is ZEROs
        for(int i = 0;i<4*n;i++)
        {
            for(int j = 0;j<4*n;j++)
            {
                KSquare(i,j) = K(i,j);
            }
        }

        LAPACKMatrix KinvMTW1,KinvMTW;
        multiply(KinvMTW1, KSquare, MT);
        multiply(KinvMTW, KinvMTW1, W);

        WMpW = KinvMTW;
    }
    else
    {
        // Real cNum = getInverseOfConditionNumber(M);
        // if(1.0/cNum > 1.0e5)
        // {
        //     printf("Condition Number of LS Matrix: %1.10e \n", 1.0/cNum);
        // }
        // if(center[0] == -1 && center[1] == -4)
        // {
        //     printf("asdf");
        // }
        M.pseudoInvertUsingSVD(10,1e-10); //M = W*M already
        multiply(WMpW, M, W); //multiply by W on the right
    }


}

void 
BCSetUp::getPinvMStenFull(const IntVect center, const LAPACKMatrix& betaCoef,
                    const LAPACKMatrix& etaCoef, const DataIndex& dit,
                    LAPACKMatrix& WMpW)

{
    CH_TIME("BCSetUp::getPinvMStenFull");
    const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit]).getFab();
    const IVSFAB<Real>& momentsFab = (*m_moments)[dit];

    const IntVectSet& BCCellsFab = momentsFab.getIVS();
    IntVectSet BCCellsSten;
    BCCellsSten.define(BCCellsFab);
    Box stenBox(center - r*IntVect::Unit, center + r*IntVect::Unit);
    BCCellsSten &= stenBox;
    int nCutCells = BCCellsSten.numPts();

    int comp;
    Real mom;
    Real sx,sy;

    Vector<Real> moments(n,0);
    Vector<Real> volMomentsLo(n,0);
    Vector<Real> volMomentsHi(n,0);

    int nCells = nStenCells;

    LAPACKMatrix Mg, Mf, M, Wg, Wf, W;

    Mg.define(nCells, n); Mg.setVal(0);
    Mf.define(nCells,n); Mf.setVal(0);
    Wg.define(nCells,nCells); Wg.setVal(0);
    Wf.define(nCells,nCells); Wf.setVal(0);

    IntVect iv; IntVect right(1,0); IntVect up(0,1);
    bool centerPos = (psiNodesFab(center) >0); //is my center cell positive phase?
    IntVect position;

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
    for(int l = -r;l<r+1;l++)
    {
        for(int c = -r;c<r+1;c++) //c++ hahhaa
        {
            position[0] = c; position[1]= l;
            iv = center + position;
            
            sx = position[0]; sy = position[1];

            d = sqrt(sx*sx + sy*sy) + 1; //we offset this to center node. now center cell will have weight 1
            wv = pow(d,-1*m_weight); //experiment with weighting. Hans says use Q+1
            
            if(BCCellsSten.contains(iv)) //its a cut cell. use only the part that is the same phase
            {
                for(int i =0;i<n;i++)
                {
                    volMomentsLo[i] = momentsFab(iv,i);
                    volMomentsHi[i] = momentsFab(iv,i+hiVol);
                }
                shiftMoments(volMomentsLo, m_Q, sx, sy); 
                shiftMoments(volMomentsHi, m_Q, sx, sy);

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

    M.define(4*nCells, 2*n); M.setVal(0);
    W.define(4*nCells, 4*nCells); W.setVal(0);

    for(int i =0;i<nStenCells;i++)
    {
        W(i,i) = Wg(i,0);
        W(i+nStenCells,i+nStenCells) = Wg(i,0);
        W(i+2*nStenCells,i+2*nStenCells) = Wf(i,0);
        W(i+3*nStenCells,i+3*nStenCells) = Wf(i,0);

        for(int j = 0;j<n;j++)
        {
            M(i,j) = Mg(i,j);
            M(i + nStenCells,j + n) = Mg(i,j);
            M(i + 2*nStenCells,j) = Mf(i,j);
            M(i + 3*nStenCells,j + n) = Mf(i,j);
        }
    }
    // if(center[0] == -3 && center[1] == 0)
    // {
    //     printf("asdf");
    // }

    M.pseudoInvertUsingSVD(10,1e-10); //M = W*M already
    
    multiply(WMpW, M, W); //multiply by W on the right
}

//you have essenitally already written this buddy
void 
BCSetUp::getFluxFaceCellFull(const LAPACKMatrix& WMpW, const LAPACKMatrix& etaCoef, 
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

    LAPACKMatrix MFxU, MFxV, MFyU, MFyV, MF;
    LAPACKMatrix MFxUGrd, MFxVGrd,  MFyUGrd,  MFyVGrd; 
    LAPACKMatrix MFxUFlt, MFxVFlt,  MFyUFlt,  MFyVFlt, SF;
    MF.define(4,4*n);

    LAPACKMatrix etaCoefGrd(n,1);
    LAPACKMatrix etaCoefFlt(n,1);  
    for(int i=0;i<n;i++)
    {
        etaCoefGrd(i,0) = etaCoef(i,0);
        etaCoefFlt(i,0) = etaCoef(i,0);
    }
    etaCoefGrd.transpose(); etaCoefFlt.transpose();

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

            MFxU.define(n, n); MFxU.setVal(0);
            MFxV.define(n, n); MFxV.setVal(0);
            MFyU.define(n, n); MFyU.setVal(0);
            MFyV.define(n, n); MFyV.setVal(0);

            iv = ll + s*normShift;

            loGrd = false;
            if(cellNodes(iv) > 0){loGrd = true;}

            row=0;
            MF.define(2,2*n); MF.setVal(0);
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

                            MFxU(row,comp) = 4 * qx * faceMom(px+qx-1,py+qy,d,s,bds) * nx + qy * faceMom(px+qx,py+qy-1,d,s,bds) * ny;
                            MFxV(row,comp) = 2 * qy * faceMom(px+qx,py+qy-1,d,s,bds) * nx + qx * faceMom(px+qx-1,py+qy,d,s,bds) * ny;
                            MFyU(row,comp) = qy * faceMom(px+qx,py+qy-1,d,s,bds) * nx     + 2 * qx * faceMom(px+qx-1,py+qy,d,s,bds) * ny;
                            MFyV(row,comp) = qx * faceMom(px+qx-1,py+qy,d,s,bds) * nx     + 4 * qy * faceMom(px+qx,py+qy-1,d,s,bds) * ny;
                        }
                    }
                    row+=1;
                }
            }
            
            MFxU.transpose();
            MFxV.transpose();
            MFyU.transpose();
            MFyV.transpose();

            if(loGrd)
            {
                multiply(MFxUGrd, etaCoefGrd, MFxU);
                multiply(MFxVGrd, etaCoefGrd, MFxV);

                multiply(MFyUGrd, etaCoefGrd, MFyU);
                multiply(MFyVGrd, etaCoefGrd, MFyV);

                for(int i = 0;i<n;i++)
                {
                    MF(0,i) = MFxUGrd(0,i); MF(0,n + i) = MFxVGrd(0,i);
                    MF(1,i) = MFyUGrd(0,i); MF(1,n + i) = MFyVGrd(0,i);
                }

            }
            else
            {
                multiply(MFxUFlt, etaCoefFlt, MFxU);
                multiply(MFxVFlt, etaCoefFlt, MFxV);

                multiply(MFyUFlt, etaCoefFlt, MFyU);
                multiply(MFyVFlt, etaCoefFlt, MFyV);

                for(int i = 0;i<n;i++)
                {
                    MF(0,i) = MFxUFlt(0,i); MF(0,n+i) = MFxVFlt(0,i);
                    MF(1,i) = MFyUFlt(0,i); MF(1,n+i) = MFyVFlt(0,i);

                }
            }

            multiply(SF, MF, WMpW);

            for(int i=0;i<4*nStenCells;i++)
            {
                fluxFaceCell[counter + i] = SF(0, i);
                fluxFaceCell[counter + i + 4*nStenCells] = SF(1, i);
            }
            counter += 2*4*nStenCells;
            
        }
    }
}

void 
BCSetUp::getFrictionCell(const LAPACKMatrix& WMpW, IntVect iv, const DataIndex& dit,
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
    MBetaU.define(1, 2*n); MBetaU.setVal(0);
    MBetaV.define(1, 2*n); MBetaV.setVal(0);

    multiply(MBeta2, MBeta, betaCoef );

    MBeta2.transpose();

    for(int i=0;i<n;i++)
    {
        MBetaU(0,i) = MBeta2(0,i);
        MBetaV(0,i + n) = MBeta2(0,i);
    }

    LAPACKMatrix betaStenU, betaStenV;
    multiply(betaStenU, MBetaU, WMpW);
    multiply(betaStenV, MBetaV, WMpW);

    for(int i=0;i<4*nStenCells;i++)
    {
        frictionCell[i] = betaStenU(0,i);
        frictionCell[i + 4*nStenCells] = betaStenV(0,i);
    }
}

void 
BCSetUp::getRegFaceStencil(int s, int d, const FArrayBox& etaFab, 
                            const IntVect iv, const DataIndex& dit,Vector<Real>& regFaceFlux )
{
    regFaceFlux.assign(0);
    Real sgn = s==0 ? -1:1; //for hi and lo side of CELL
    Real nx = ((d+1)%2)*sgn;//if d is 0, nx is 1
    Real ny = d*sgn; //if d is 0, ny is 0

    int Q = m_Q;
    int nCells = nStenCells;

    LAPACKMatrix MFxU, MFxV, MFyU, MFyV;
    LAPACKMatrix MpEta, MvU, MvV, STemp;

    MFxU.define(n, n); MFxU.setVal(0);
    MFxV.define(n, n); MFxV.setVal(0);
    MFyU.define(n, n); MFyU.setVal(0);
    MFyV.define(n, n); MFyV.setVal(0);
    MvU.define(nCells,n);   MvU.setVal(0);
    MvV.define(nCells,n);   MvV.setVal(0);
    MpEta.define(nCells,n); MpEta.setVal(0);

    LAPACKMatrix MFx, MFy;

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

            eta(0,row) = etaFab(iv+position);


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
                    MvV(row,comp) = mom;

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

                    MFxU(row,comp) = 4 * qx * faceMom(px+qx-1,py+qy,d,s) * nx + qy * faceMom(px+qx,py+qy-1,d,s) * ny;
                    MFxV(row,comp) = 2 * qy * faceMom(px+qx,py+qy-1,d,s) * nx + qx * faceMom(px+qx-1,py+qy,d,s) * ny;
                    MFyU(row,comp) = qy * faceMom(px+qx,py+qy-1,d,s) * nx     + 2 * qx * faceMom(px+qx-1,py+qy,d,s) * ny;
                    MFyV(row,comp) = qx * faceMom(px+qx-1,py+qy,d,s) * nx     + 4 * qy * faceMom(px+qx,py+qy-1,d,s) * ny;
                    
                    
                }
            }
            row+=1;
        }
    }


    //SEta += pinv(MpEta)^T MFx^T pinv(Mv)
    MvU.pseudoInvertUsingSVD(10,1e-10);
    MvV.pseudoInvertUsingSVD(10,1e-10);
    LAPACKMatrix pM(2*n, 4*nStenCells); 
    pM.setVal(0);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<nStenCells;j++)
        {
            pM(i,j) = MvU(i,j);
            pM(i+n, j + nStenCells) = MvV(i,j);
        }
    }

    MpEta.pseudoInvertUsingSVD(10,1e-10);
    MpEta.transpose();
    LAPACKMatrix etaCoef, etaCoefU, etaCoefV;

    multiply(etaCoef, eta, MpEta);

    MFxU.transpose();
    MFxV.transpose();
    MFyU.transpose();
    MFyV.transpose();

    LAPACKMatrix SFxU, SFxV, SFyU, SFyV,SF;
    SF.define(2,2*n); SF.setVal(0);

    multiply(SFxU, etaCoef, MFxU);
    multiply(SFxV, etaCoef, MFxV);
    multiply(SFyU, etaCoef, MFyU);
    multiply(SFyV, etaCoef, MFyV);

    for(int i=0;i<n;i++)
    {
        SF(0,i) = SFxU(0,i);
        SF(0,i+n) = SFxV(0,i);
        SF(1,i) = SFyU(0,i);
        SF(1,i+n) = SFyV(0,i);
    }



    multiply(STemp, SF, pM);
    for(int i=0;i<4*nStenCells;i++)
    {
        regFaceFlux[i] = STemp(0,i);
        regFaceFlux[i + 4*nStenCells] = STemp(1,i);
    }

}

 */

//just set beta to whatever Ccoef is there
void 
BCSetUp::updateBetaAndEta(LevelData<FArrayBox>& beta, 
                    LevelData<FArrayBox>& eta,
                    const Real AFace,
                    const LevelData<FArrayBox>& Ccoef,
                    const LevelData<FArrayBox>& hCoef,
                    const LevelData<IVSFAB<Real>>& phiCut,
                    const LevelData<FArrayBox>& phi,
                    const LevelData<FArrayBox>& a_H,
                    bool constFriction, bool constMu, Real muCGrd, Real muCFlt, Real betaScale)
{
    for(DataIterator dit = beta.dataIterator();dit.ok();++dit)
    {
        const FArrayBox& phiFab = phi[dit()];
        FArrayBox& betaFab = beta[dit()];
        FArrayBox& etaFab = eta[dit()];
        const FArrayBox& HFab = a_H[dit()];
        const FArrayBox& CFab = Ccoef[dit()];
        const FArrayBox& hCoefFab = hCoef[dit()]; //grd then floating

        Real C;

        const IVSFAB<Real>& phiCutFab = phiCut[dit()];

        const FArrayBox& psiNodesFab = ((*m_psiNodes)[dit()]).getFab();
        const IVSFAB<Real>& momentsFab = (*m_moments)[dit()];

        const IntVectSet& BCCellsFab = momentsFab.getIVS();

        const IntVectSet& taggedCellsFab = ((*m_taggedFluxCells)[dit()]).getIVS();

        Box valid = (beta.disjointBoxLayout() ) [dit()];
        // valid.enclosedCells();
        // valid.grow(-r*IntVect::Unit);

        IntVectSet validBC = BCCellsFab;
        validBC&=valid;

        IntVectSet validTagged = taggedCellsFab;
        validTagged&=valid;

        IntVect iv, position, sten;
        int comp;
        Real mom;
        Real sx, sy;
        Vector<Real> momentsGrd(n,0);
        Vector<Real> momentsFlt(n,0);
        Real m = 1.0 / 3.0;
        int cx = 1*m_Q + (3 - 1)/2 + 0;
        int cy = 1;
        Real eps = 1e-12; Real delta = 0;


        LAPACKMatrix MvGrd(nStenCells,2*n); MvGrd.setVal(0);
        LAPACKMatrix uPhiGrd(nStenCells,1); uPhiGrd.setVal(0);
        LAPACKMatrix vPhiGrd(nStenCells,1); vPhiGrd.setVal(0);

        LAPACKMatrix MvFlt(nStenCells,n); MvFlt.setVal(0);
        LAPACKMatrix uPhiFlt(nStenCells,1); uPhiFlt.setVal(0);
        LAPACKMatrix vPhiFlt(nStenCells,1); vPhiFlt.setVal(0);

        LAPACKMatrix uGrd, vGrd,uxGrd,uyGrd,vxGrd,vyGrd, uFlt, vFlt,uxFlt,uyFlt,vxFlt,vyFlt;
        Real ctdXGrd, ctdYGrd, ctdXFlt, ctdYFlt;
        Real SRIGrd, SRIFlt, mu0, muFlt, muGrd;


        LAPACKMatrix SGrd, SxGrd, SyGrd, SFlt, SxFlt, SyFlt;
        LAPACKMatrix uCoefGrd, vCoefGrd, uCoefFlt, vCoefFlt;
        LAPACKMatrix hCoefGrd, hCoefFlt, hGrd, hFlt;

        LAPACKMatrix mA, mdA, uFaceG, uFaceF, uFaceGx, uFaceFx;

        int compx;

        for(IVSIterator ivIt(validBC); ivIt.ok(); ++ivIt)
        {
            iv = ivIt();
            //Make interpolation matrix and vector of values to interpolate with 
            int row = 0;
            MvGrd.define(nStenCells,n); MvGrd.setVal(0);
            uPhiGrd.define(nStenCells,1); uPhiGrd.setVal(0);
            vPhiGrd.define(nStenCells,1); vPhiGrd.setVal(0);

            MvFlt.define(nStenCells,n); MvFlt.setVal(0);
            uPhiFlt.define(nStenCells,1); uPhiFlt.setVal(0);
            vPhiFlt.define(nStenCells,1); vPhiFlt.setVal(0);

            hCoefGrd.define(n,1); 
            hCoefFlt.define(n,1);
            for(int i =0;i<n;i++)
            {
                hCoefGrd(i,0) = hCoefFab(iv,i);
                hCoefFlt(i,0) = hCoefFab(iv, i + n);
            }

            for(int l = -r;l<r+1;l++)
            {
                for(int c = -r;c<r+1;c++) //c++ hahhaa
                {
                    position[0] = c; position[1]= l;
                    sten = iv + position;

                    sx = position[0];
                    sy = position[1];

                    if(BCCellsFab.contains(sten))
                    {
                        if(psiNodesFab(sten) > 0)
                        {
                            for(int i=0;i<n;i++)
                            {
                                momentsGrd[i] = momentsFab(sten,i);
                                momentsFlt[i] = momentsFab(sten,i+hiVol);
                            }
                        }
                        else //hi side grd lo side flt
                        {
                            for(int i=0;i<n;i++)
                            {
                                momentsGrd[i] = momentsFab(sten,i + hiVol);
                                momentsFlt[i] = momentsFab(sten,i);
                            }
                        }

                        shiftMoments(momentsGrd,m_Q, sx, sy);
                        shiftMoments(momentsFlt,m_Q, sx, sy);

                        for(int i=0;i<n;i++)
                        {
                            MvGrd(row,i) = momentsGrd[i] / momentsGrd[0];
                            MvFlt(row,i) = momentsFlt[i] / momentsFlt[0];
                        }
                        uPhiGrd(row,0) = phiCutFab(sten,0);
                        vPhiGrd(row,0) = phiCutFab(sten,1);

                        uPhiFlt(row,0) = phiCutFab(sten,2);
                        vPhiFlt(row,0) = phiCutFab(sten,3);

                        row+=1;
                    }
                    else //only use grd full cells
                    {
                        sx = position[0];
                        sy = position[1];
                        if(psiNodesFab(sten) > 0)
                        {
                            for(int px =0; px<=m_Q; px++)
                            {
                                for(int py = 0; py<=m_Q-px;py++)
                                {
                                    comp=  px*m_Q + (3*px - px*px)/2 + py;
                                    mom = (pow(sx+.5,px+1) - pow(sx-.5,px+1) ) / (px+1);
                                    mom*= (pow(sy+.5,py+1) - pow(sy-.5,py+1) ) / (py+1);
                                    MvGrd(row,comp) = mom;

                                }
                            }
                            uPhiGrd(row,0) = phiFab(sten,0);
                            vPhiGrd(row,0) = phiFab(sten,1);
                        }
                        else
                        {
                            for(int px =0; px<=m_Q; px++)
                            {
                                for(int py = 0; py<=m_Q-px;py++)
                                {
                                    comp=  px*m_Q + (3*px - px*px)/2 + py;
                                    mom = (pow(sx+.5,px+1) - pow(sx-.5,px+1) ) / (px+1);
                                    mom*= (pow(sy+.5,py+1) - pow(sy-.5,py+1) ) / (py+1);
                                    MvFlt(row,comp) = mom;

                                }
                            }
                            uPhiFlt(row,0) = phiFab(sten,0);
                            vPhiFlt(row,0) = phiFab(sten,1);
                        }
                        row+=1;
                    }

                }
            }
        
            MvGrd.pseudoInvertUsingSVD(10,1e-10);
            MvFlt.pseudoInvertUsingSVD(10,1e-10);
            multiply(uCoefGrd, MvGrd, uPhiGrd);
            multiply(vCoefGrd, MvGrd, vPhiGrd);
            multiply(uCoefFlt, MvFlt, uPhiFlt);
            multiply(vCoefFlt, MvFlt, vPhiFlt);

            //use these to get u,v,ux,uy,vx,vy at cell centroid
            if(psiNodesFab(iv) > 0)
            {
                ctdXGrd = momentsFab(iv,cx) / momentsFab(iv,0);
                ctdYGrd = momentsFab(iv,cy) / momentsFab(iv,0);

                ctdXFlt = momentsFab(iv,cx + hiVol) / momentsFab(iv,hiVol);
                ctdYFlt = momentsFab(iv,cy + hiVol) / momentsFab(iv,hiVol);
            }
            else
            {
                ctdXFlt = momentsFab(iv,cx) / momentsFab(iv,0);
                ctdYFlt = momentsFab(iv,cy) / momentsFab(iv,0);

                ctdXGrd = momentsFab(iv,cx + hiVol) / momentsFab(iv,hiVol);
                ctdYGrd = momentsFab(iv,cy + hiVol) / momentsFab(iv,hiVol);
            }
            

            SGrd.define(1,n); SGrd.setVal(0);
            SFlt.define(1,n); SFlt.setVal(0);

            SxGrd.define(1,n); SxGrd.setVal(0);
            SyGrd.define(1,n); SyGrd.setVal(0);

            SxFlt.define(1,n); SxFlt.setVal(0);
            SyFlt.define(1,n); SyFlt.setVal(0);

            if(ctdXFlt==0){ctdXFlt=1e-300;} if(ctdXGrd==0){ctdXGrd=1e-300;}
            if(ctdYFlt==0){ctdYFlt=1e-300;} if(ctdYGrd==0){ctdYGrd=1e-300;}

            mA.define(1,n); mdA.define(1,n);
            for(int i =0;i<n;i++)
            {
                mA(0,i) = momentsFab(iv, area + i);
            }
            // comp=0;
            // for(int px =0; px<=m_Q; px++)
            // {
            //     for(int py = 0; py<=m_Q-px;py++)
            //     {
            //         compx=  (px-1)*m_Q + (3*(px-1) - (px-1)*(px-1))/2 + py;
            //         mdA(0,comp) =  px >=1 ? px * momentsFab(iv, area + compx ) : 0;

            //         comp+=1;
            //     }
            // }

            mdA(0,3) = 4; mdA(0,5) = -0.536059288191729;

            for(int px =0; px<=m_Q; px++)
            {
                for(int py = 0; py<=m_Q-px;py++)
                {
                    comp=  px*m_Q + (3*px - px*px)/2 + py;
                    SGrd(0,comp) = pow(ctdXGrd, px) * pow(ctdYGrd, py);
                    SFlt(0,comp) = pow(ctdXFlt, px) * pow(ctdYFlt, py);
                    SxGrd(0,comp) = px*pow(ctdXGrd, px-1) * pow(ctdYGrd, py);
                    SxFlt(0,comp) = px*pow(ctdXFlt, px-1) * pow(ctdYFlt, py);
                    SyGrd(0,comp) = py*pow(ctdXGrd, px) * pow(ctdYGrd, py-1);
                    SyFlt(0,comp) = py*pow(ctdXFlt, px) * pow(ctdYFlt, py-1);
                }
            }

            multiply(uGrd, SGrd, uCoefGrd); multiply(vGrd, SGrd, vCoefGrd);
            multiply(uFlt, SFlt, uCoefFlt); multiply(vFlt, SFlt, vCoefFlt);

            multiply(uxGrd, SxGrd, uCoefGrd); multiply(vxGrd, SxGrd, vCoefGrd); uxGrd*=(1.0/m_dx); vxGrd*=(1.0/m_dx);
            multiply(uxFlt, SxFlt, uCoefFlt); multiply(vxFlt, SxFlt, vCoefFlt); uxFlt*=(1.0/m_dx); vxFlt*=(1.0/m_dx);

            multiply(uyGrd, SyGrd, uCoefGrd); multiply(vyGrd, SyGrd, vCoefGrd); uyGrd*=(1.0/m_dx); vyGrd*=(1.0/m_dx);
            multiply(uyFlt, SyFlt, uCoefFlt); multiply(vyFlt, SyFlt, vCoefFlt); uyFlt*=(1.0/m_dx); vyFlt*=(1.0/m_dx);

            multiply(hGrd, SGrd, hCoefGrd);
            multiply(hFlt, SFlt, hCoefFlt);

            multiply(uFaceG, mA, uCoefGrd);
            multiply(uFaceF, mA, uCoefFlt);
            multiply(uFaceGx, mdA, uCoefGrd);
            multiply(uFaceFx, mdA, uCoefFlt);

            printf("jump condition dif: %1.6e, dx jump: %1.6e \n", uFaceG(0,0) - uFaceF(0,0), uFaceGx(0,0) - uFaceFx(0,0) );


            if(constFriction)
            {
                betaFab(iv) = CFab(iv) * betaScale;
            }
            else
            {
                betaFab(iv) = CFab(iv) * pow( abs( pow(uGrd(0,0),2) + pow(vGrd(0,0),2) ), (m-1)/2.0);
            }

            SRIGrd = uxGrd(0,0)*uxGrd(0,0) + vyGrd(0,0)*vyGrd(0,0) + pow(uxGrd(0,0)+vyGrd(0,0),2) + .5*pow(uyGrd(0,0)+vxGrd(0,0),2);
            SRIFlt = uxFlt(0,0)*uxFlt(0,0) + vyFlt(0,0)*vyFlt(0,0) + pow(uxFlt(0,0)+vyFlt(0,0),2) + .5*pow(uyFlt(0,0)+vxFlt(0,0),2);

            mu0 = .5*pow(AFace, -1.0/3.0);

            muGrd =  mu0*( pow( (eps + SRIGrd), -1.0/3.0) + delta ); 
            muFlt =  mu0*( pow( (eps + SRIFlt), -1.0/3.0) + delta ); 

            // if(iv[0]==20 && iv[1]==18)
            // {
            //     printf("qwer");
            // }

            if(constMu)
            {
                etaFab(iv,0) = muCGrd * hGrd(0,0);
                etaFab(iv,1) = muCFlt * hFlt(0,0);
            }
            else
            {
                etaFab(iv,0) = muGrd * hGrd(0,0);
                etaFab(iv,1) = muFlt * hFlt(0,0);
            }

            

        
        }
        
        
        bool grd;
        
        LAPACKMatrix Mv(nStenCells,n); Mv.setVal(0);
        LAPACKMatrix uPhi(nStenCells,1); uPhi.setVal(0);
        LAPACKMatrix vPhi(nStenCells,1); vPhi.setVal(0);

        LAPACKMatrix u, v,ux,uy,vx,vy;
        LAPACKMatrix uCoef, vCoef;
        Real ctdX, ctdY;
        Real SRI, mu;


        LAPACKMatrix S, Sx, Sy;
        LAPACKMatrix hCoef, h;

        Vector<Real> moments(n,0);

        for(IVSIterator ivIt(validTagged); ivIt.ok(); ++ivIt)
        {
            iv = ivIt();
            grd = (psiNodesFab(iv) > 0);
            //only use full cells from this side of the BC
            int row = 0;
            Mv.define(nStenCells,n); Mv.setVal(0);
            uPhi.define(nStenCells,1); uPhi.setVal(0);
            vPhi.define(nStenCells,1); vPhi.setVal(0);

            hCoef.define(n,1); 
            for(int i =0;i<n;i++)
            {
                hCoef(i,0) = hCoefFab(iv,i);
            }

            for(int l = -r;l<r+1;l++)
            {
                for(int c = -r;c<r+1;c++) //c++ hahhaa
                {
                    position[0] = c; position[1]= l;
                    sten = iv + position;

                    if(BCCellsFab.contains(sten))
                    {
                        if( (psiNodesFab(sten) > 0) == grd)
                        {
                            for(int i=0;i<n;i++)
                            {
                                moments[i] = momentsFab(sten,i);
                            }
                        }
                        else 
                        {
                            for(int i=0;i<n;i++)
                            {
                                moments[i] = momentsFab(sten,i + hiVol);
                            }
                        }

                        shiftMoments(moments,m_Q, position[0], position[1]);

                        for(int i=0;i<n;i++)
                        {
                            Mv(row,i) = moments[i] / moments[0];
                        }
                        if(grd)
                        {
                            uPhi(row,0) = phiCutFab(sten,0);
                            vPhi(row,0) = phiCutFab(sten,1);
                        }
                        else
                        {
                            uPhi(row,0) = phiCutFab(sten,2);
                            vPhi(row,0) = phiCutFab(sten,3);
                        }

                        row+=1;
                    }
                    else //only use full cells in same phase
                    {
                        sx = position[0];
                        sy = position[1];
                        if( (psiNodesFab(sten) > 0) == grd)
                        {
                            for(int px =0; px<=m_Q; px++)
                            {
                                for(int py = 0; py<=m_Q-px;py++)
                                {
                                    comp=  px*m_Q + (3*px - px*px)/2 + py;
                                    mom = (pow(sx+.5,px+1) - pow(sx-.5,px+1) ) / (px+1);
                                    mom*= (pow(sy+.5,py+1) - pow(sy-.5,py+1) ) / (py+1);
                                    Mv(row,comp) = mom;

                                }
                            }
                            uPhi(row,0) = phiFab(sten,0);
                            vPhi(row,0) = phiFab(sten,1);
                        }
                        
                        row+=1;
                    }

                }
            }
        
            Mv.pseudoInvertUsingSVD(10,1e-10);
            multiply(uCoef, Mv, uPhi);
            multiply(vCoef, Mv, vPhi);
            ctdX = ctdY = 0;
            //use these to get u,v,ux,uy,vx,vy at cell centroid
            S.define(1,n); S.setVal(0);

            Sx.define(1,n); Sx.setVal(0);
            Sy.define(1,n); Sy.setVal(0);

            for(int px =0; px<=m_Q; px++)
            {
                for(int py = 0; py<=m_Q-px;py++)
                {
                    comp=  px*m_Q + (3*px - px*px)/2 + py;
                    S(0,comp) = pow(ctdX, px) * pow(ctdY, py);
                    Sx(0,comp) = (px>=1) ? px*pow(ctdX, px-1) * pow(ctdY, py) : 0;
                    Sy(0,comp) = (py>=1) ? py*pow(ctdX, px) * pow(ctdY, py-1) : 0;
                }
            }

            multiply(u, S, uCoef); multiply(v, S, vCoef);

            multiply(ux, Sx, uCoef); multiply(vx, Sx, vCoef); ux*=(1.0/m_dx); vx*=(1.0/m_dx);

            multiply(uy, Sy, uCoef); multiply(vy, Sy, vCoef); uy*=(1.0/m_dx); vy*=(1.0/m_dx);

            multiply(h, S, hCoef);
            if(psiNodesFab(iv) > 0)
            {
                if(constFriction)
                {
                    betaFab(iv) = CFab(iv) * betaScale;
                }
                else
                {
                    betaFab(iv) = CFab(iv) * pow( abs( pow(uGrd(0,0),2) + pow(vGrd(0,0),2) ), (m-1)/2.0);
                }
                
            }

            SRI = ux(0,0)*ux(0,0) + vy(0,0)*vy(0,0) + pow(ux(0,0)+vy(0,0),2) + .5*pow(uy(0,0)+vx(0,0),2);
            
            mu0 = .5*pow(AFace, -1.0/3.0);

            mu =  mu0*( pow( (eps + SRI), -1.0/3.0) + delta ); 

            // if(iv[0]==20 && iv[1]==18)
            // {
            //     printf("qwer");
            // }

            if(constMu)
            {
                if(grd)
                {
                    etaFab(iv,0) = muCGrd * h(0,0);
                }
                else
                {
                    etaFab(iv,0) = muCFlt * h(0,0);
                }
            }
            else
            {
                etaFab(iv,0) = mu* h(0,0);
            }

        }


        

    }
    beta.exchange();
    eta.exchange();




}




#include "NamespaceFooter.H"