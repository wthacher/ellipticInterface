#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBSolveUtilities.H"

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "NamespaceHeader.H"



void makeRegularCellStencil(LAPACKMatrix& SBeta,
                            Vector <LAPACKMatrix>& SEta,
                            int order)
{
    int Q = order;
    int n = ((Q+1)*(Q+2) )/2;
    int r = (Q+1)/2;
    int row=0;
    int comp;
    int nCells = (2*r+1)*(2*r+1);
    Vector<Real> moments(n,0); Vector<Real> grads(2,0);
    Vector<Real> hCoefVect(n,0); Vector<Real> zbCoefVect(n,0); 
    Real mom;

    
    LAPACKMatrix SBeta1, SEtaTemp1, SEtaTemp2;
    IntVect position,sten;

    //can skip cells that are further than L1 distance of r-1 away

    Real dist, weight,w;
    w = Q + 1;

    LAPACKMatrix W, MvU2, MvV2;

    W.define(nCells, nCells); W.setVal(0);

    Real sx,sy;

    //bc of symmetry all you need is two orders less than Q for this term
    int bQ = Q - 2;

    LAPACKMatrix Mv( nCells, n); Mv.setVal(0); //interpolation matrix of VOLUMES
    LAPACKMatrix Mp( nCells, n); Mp.setVal(0); //interpolation matrix of POINTS
    LAPACKMatrix MBeta(n,n);

    for(int l = -r;l<r+1;l++)
    {
        for(int c = -r;c<r+1;c++) //c++ hahhaa
        {
            if(abs(l) + abs(c) > r-1 ){row+=1; continue;} //dont need these boys
            for(int px =0; px<=Q; px++)
            {
                for(int py = 0; py<=Q-px;py++)
                {
                    if(px + py > Q - 2 || ( ( (px + py)==(Q-2) ) && ( (px%2)!=0 || (py%2)!=0) ) ){continue;}
                    comp=  px*Q + (3*px - px*px)/2 + py;
                    mom = (pow(c+.5,px+1) - pow(c-.5,px+1) ) / (px+1);
                    mom*= (pow(l+.5,py+1) - pow(l-.5,py+1) ) / (py+1);
                    Mv(row,comp) = mom;

                    Mp(row,comp) = pow(c,px) * pow(l,py);
                }
            }
            row+=1;
        }
    }

    Mv.pseudoInvertUsingSVD(10,1e-10);
    Mp.pseudoInvertUsingSVD(10,1e-10);
    Mp.transpose();

    row =  0; int pxx,pyy;

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

    MBeta.transpose();

    multiply(SBeta1, MBeta, Mv);
    multiply(SBeta, Mp, SBeta1); //multiply from the left by point values of beta, and on the right by cell averaged values of u


    //now do the flux stencils

    Real nx,ny,sgn;
    LAPACKMatrix MF;
    LAPACKMatrix MpEta, MvU;
  
    Real xCtd, yCtd;

    for(int i=0;i<1;i++)
    {
        SEta[i].define(nCells,nCells);
        SEta[i].setVal(0);
    }

    
    
    for(int d = 0;d<2;d++)
    {
        for(int s=0;s<2;s++)
        {
            sgn = s==0 ? -1:1; //for hi and lo side of CELL
            nx = ((d+1)%2)*sgn;//if d is 0, nx is 1
            ny = d*sgn; //if d is 0, ny is 0

            MF.define(n, n); MF.setVal(0);

            MvU.define(nCells,n);   MvU.setVal(0);

            MpEta.define(nCells,n); MpEta.setVal(0);
            W.define(nCells, nCells); W.setVal(0);

            //build up interpolation matrices. dont need the term such that p_d >= Q, which zeros out columns
            //also dont need some of the points; this zeros out rows
            // for d=0,s=0, dont need c = r.   for d=0,s=1, dont need c = -r.
            // for d=1,s=0, dont need l = r.   for d=1,s=1, dont need l = -r.

            row=0;
            for(int l = -r;l<r+1;l++)
            {
                for(int c = -r;c<r+1;c++) 
                {
                    //skip these rows bc they are not in the stencil for this face
                    if(d==0 && s==0 && c==r){row+=1; continue;}
                    if(d==0 && s==1 && c==-r){row+=1; continue;}
                    if(d==1 && s==0 && l==r){row+=1; continue;}
                    if(d==1 && s==1 && l==-r){row+=1; continue;}

                    sx = c; sy=l;

                    dist = sqrt(sx*sx + sy*sy) + 1; //we offset this to center node. now center cell will have weight 1
                    weight= pow(dist,-1*w); //experiment with weighting. Hans says use Q+1

                    W(row,row) = weight;


                    for(int px =0; px<=Q; px++)
                    {
                        for(int py = 0; py<=Q-px;py++)
                        {
                            if(d == 0 && px == Q){continue;} //dont need this columns bc we dont have the rank for them
                            if(d == 1 && py == Q){continue;}

                            comp=  px*Q + (3*px - px*px)/2 + py;
                            mom = (pow(c+.5,px+1) - pow(c-.5,px+1) ) / (px+1);
                            mom*= (pow(l+.5,py+1) - pow(l-.5,py+1) ) / (py+1);
                            MvU(row,comp) = weight*mom;

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


            //SEta += pinv(MpEta)^T MFx^T pinv(W*Mv)*W
            MvU.pseudoInvertUsingSVD(10,1e-10);

            multiply(MvU2, MvU, W);

            MpEta.pseudoInvertUsingSVD(10,1e-10);
            MpEta.transpose();

            MF.transpose();

            multiply(SEtaTemp1, MF, MvU2);
            multiply(SEtaTemp2, MpEta, SEtaTemp1);
            SEta[0] += SEtaTemp2; //L(u)

        }
    }

}


//writes errors for all comps into rows of err matrix
//first column is L1, second is L2, third is max norm
//all of LevelDatas these should be defined on a single box
//they will contain: [uGrd, uFlt, rhsGrd, rhsFlt, 
//betaGrd, betaFlt, etaGrd, etaFlt, psiNodes, grdMoment, ]

void convergenceTest(const LevelData<FArrayBox>& solCoarse,
                     const LevelData<FArrayBox>& solFine,
                     Real dxC, Real dxF, int Q,
                     LAPACKMatrix& err)
{
    int nComps = 4;
    int mom=9; int psi=8;
    Real dxC2 = dxC * dxC; Real dxF2 = dxF * dxF;

    err.define(nComps, 3); err.setVal(0);
    Real ratio = dxC / dxF;

    for(DataIterator dit = solCoarse.dataIterator();dit.ok();++dit)
    {
        const FArrayBox& solFabC = solCoarse[dit()];

        Box validC = (solCoarse.disjointBoxLayout() )[dit()];
        
        for(DataIterator dit = solFine.dataIterator();dit.ok();++dit)
        {
            const FArrayBox& solFabF = solFine[dit()];

            Box validF = (solFine.disjointBoxLayout() )[dit()];

            
            FArrayBox dif(validC, 2*nComps);
            dif.setVal(sqrt(-1.0)); //nans are good for plotting

            Real volF, volC;
            Vector<Real> coarseSol(2*nComps,0);
            Vector<Real> fineSol(2*nComps, 0);

            int nCut = 0;
            int c;

            Box fineBox; IntVect iv;

            for(BoxIterator bit(validC); bit.ok();++bit)
            {
                nCut += 1;
                iv = bit();


                if(solFabC(iv,mom) != -1) //cut cell. mom is the GRD volume fraction
                {
                    nCut += 1;
                    volC = solFabC(iv,mom);
                    fineBox.define(ratio*iv,ratio*iv + (ratio-1)*IntVect::Unit);
                    fineSol.assign(0); 

                    for(BoxIterator stenIt(fineBox); stenIt.ok(); ++stenIt)
                    {
                        if(solFabF(stenIt(),mom) != -1) //cut cell on fine level
                        {
                            volF = solFabF(stenIt(),mom );
                            for(int i = 0;i<nComps;i++)
                            {
                                fineSol[2*i] += solFabF(stenIt(), 2*i) * volF * dxF2;
                                fineSol[2*i+ 1] += solFabF(stenIt(), 2*i + 1) * (1.0-volF) * dxF2;
                            }
                        }
                        else
                        {
                            if(solFabF(stenIt(), psi) > 0)
                            {
                                c = 0;
                            }
                            else
                            {
                                c=1;
                            }
                            for(int i = 0;i<nComps;i++)
                            {
                                fineSol[2*i+c] += solFabF(stenIt(), 2*i + c) * dxF2;
                            }

                        }
                    }

                    for(int i=0;i<nComps;i++)
                    {
                        coarseSol[2*i] = solFabC(iv,2*i) * volC * dxC2;
                        coarseSol[2*i + 1] = solFabC(iv,2*i + 1) * (1.0-volC) * dxC2;
                    }

                    for(int i=0;i<nComps;i++)
                    {
                        err(i,0) += abs(coarseSol[2*i] - fineSol[2*i]) ;
                        err(i,0) += abs(coarseSol[2*i+1] - fineSol[2*i+1]) ;

                        err(i,1) += pow(coarseSol[2*i] - fineSol[2*i], 2);
                        err(i,1) += pow(coarseSol[2*i+1] - fineSol[2*i+1], 2);

                        err(i,2) = max(err(i,2), abs(coarseSol[2*i] - fineSol[2*i])/(dxC2) );
                        err(i,2) = max(err(i,2), abs(coarseSol[2*i+1] - fineSol[2*i+1])/(dxC2) );

                        dif(iv,2*i) = (fineSol[2*i] - coarseSol[2*i] ) / (volC * dxC2);
                        dif(iv,2*i + 1) = (fineSol[2*i+1] - coarseSol[2*i+1]) / ((1.0-volC) * dxC2);
                    
                    }

                    

                }
                else //full celll
                {
                    if(solFabC(iv,psi) > 0) //grd
                    {
                        c = 0;
                    }
                    else //flt
                    {
                        c = 1;
                    }
                    fineBox.define(ratio*iv,ratio*iv + (ratio-1)*IntVect::Unit);
                    fineSol.assign(0); 

                    for(int i=0;i<nComps;i++)
                    {
                        coarseSol[i] = solFabC(iv, 2*i + c);
                    }

                    for(BoxIterator stenIt(fineBox); stenIt.ok(); ++stenIt)
                    {
                        for(int i = 0;i<nComps;i++)
                        {
                            fineSol[i] += (1.0)/(ratio*ratio) * solFabF(stenIt(), 2*i + c);
                        }
                    }

                    for(int i=0;i<nComps;i++) //L^1 err is | (coarseAvg - fineAvg)*dxC^2 |
                    {
                        err(i,0) += abs(fineSol[i] - coarseSol[i])*dxC2;
                        err(i,1) += pow( (fineSol[i] - coarseSol[i])*dxC2,2);
                        err(i,2) = max(abs(fineSol[i] - coarseSol[i]), err(i,2) );

                        dif(iv,2*i + 0) = (fineSol[i] - coarseSol[i]); //average error over this cell
                        dif(iv,2*i + 1) = (fineSol[i] - coarseSol[i]);
                    }


                }


                
            }

            int bd;
            bd=0;
            
        }
    }

    for(int i=0;i<nComps;i++)
    {
        err(i,1) = sqrt(err(i,1) ); //L^2 norm
    }
}


#include "NamespaceFooter.H"