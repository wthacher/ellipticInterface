#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef __SOLUTIONSETUP_H_
#define __SOLUTIONSETUP_H_

#include "REAL.H"
#include "Vector.H"
#include <cmath>
#include <math.h>

#include "NamespaceHeader.H" 

//use this to fill in RHS and exact solution
//will probably need derivative info in order to take integrals
//this will get passed as an argument to set things up

//NOTE - THINK ABOUT UNITS!!!!!!!!!!! WHAT IS EASIEST WAY TO DEAL WITH THEM?????/
//MOMENTS ARE IN DX = 1 UNITS



class RealFunc
{
public:
    LAPACKMatrix pM,moments;
     //pM is for making coefficients, and M is for full cell integral

    RealFunc()
    {
    }

    virtual Real val(Real x, Real y)
    {
        return 0;
    }

    virtual Real derivVal(Real x, Real y, int dx, int dy)
    {
        return 0;
    }

    void makeMatrices(int Q)
    {
        int r = ( Q + 1 )/2;
        int n=  ((Q+1)*(Q+2) )/2;

        
        pM.define((2*r+1)*(2*r+1), n);

        Real mom;
        int comp;

        int row=0;
        for(int l = -r;l<r+1;l++)
        {
            for(int c = -r;c<r+1;c++) 
            {
                comp=0;
                for(int px =0; px<=Q; px++)
                {
                    for(int py = 0; py<=Q-px;py++)
                    {
                        mom = pow(c, px) * pow(l, py);
                        pM(row,comp) = mom;
                        comp+=1;
                    }
                }
                row+=1;
            }
        }

        pM.pseudoInvertUsingSVD(10,1e-10);

        moments.define(1,n); //moments for normal square cell
        int c=0;

        for(int px =0; px<=Q; px++)
        {
            for(int py = 0; py<=Q-px;py++)
            {
                mom = (pow(.5,px+1) - pow(-.5,px+1) ) / (px+1);
                mom *= (pow(.5,py+1) - pow(-.5,py+1) ) / (py+1);
                moments(0,c) = mom;
                c+=1;
            }
        }

    }

    //generate Taylor Series coefficients up to order Q about xc,yc
    //generates this with unit spacing in M with center (0,0)
    void generateCoefs(int Q, Real xc, Real yc, Real dx,
                        LAPACKMatrix& coefs)
    {
        //interpolate uniform grid with spacing dx - there is probably a more intelligent option
        int r = ( Q + 1 )/2;
        int n=  ((Q+1)*(Q+2) )/2;

        LAPACKMatrix Phi;
        Phi.define((2*r+1)*(2*r+1), 1);

        int row=0;
        for(int l = -r;l<r+1;l++)
        {
            for(int c = -r;c<r+1;c++) //c++ hahhaa
            {
                Phi(row,0) = val(xc + c*dx, yc + l*dx);
                row+=1;
            }
        }

        multiply(coefs, pM, Phi); //pinv(M) * Phi
    }

    //integral over square cell. x1 y1 is lower left node
    Real cellIntegral(int Q, Real x1, Real y1, Real dx)
    {
        LAPACKMatrix v,coefs;
        generateCoefs(Q, x1+.5*dx, y1+.5*dx, dx, coefs);
        multiply(v, moments, coefs);
        return v(0,0);
    }

    //integral over region defined by moments - could be area or volume
    //xc yc is the center of the square cell in physical space 
    Real integral(int Q, LAPACKMatrix& a_moments, Real xc, Real yc, Real dx)
    {
        LAPACKMatrix v, coefs;
        generateCoefs(Q, xc,yc, dx, coefs); //gets coefs with spacing dx
        multiply(v, a_moments, coefs);
        return v(0,0);

    }

    //integral of eta grad(u) cdot n over EB
    // form matrix M which has px m_x^q+p-ex + py m_y^q+p-_ey
    Real integrateDuDn(int Q, const LAPACKMatrix& nXmoments, const LAPACKMatrix& nYmoments, 
                            RealFunc& etaFunc, Real xc, Real yc, Real dx)
    {
        LAPACKMatrix etaCoefs, uCoefs,v;
        etaFunc.generateCoefs(Q,xc,yc,dx,etaCoefs);
        etaCoefs.transpose();

        generateCoefs(Q,xc,yc,dx,uCoefs);

        int row=0;
        int n =  ((Q+1)*(Q+2) )/2;
        int comp,compx, compy; Real nXx, nYy;

        LAPACKMatrix MF(n,n); MF.setVal(0);

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

                        nXx = (qx>=1) ? nXmoments(0,compx): 0;

                        nYy = (qy>=1) ? nYmoments(0,compy): 0;

                        MF(row,comp) = qx * nXx + qy * nYy;
                    }
                }
                row+=1;
            }
        }
        
        MF.transpose();

        multiply(v, etaCoefs, MF);
        LAPACKMatrix v2;
        multiply(v2, v, uCoefs);

        return v2(0,0);

    }

    
    
};


class SolutionSetUp
{
public:
    Real x0,y0,dx; //x0, y0 is lower left corner in physical space
    int Q,n;

    SolutionSetUp(){}

    SolutionSetUp(Real a_x0, Real a_y0, Real a_dx, int a_Q)
    {
        x0 = a_x0; y0 =a_y0; dx = a_dx; Q = a_Q;
        n = ( (Q+1)*(Q+2) )/2;
    }

    //if avg = false, fill in point values at centroids
    void setUpFullCellField(LevelData<FArrayBox>& cellPhi, 
                            const LevelData<NodeFArrayBox>& psiNodes,
                            RealFunc& funcIn, RealFunc& funcOut, 
                            bool avg,
                            Real dx)
    {
        for(DataIterator dit(cellPhi.dataIterator()); dit.ok(); ++dit)
        {
            FArrayBox& cellPhiFab = cellPhi[dit()];
            const FArrayBox& psiNodesFab = psiNodes[dit()].getFab();
            Real x1,y1,x2,y2,xc,yc;
            IntVect iv;
            Real dx2 = dx*dx;

            funcIn.makeMatrices(Q);
            funcOut.makeMatrices(Q);

            for(BoxIterator bit(cellPhiFab.box()); bit.ok(); ++bit)
            {
                iv = bit();
                x1 = iv[0] * dx + x0;
                y1 = iv[1] * dx + y0;

                if(psiNodesFab(iv) > 0)
                {
                    if(avg){cellPhiFab(iv) = funcIn.cellIntegral(Q,x1,y1,dx) ;}
                }
                else
                {
                    if(avg){cellPhiFab(iv) = funcOut.cellIntegral(Q,x1,y1,dx) ;}
                }

                //if filling out coefficients, fill in both vals
                if(!avg)
                {
                    cellPhiFab(iv,0) = funcIn.val(x1 + .5*dx, y1 + .5*dx);
                    cellPhiFab(iv,1) = funcOut.val(x1 + .5*dx, y1 + .5*dx);
                }
                
            }
        }
    }

    void setUpCutCellField(LevelData<IVSFAB<Real>>& cellPhiCut, 
                            const LevelData<IVSFAB<Real>>& moments,
                            const LevelData<NodeFArrayBox>& psiNodes,
                            RealFunc& funcIn,
                            RealFunc& funcOut,
                            bool avg,
                            Real dx)
    {
        for(DataIterator dit(moments.dataIterator()); dit.ok(); ++dit)
        {
            IVSFAB<Real>& phiCutFab = cellPhiCut[dit()];
            const FArrayBox& psiNodesFab = psiNodes[dit()].getFab();
            const IVSFAB<Real>& momentsFab = moments[dit()];
            IntVectSet EBCellsFab = momentsFab.getIVS();

            funcIn.makeMatrices(Q);
            funcOut.makeMatrices(Q);

            phiCutFab.define(EBCellsFab, 2); //1 for each phase

            Real xc,yc, xcHi, xcLo, ycHi, ycLo;
            IntVect iv;
            Real dx2 = dx*dx;
            int cx = Q + 1;

            LAPACKMatrix momentsLo(1,n);
            LAPACKMatrix momentsHi(1,n);
            Real cVIn, cVOut, cVInP, cVOutP;

            for(IVSIterator ivIt(EBCellsFab); ivIt.ok(); ++ivIt)
            {
                iv = ivIt();
                xc = iv[0] * dx + x0 + .5*dx;
                yc = iv[1] * dx + y0 + .5*dx;

                for(int i=0;i<n;i++)
                {
                    momentsLo(0,i) = momentsFab(iv,i);
                    momentsHi(0,i) = momentsFab(iv,i + n);
                }

                if(psiNodesFab(iv) > 0)
                {
                    if(avg)
                    {
                        phiCutFab(iv,0) = (funcIn.integral(Q,momentsLo, xc,yc,dx) ) / (momentsLo(0,0) );
                        phiCutFab(iv,1) = (funcOut.integral(Q,momentsHi, xc,yc,dx) ) / (momentsHi(0,0) );
                    }
                    else //get cell centroid
                    {
                        xcLo = dx*(momentsFab(iv,cx) / momentsFab(iv,0) ) + xc;
                        ycLo = dx*(momentsFab(iv,1) / momentsFab(iv,0) ) + yc;

                        xcHi = dx*(momentsFab(iv,cx+n) / momentsFab(iv,n) ) + xc;
                        ycHi = dx*(momentsFab(iv,1+n) / momentsFab(iv,n) ) + yc;

                        phiCutFab(iv,0) = funcIn.val(xcLo, ycLo);
                        phiCutFab(iv,1) = funcOut.val(xcHi, ycHi);


                    }
                }
                else
                {
                    if(avg)
                    {
                        phiCutFab(iv,0) = (funcIn.integral(Q,momentsHi, xc,yc,dx) ) / (momentsHi(0,0) );
                        phiCutFab(iv,1) = (funcOut.integral(Q,momentsLo, xc,yc,dx) ) / ( momentsLo(0,0) );
                        
                    }
                    else //get cell centroid
                    {
                        xcLo = dx*(momentsFab(iv,cx) / momentsFab(iv,0) ) + xc;
                        ycLo = dx*(momentsFab(iv,1) / momentsFab(iv,0) ) + yc;

                        xcHi = dx*(momentsFab(iv,cx+n) / momentsFab(iv,n) ) + xc;
                        ycHi = dx*(momentsFab(iv,1+n) / momentsFab(iv,n) ) + yc;

                        phiCutFab(iv,0) = funcIn.val(xcHi, ycHi);
                        phiCutFab(iv,1) = funcOut.val(xcLo, ycLo);
                    }
                }
                
            }
        }
    }

    //jump is always in - out
    void setUpJumps(LevelData<IVSFAB<Real>>& jumps,
                    const LevelData<IVSFAB<Real>>& moments,
                    const LevelData<NodeFArrayBox>& psiNodes,
                    RealFunc& funcIn,
                    RealFunc& funcOut,
                    RealFunc& etaFuncIn,
                    RealFunc& etaFuncOut,
                    Real dx)
    {
        for(DataIterator dit(moments.dataIterator()); dit.ok(); ++dit)
        {
            IVSFAB<Real>& jumpFab = jumps[dit()];
            const FArrayBox& psiNodesFab = psiNodes[dit()].getFab();
            const IVSFAB<Real>& momentsFab = moments[dit()];
            IntVectSet EBCellsFab = momentsFab.getIVS();

            jumpFab.define(EBCellsFab, 2); //1 for each phase 

            Real xc,yc, xcHi, xcLo, ycHi, ycLo;
            IntVect iv;
            Real dx2 = dx*dx;
            int cx = Q + 1;

            LAPACKMatrix areaMoments(1,n);
            LAPACKMatrix  nXMoments(1,n);
            LAPACKMatrix  nYMoments(1,n);

            Real uJump, duJump;

            funcIn.makeMatrices(Q);
            funcOut.makeMatrices(Q);

            etaFuncIn.makeMatrices(Q);
            etaFuncOut.makeMatrices(Q);

            for(IVSIterator ivIt(EBCellsFab); ivIt.ok(); ++ivIt)
            {
                iv = ivIt();
                xc = iv[0] * dx + x0 + .5*dx;
                yc = iv[1] * dx + y0 + .5*dx;

                
                for(int i=0;i<n;i++)
                {
                    areaMoments(0,i) = momentsFab(iv,i + 2*n);
                    nXMoments(0,i) = momentsFab(iv,i + 3*n);
                    nYMoments(0,i) = momentsFab(iv,i + 4*n);
                }
                
                uJump = funcIn.integral(Q,areaMoments, xc, yc,dx) - funcOut.integral(Q,areaMoments, xc, yc,dx);

                jumpFab(iv,0) = uJump / areaMoments(0,0); //average over face

                duJump = funcIn.integrateDuDn(Q,nXMoments, nYMoments, etaFuncIn, xc, yc,dx) ;
                duJump -= funcOut.integrateDuDn(Q,nXMoments, nYMoments, etaFuncOut, xc, yc,dx) ;

                jumpFab(iv,1) = duJump / areaMoments(0,0);                
                
            }
        }
    }

    //these are moments in boundary cut cells for irregular geometry
    void setUpCutCellBCs(LevelData<IVSFAB<Real>>& BCs,
                    const LevelData<IVSFAB<Real>>& moments,
                    const LevelData<NodeFArrayBox>& psiNodes,
                    RealFunc& func,
                    Real dx)
    {

    }

    //bcType 0 is dirchlet, 1 is neumann
    //how do we want to deal with corner cells?
    void setUpFullCellBCs(LevelData<IVSFAB<Real>>& BCs,
                    RealFunc& func,
                    int bcType,
                    Real dx)
    {
        
    }
};

//different real functions for different probs:

class constFunc : public RealFunc
{
public:
    Real c;

    constFunc(Real a_c){c =  a_c;}
    Real val(Real x, Real y)
    {
        return c;
    }
};

//u = sin^2(wx) * sin^2(wy)
class sin2Func : public RealFunc
{
public:
    Real w = M_PI/2.0;
    
    sin2Func(){}
    
    Real val(Real x, Real y)
    {
        return sin(w * x) * sin(w * x) * sin(w * y) * sin(w * y) ; 
        //return 1 + x;// + y +x*y +x*x + y*y;// + y;
    }
};

//f = -beta u + eta nabla^2 u
//with u = sin^2(wx) * sin^2(wy)
//beta, eta const
class sin2RhsFunc : public RealFunc
{
public:
    Real beta, eta;
    Real w = M_PI/2.0;


    sin2RhsFunc(Real a_beta, Real  a_eta)
    {
        beta = a_beta; 
        eta = a_eta; 
    }

    Real val(Real x, Real y)
    {
        Real sin2x = sin(w*x) * sin(w*x); Real sin2y = sin(w*y) * sin(w*y);
        Real cos2x = cos(w*x) * cos(w*x); Real cos2y = cos(w*y) * cos(w*y);


        Real z =  -1 * beta * (sin2x * sin2y ) ; 
        z += eta * (2*w*w * (cos2x - sin2x ) * sin2y );
        z += eta * (2*w*w * (cos2y - sin2y ) * sin2x );

        return z;
    }
};

class polyFunc : public RealFunc
{
public:
    int Q,n;
    Vector<Real> m_coefs;

    polyFunc(){}

    virtual ~polyFunc(){}

    polyFunc(int a_Q, const Vector<Real>& a_coefs)
    {
        Q = a_Q;
        n = ((Q+1)*(Q+2) )/2;
        m_coefs = a_coefs;
    }

    Real val(Real x, Real y)
    {
        Real z=0;
        int c=0;
        for(int px =0; px<=Q; px++)
        {
            for(int py = 0; py<=Q;py++)
            {
                z+= m_coefs[c] * pow(x,px) * pow(y,py);
                c+=1;
            }
        }
        return z;
    }

    Real derivVal(Real x, Real y ,int dx, int dy)
    {
        Real z=0;
        int c=0;
        int cx, cy;
        for(int px =0; px<=Q; px++)
        {
            for(int py = 0; py<=Q;py++)
            {
                if( px >= dx && py >= dy)
                {
                    cx = 1; cy=1;
                    for(int i=0; i<dx;i++ ){cx*= (px-i);}
                    for(int i=0; i<dy;i++ ){cy*= (py-i);}
                    z += m_coefs[c] * cx * pow(x,px - dx) * cy* pow(y,py - dy);
                }
                
                c+=1;
            }
        }
        return z;
    }

        
};

//f = -beta u + eta nabla^2 u
//with u = P(x,y)
//beta, eta const
class RhsFunc : public RealFunc
{
public:
    RealFunc *betaFunc, *etaFunc, *uFunc;

    RhsFunc(){}

    virtual ~RhsFunc(){}

    RhsFunc(RealFunc *a_betaFunc, RealFunc *a_etaFunc, RealFunc *a_uFunc)
    {
        betaFunc = a_betaFunc;
        etaFunc = a_etaFunc;
        uFunc = a_uFunc;
    }

    Real val(Real x, Real y)
    {
        Real z=0;
        
        z -= betaFunc->val(x,y);
        z += (etaFunc->derivVal(x,y,1,0) ) * ( uFunc->derivVal(x,y,1,0) );
        z += (etaFunc->val(x,y) ) * ( uFunc->derivVal(x,y,2,0) );

        z += (etaFunc->derivVal(x,y,0,1) ) * ( uFunc->derivVal(x,y,0,1) );
        z += (etaFunc->val(x,y) ) * ( uFunc->derivVal(x,y,0,2) );

        return z;
    }

};

class gibouUFunc : public RealFunc
{
public:
    gibouUFunc(){}

    Real val(Real x, Real y)
    {
        Real z = (16.0 * pow((y-x)/5.0,5) - 20.0 * pow(y-x/3.0, 3) + 5.0 * ((y-x)/3.0 ) );
        z *= log(x + y + 3);
        return z;
    }
};

//linear combination of basis functions Re(e^(i nx w x) e^(i ny w y) )
//where nx, ny from -N:N
class periodicFunc : public RealFunc
{
public:
    Vector<Real> coef;
    int N;
    Real w;

    periodicFunc(){}

    virtual ~periodicFunc(){}

    periodicFunc(Vector<Real>& a_coef, int a_N , Real a_w)
    {
        coef = a_coef; N =a_N; w = a_w;
    }
    Real val(Real x, Real y)
    {
        Real z = 0;
        int c=0;
        for(int nx=-N; nx <= N; nx++)
        {
            for(int ny=-N; ny <= N; ny++)
            {
                z += coef[c]* (cos(w*nx*x)*cos(w*ny*y) - sin(w*nx*x)*sin(w*ny*y) );
                c++;
            }
        }

        return z;
    }

    Real derivVal(Real x, Real y, int dx, int dy)
    {
        Real z = 0;
        Real c1,c2,s1,s2;
        int c=0;
        for(int nx=-N; nx <= N; nx++)
        {
            for(int ny=-N; ny <= N; ny++)
            {
                if(dx == 0){c1 = cos(w*nx*x); s1 = sin(w*nx*x);}
                else if(dx == 1){c1 = w*nx*-1*sin(w*nx*x);  s1 = w*nx*cos(w*nx*x);}
                else{c1 = w*w*nx*nx*-1*cos(w*nx*x); s1 = w*w*nx*nx*-1*sin(w*nx*x); }

                if(dx == 0){c2 = cos(w*ny*y); s2 = sin(w*ny*y);}
                else if(dx == 1){c2 = w*nx*-1*sin(w*ny*y);  s2 = w*nx*cos(w*ny*y);}
                else{c2 = w*w*ny*ny*-1*cos(w*ny*y); s2 = w*w*ny*ny*-1*sin(w*ny*y); }

                z+= coef[c]*( c1*c2 - s1*s2 );
                c++;
            }
        }

        return z;
    }

};

#include "NamespaceFooter.H"

#endif
