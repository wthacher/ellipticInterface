#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef __ROMBERG_H_
#define __ROMBERG_H_

#include "REAL.H"
#include "Vector.H"
#include <cmath>
#include <math.h>
#include <iostream>
#include <fstream>

#include "NamespaceHeader.H"
/*TO DO:
1. Make function that calculates moments on other side of EB? -
this will just involve grabbing the other corners
3. What if interpolated GL creates a third intersection?--> i guess AMR
4. 
*/

class momentConcept
{
    // calculates moment on shape defined by EB and cell corners
    // eval(n) means split up our EB into 2^n segments
    virtual Real eval(const int n);
    virtual Real areaEval(const int n, const int dir);
    virtual ~momentConcept() {};
};


template <class moment> class romberg //template just assumes you have an eval(N,x,y,p,q) function
{
public:
//
//uses romberg integration to approximate p,q moment
    Real
    integrate(moment a_func, int p, int q, Real a_tol = 1e-12, int a_max_iter = 8)
    {
        CH_TIME("romberg::integrate");
        Real err = 2*a_tol;
        Real errPrevious = err; //we dont want to run into case where we stop early accidentally
        Real sln = 0.0;
        
        Vector<Real> row(a_max_iter,0.0);//thats the longest a row could be
        Vector<Vector<Real>> R (a_max_iter,row);//square matrix
        //initialize with n=0
        R[0][0] = a_func.eval(0,p,q);
        int n=1;
        Real c;
        while(n<a_max_iter&& (errPrevious>a_tol || err>a_tol) )//we need to have this iterate and last iterate ok
        {
            R[n][0] = a_func.eval(n,p,q);
            int m = 1;
            while (m<=n)
            {
                c = pow(2.0,m+1);
                R[n][m] = (c*R[n][m-1] - R[n-1][m-1])/(c-1.0);
                m = m+1;
            }
            errPrevious = err;
            err = abs(R[n][n] - R[n][n-1])/abs(R[n][n-1]);
            sln = R[n][n];
            n= n+1;

        }
        //std::cout<<"Iterations: "<<n-1<<"\n";
        return sln;
         
    };

    Real
    integrateArea(moment a_func, int p, int q, int dir, Real a_tol = 1e-15, int a_max_iter = 12)
    {
        CH_TIME("romberg::integrateArea");
        Real err = 2*a_tol;
        Real errPrevious = err;
        Real sln = 0.0;
        
        Vector<Real> row(a_max_iter,0.0);//thats the longest a row could be
        Vector<Vector<Real>> R (a_max_iter,row);//square matrix
        //initialize with n=0
        R[0][0] = a_func.areaEval(0,p,q,dir);
        int n=1;
        Real c;
        while(n<a_max_iter&& (errPrevious>a_tol || err>a_tol) ) 
        {
            R[n][0] = a_func.areaEval(n,p,q,dir);
            int m = 1;
            while (m<=n)
            {
                c = pow(2.0,m+1);
                R[n][m] = (c*R[n][m-1] - R[n-1][m-1])/(c-1.0);
                m = m+1;
            }
            errPrevious = err;
            err = abs(R[n][n] - R[n][n-1])/abs(R[n][n-1]);
            sln = R[n][n];
            n= n+1;

        }
        //std::cout<<"Iterations: "<<n-1<<"\n";
        // if(sln == 1)
        // {
        //     n+=1;
        // }
        return sln;
         
    };

    

    //approximates all moment combos up to P,Q. returns in order: p=0, all q in order, p=1, all q etc
    Vector<Real>
    integrate_all(moment a_func, int P, int Q, Real dx, Real a_tol = 1e-12, int a_max_iter = 8)
    {
        Real err = 2*a_tol;
        Vector<Real> out((P+1)*(Q+1),0.0);
        Real sln;
        
        Vector<Real> row(a_max_iter,0.0);//thats the longest a row could be
        Vector<Vector<Real>> R (a_max_iter,row);//square matrix
        
        for(int p =0;p<=P;p++)
        {
            for(int q = 0;q<=Q;q++)
            {
                err = 2*a_tol;//reset error!
                R.assign(row);
                R[0][0] = a_func.eval(0,p,q);
                int n=1;
                Real c;
                while(err>a_tol && n<=a_max_iter)
                {
                    R[n][0] = a_func.eval(n,p,q);
                    int m = 1;
                    while (m<=n)
                    {
                        c = pow(4.0,m);
                        R[n][m] = (c*R[n][m-1] - R[n-1][m-1])/(c-1.0);
                        m = m+1;
                    }
                    err = abs(R[n][n] - R[n][n-1])/abs(R[n][n-1]);
                    sln = R[n][n];
                    n= n+1;

                }
                out[p*(Q+1)+q] = sln*pow(dx,q+p+2);//need to scale back from unit dx
            }
        }
        return out;

    };

};


#include "NamespaceFooter.H"

#endif
