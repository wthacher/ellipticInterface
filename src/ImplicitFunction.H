#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef __IMPLICITFUNCTION_H_
#define __IMPLICITFUNCTION_H_

#include "REAL.H"
#include "Vector.H"
#include <cmath>
#include <math.h>


#include "NamespaceHeader.H" 

/*This class wraps an ivsfab<real> which represents a stencil. It is used in 2 ways:
1) as a stencil for gradients on full and partial faces
2) as a stencil for flux divergence on EB and tagged cells

The 'add' functions modify the ivs, but not the data
The 'add_' functions where _ is some other data structure adds that type of data to the stencil
The 'add_ToFlux' functions adds VTO components using gradients to the flux divergence
 */
class ImplicitFunction
{
public:

    Real tol; 
    int maxIter;

    ImplicitFunction()
    {
        tol = 1.0e-14;
        maxIter = 30;
    }

    virtual Real val(Real x, Real y)
    {
        return 0;
    }

    //find root between (x1,y1) and (x2,y2)
    //answer goes back into x2,y2
    void secantMethod(Real x1, Real y1, Real& x2, Real& y2 )
    {
        CH_TIME("IF:SecantMethod");
        Real xn, yn, xn1, yn1; //x_n, x_{n+1}
        xn = x1; yn = y1; 
        xn1 = x2; yn1 = y2;

        Real dx, dy, rho, f, f1;

        Real err = 2*tol;

        int iter=0;
        
        //secant method:
        //dx = x_{n+1} - x_n
        //rho = f(x_{n+1}) / (f(x_{n+1}) - f(x_n) )
        //x_{n+1} = x_{n+1} - rho * dx

        dx = xn1 - xn;
        dy = yn1 - yn;
        err = sqrt( dx*dx + dy*dy  );

        while(err > tol && iter < maxIter)
        {
            f = val(xn,yn);
            f1 = val(xn1, yn1);

            rho = f1 / (f1-f);

            xn = xn1;
            yn = yn1;

            xn1 = xn1 - rho * dx;
            yn1 = yn1 - rho * dy;

            dx = xn1 - xn;
            dy = yn1 - yn;

            err = sqrt( dx*dx + dy*dy  ) ;

            iter = iter+1;
        }

        x2 = xn1; y2 = yn1;
    }

};

//psi(x,y) = 1 - a(x-x0)^2 - b(y-y0)^2
class EllipsoidIF : public ImplicitFunction
{
public:
    Real a,b;
    Real x0, y0;

    EllipsoidIF()
    {
        a = b =1.0;
        x0 = y0 = 0;
    }

    EllipsoidIF(Real a_a, Real a_b, Real a_x0, Real a_y0)
    {
        a = a_a;
        b = a_b;
        x0 = a_x0;
        y0 = a_y0;
    }

    Real val(Real x, Real y)
    {

        return 1.0 - a * (x-x0) * (x-x0) - b * (y-y0) * (y-y0);
    }


};

//from https://arxiv.org/pdf/1905.08718.pdf 
class annulusIF : public ImplicitFunction
{
public:
    Real r0 = .483;
    Real p1[3] = {3, .1, .5};
    Real p2[3] = {4, -.1, 1.8};
    Real p3[3] = {7, .15, 0};

    annulusIF()
    {

    }

    Real val(Real x, Real y)
    {
        Real z1 = 0;
        z1 = p1[1] * cos( p1[0] * (atan2(y,x) - p1[2]) ) ;
        z1 += p2[1] * cos( p2[0] * (atan2(y,x) - p2[2]) ) ;
        z1 += p3[1] * cos( p3[0] * (atan2(y,x) - p3[2]) ) ;

        return sqrt(x*x + y*y) - r0 *(1+ z1);
    }

};

//from crockett paper
class rhodoneaIF : public ImplicitFunction
{
public:
    Real r0 = .5;
    Real r1 = .2;
    Real w = 5;
    Real xc = sqrt(2.0) / 20;
    Real yc = sqrt(2.0) / 20;

    rhodoneaIF(){}

    Real val(Real x, Real y)
    {
        //convert to polar
        Real r = sqrt( (x-xc)*(x-xc) + (y-yc)*(y-yc) );
        Real theta = atan2(y-yc, x-xc);

        return -r + (r0 + r1 * sin(w * theta) );
    }
};

//so that we can have no EB
class constIF : public ImplicitFunction
{
public:
    constIF(){}

    Real val(Real x, Real y)
    {
        return 1;
    }
};

#include "NamespaceFooter.H"

#endif