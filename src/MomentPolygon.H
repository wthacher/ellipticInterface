#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef __MOMENT_POLYGON_H_
#define __MOMENT_POLYGON_H_


#include "REAL.H"
#include "Vector.H"
#include <cmath>
#include <math.h>
#include <assert.h> 
#include "ImplicitFunction.H"

#include <iostream>

//#include "Intersection_utils.H"

#include "NamespaceHeader.H"
//Implementing this as a class so it can be passed as a function to Romberg

//helper functions
inline int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}

inline int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}

inline Real r8_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in R8 arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, the values of N and K.
//
//    Output, double R8_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  Real value;

  mn = i4_min ( k, n - k );

  if ( mn < 0 )
  {
    value = 0.0;
  }
  else if ( mn == 0 )
  {
    value = 1.0;
  }
  else
  {
    mx = i4_max ( k, n - k );
    value = ( Real ) ( mx + 1 );

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( Real ) ( mx + i ) ) / ( Real ) i;
    }
  }

  return value;
}


//just uses implicit function info
class IFPolygon
{
public:
    ImplicitFunction* IFPtr;

    Vector<Real> xB;//intersection pairs, right then left
    Vector<Real> yB;//corresponding y pts
    Vector<Real> xc; //cell corners -1, 2 or 3 points
    Vector<Real> yc; ///corresponding y pts
    Real xCtd, yCtd; //so that we can translate info from IF
    Real dx;

    bool xdiv; //divide using x axis?

    IFPolygon(){};

    IFPolygon(Vector<Real>& a_xB, Vector<Real>& a_yB, 
                      Vector<Real>& a_xc, Vector<Real>& a_yc,
                      Real a_xCtd, Real a_yCtd, Real a_dx,
                      ImplicitFunction& a_IF)
    {
        xB  = a_xB; 
        yB =a_yB; 
        xc = a_xc;
        yc = a_yc;
        xCtd = a_xCtd;
        yCtd = a_yCtd;
        dx = a_dx;
        IFPtr = &a_IF;

        xdiv=false;
        if(abs(xB[0]-xB[1]) > abs(yB[0]-yB[1])) {xdiv=true;}
    }

    void define(Vector<Real>& a_xB, Vector<Real>& a_yB, 
                      Vector<Real>& a_xc, Vector<Real>& a_yc,
                      Real a_xCtd, Real a_yCtd, Real a_dx,
                      ImplicitFunction& a_IF)
    {
        xB  = a_xB; 
        yB =a_yB; 
        xc = a_xc;
        yc = a_yc;
        xCtd = a_xCtd;
        yCtd = a_yCtd;
        dx = a_dx;
        IFPtr = &a_IF;

        xdiv=false;
        if(abs(xB[0]-xB[1]) > abs(yB[0]-yB[1])) {xdiv=true;}
    }

    void split(int n, Vector<Real>& x_n, Vector<Real>& y_n)
    {
        CH_TIME("IFPolygon::split");
        if(xdiv) //divide along x axis, meaning fix x and solve for y
        {
            Real N = pow(2.0,n);//N+1 represents number of points along EB INCLUDING endpoints
            x_n.resize(N+1+xc.size() );
            y_n.resize(N+1+xc.size() );
            Real dif = xB[0] - xB[1]; //right minus left
            Real root,x1,y1,x2,y2;
            //assert(dif != 0.0); //will have to figure out how to deal with this later
            for (int i =0;i<=N;i++) //this is counter clockwise. fills in N+1 points to create N segments
            {
                x_n[i] = xB[0] - (dif/N)*i;  //right boundary minus step size down to left boundary
            }
            for(int j=N+1;j<x_n.size();j++) //fill in our remaining x points
            {
                x_n[j] = xc[j-N-1]; //starts from zero
                y_n[j] = yc[j-N-1];
            }
            y_n[0] = yB[0];
            y_n[N] = yB[1];

            for (int i =1;i<=N-1;i++) //this is counter clockwise. fills in N+1 points to create N segments
            {
                x1 = x_n[i]*dx + xCtd;
                y1 = -.5*dx + yCtd;
                x2 = x_n[i]*dx + xCtd;
                y2 = .5*dx + yCtd;

                IFPtr->secantMethod(x1,y1,x2,y2);

                y_n[i] = (y2 - yCtd)/dx;
                
            }
        }
        else
        {
            Real N = pow(2.0,n);//N+1 represents number of points along EB INCLUDING endpoints
            x_n.resize(N+1+xc.size() );
            y_n.resize(N+1+xc.size() );
            Real dif = yB[0] - yB[1]; //right minus left
            Real root,x1,y1,x2,y2;
            //assert(dif != 0.0); //will have to figure out how to deal with this later
            for (int i =0;i<=N;i++) //this is counter clockwise. fills in N+1 points to create N segments
            {
                y_n[i] = yB[0] - (dif/N)*i;  //right boundary minus step size down to left boundary
            }
            for(int j=N+1;j<x_n.size();j++) //fill in our remaining x points
            {
                y_n[j] = yc[j-N-1]; //starts from zero
                x_n[j] = xc[j-N-1];
            }
            x_n[0] = xB[0];
            x_n[N] = xB[1];

            for (int i =1;i<=N-1;i++) //this is counter clockwise. fills in N+1 points to create N segments
            {
                y1 = y_n[i]*dx + yCtd;
                x1 = -.5*dx + xCtd;
                y2 = y_n[i]*dx + yCtd;
                x2 = .5*dx + xCtd;

                IFPtr->secantMethod(x1,y1,x2,y2);

                x_n[i] = (x2 - xCtd)/dx;
            }
        }
    };

    Real eval(int n, int p, int q)
    {
        CH_TIME("IFPolygon::eval");
        Vector<Real> x,y; 
        split(n,x,y);
        
        x.push_back(xB[0]); //so it wraps around
        y.push_back(yB[0]); //so it wraps around
        int nSeg = x.size() - 1; //number of segments
        Real tot = 0;
        Real dx, dy, xk, yk;
        for(int s =0;s<nSeg;s++)
        { 
            xk = x[s];
            yk = y[s];
            dx = x[s+1] - xk;
            dy = y[s+1] - yk;
            for(int i =0; i<= p+1; i++)
            {
                for(int j = 0;j<=q;j++)
                {
                    tot+= r8_choose(p+1,i) * r8_choose(q,j) * pow(xk,i) * pow(yk,j)*pow(dx,p+1-i) * pow(dy, q+1-j) * (1.0/(p+2-i+q-j));
                }
            }
        }

        tot*=1.0 / (p+1);
        return tot;
    }

    Real areaEval (int n,const int p,const int q, const int dir)
    {   
        CH_TIME("IFPolygon::areaEval");
        Vector<Real> x,y; 
        split(n,x,y);
        //cut off the corners
        for(int i =0;i<xc.size();i++)
        {
            x.pop_back();
            y.pop_back();
        }

        Real tot = 0;
        Real xk, yk, dx, dy, dc, t;

        int nSeg = x.size() - 1;

        for(int s =0; s< nSeg; s++)
        {
            xk = x[s];
            dx = x[s+1] - x[s];
            yk = y[s];
            dy = y[s+1] - y[s];
            dc = pow(dx*dx + dy*dy, .5);

            for(int i =0; i<= p; i++)
            {
                for(int j = 0;j<=q;j++)
                {
                    t = r8_choose(p,i) * r8_choose(q,j) * pow(xk,i) * pow(yk,j)*pow(dx,p-i) * pow(dy, q-j) * (1.0/(p+1-i+q-j));
                    if(dir ==0) //x normal - multiply by dy/dc * dc
                    {
                        tot+= t*dy;
                    }
                    else if(dir ==1) //y normal - multiply by -dx/dc *dc
                    {
                        tot+= -1.0*t*dx;
                    }
                    else //do area integral, multiply by dC
                    {
                        tot+= t*dc;
                    }
                    
                    
                }
            }
        }

        return tot;
    };




};

#include "NamespaceFooter.H"

#endif
