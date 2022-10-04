
#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "IntersectionUtils.H"
#include "math.h"
#include <iomanip>

#include "NamespaceHeader.H"



//this func gives correct pts for INTEGRAL CONTAINING LOWER LEFT NODE AS BOUNDARY POINT
void
get_points(Vector<Real>& x, Vector<Real>& y, Vector<Real>& x_int, Vector<Real>& y_int, Real dx)
{
    CH_TIME("get_points");
    Real y_l = y[0];
    Real x_l = x[0];
    //Real dx =1;
    y.clear(); x.clear(); //we'll refill these
    assert(x_int.size()==2);

    //re order so we everything is counterclockwise
    //two special cases if we have x0 > x1 and y1 = yl,y0>yl
    if( ( x_int[0]<x_int[1]  && !(y_int[0]==y_l && y_int[1]>y_l) ) || (x_int[0]>x_int[1] && y_int[0]>y_l && y_int[1]==y_l ))
    {
        Real x_place = x_int[0];
        Real y_place = y_int[0];
        x_int[0] = x_int[1]; x_int[1] = x_place;
        y_int[0] = y_int[1]; y_int[1] = y_place;
    }
    //set some corner points based on where intersections are: 6 cases 
    if(y_int[0] == y_l)
    {
        //lower left triangle
        if(y_l<y_int[1] && y_int[1]<y_l+dx && x_int[1]==x_l)
        {
            x.push_back(x_l);y.push_back(y_l);
        }
        //lower right triangle
        else if(y_l<y_int[1] && y_int[1]<y_l+dx && x_int[1]==x_l+dx)
        {
            x.push_back(x_l+dx);y.push_back(y_l+dx);
            x.push_back(x_l);y.push_back(y_l+dx);
            x.push_back(x_l);y.push_back(y_l);
        }
        //vertical trapezoid
        else
        {
            x.push_back(x_l);y.push_back(y_l+dx);
            x.push_back(x_l);y.push_back(y_l);
        }
    }
    //upper left triangle
    else if (y_int[0] == y_l +dx)
    {
        x.push_back(x_l);y.push_back(y_l);
        x.push_back(x_l+dx);y.push_back(y_l);
        x.push_back(x_l+dx);y.push_back(y_l+dx);
    }
    else
    {
        //horizontal trapezoid 
        if(y_l<y_int[1] && y_int[1]<y_l+dx)
        {
            x.push_back(x_l);y.push_back(y_l);
            x.push_back(x_l+dx);y.push_back(y_l);
        }
        //uppper right triangle
        else
        {
            x.push_back(x_l);y.push_back(y_l+dx);
            x.push_back(x_l);y.push_back(y_l);
            x.push_back(x_l+dx);y.push_back(y_l);
        }
    }
    
}

void
writeLevelIVSFAB(const LevelData<IVSFAB<Real>>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  LevelData<FArrayBox> LD(a_dataPtr->disjointBoxLayout(),a_dataPtr->nComp(),a_dataPtr->ghostVect());
  for(DataIterator dit = LD.dataIterator();dit.ok();++dit)
  {
    const IVSFAB<Real>& currIVSF = (*a_dataPtr)[dit()];
    const IntVectSet& ivs = currIVSF.getIVS();
    const Box& bx = ivs.minBox();
    int ncomp = a_dataPtr->nComp();

    FArrayBox& fab = LD[dit()];
    fab.setVal(sqrt(-1.0));//or whatever other null value
    for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit)
        {
        IntVect iv = ivsit();
        for (int ivar = 0; ivar < ncomp; ivar++)
            {
            fab(iv, ivar) = currIVSF(iv, ivar);
            }
        }
  }

  writeLevel(&LD);
}

void print_matrix(LAPACKMatrix& A,int digits)
{
    std::pair<int,int> dims(A.dims());
    std::cout<<"[";
    for(int i=0;i<dims.first;i++)
    {
        for(int j = 0;j<dims.second;j++)
        {
            if (A(i,j) ==0)
            {
                std::cout<<0<<" ";
            }
            else
            {
                //std::cout<<A(i,j)<<" ";
                std::cout<<setprecision(digits)<<scientific<<A(i,j)<<" ";
            }
        }
        std::cout<<"; ... \n";
    }
    std::cout<<"];\n";
}

void writeIVSFAB(const IVSFAB<Real>* a_dataPtr)
{
    const IVSFAB<Real>& currIVSF = (*a_dataPtr);
    const IntVectSet& ivs = currIVSF.getIVS();
    const Box& bx = ivs.minBox();
    int ncomp = a_dataPtr->nComp();

    FArrayBox fab(bx,ncomp);
    fab.setVal( sqrt(-1.0) );//or whatever other null value
    for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit)
    {
        IntVect iv = ivsit();
        for (int ivar = 0; ivar < ncomp; ivar++)
            {
            fab(iv, ivar) = currIVSF(iv, ivar);
            }
    }
    writeFAB(&fab);
}

void writeIVS(const IntVectSet* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  Vector<Box> boxes = a_dataPtr->boxes();
  Vector<int> procs(boxes.size(), 0);
  DisjointBoxLayout dbl(boxes, procs);
  writeDBL(&dbl);
}

void shiftMoments(Vector<Real>& moments,int Q, Real sx, Real sy)
{
    CH_TIME("shiftMoments");
    Vector<Real> momentsCopy = moments; //need a copy so we can use old moments
    int compOuter, compInner;
    Real m;
    
    for(int px = 0; px<=Q; px++) //(px,py) is the moment that we want to shift
    {
        for(int py = 0; py <= Q-px; py++)
        {
            compOuter = px*Q + (3*px - px*px)/2 + py;
            m = 0;
            //now we need to have an inner loop that does the combinatorics of the shift
            for(int qx = 0;qx <= px; qx++)
            {
                for(int qy = 0;qy <= py; qy++)
                {
                    compInner = qx*Q + (3*qx - qx*qx)/2 + qy;
                    m += r8_choose(px,qx) * pow(sx,px-qx) * r8_choose(py,qy) * pow(sy,py-qy) * momentsCopy[compInner];
                }
            }

            moments[compOuter] = m; //fill in the px,py moment
        }
    }
}

void writeLDFAB(const LevelData<FArrayBox>* dataPtr)
{
    for(DataIterator dit = dataPtr->dataIterator();dit.ok();++dit)
    {
        writeFAB( &((*dataPtr)[dit()]) );
    }
}

//integral over regular face
Real faceMom(int px,int py, int d, int s)
{
    Real m = 0; Real sgn;
    sgn = s==0 ? -1:1; //for hi and lo side of CELL
    
    if(px < 0 || py <0)
    {
        return m;
    }

    if(d == 0)
    {
        m = pow(.5*sgn, px) * (1.0/(py+1.0)) * (pow(.5,py+1) - pow(-.5,py+1) );
    }
    if(d==1)
    {
        m = pow(.5*sgn, py) * (1.0/(px+1.0)) * (pow(.5,px+1) - pow(-.5,px+1) );
    }

    return m;
}

Real faceMom(int px,int py, int d, int s, Vector<Real> bds)
{
    Real m = 0; Real sgn;
    sgn = s==0 ? -1:1; //for hi and lo side of CELL

    if(px < 0 || py <0)
    {
        return m;
    }

    if(d == 0)
    {
        m = pow(.5*sgn, px) * (1.0/(py+1.0)) * (pow(bds[1],py+1) - pow(bds[0],py+1) );
    }
    if(d==1)
    {
        m = pow(.5*sgn, py) * (1.0/(px+1.0)) * (pow(bds[1],px+1) - pow(bds[0],px+1) );
    }

    return m;
}

void writePETSC(Vec u, const LevelData<FArrayBox>& cellIDs, LevelData<FArrayBox>& cellU)
{
    const PetscScalar *avec;
    VecGetArrayRead(u,&avec); 
    IntVect iv;
    for(DataIterator dit(cellU.dataIterator());dit.ok();++dit)
    {
        const FArrayBox& cellIDFab = cellIDs[dit()];
        FArrayBox& uFab = cellU[dit()];

        Box region = (cellU.disjointBoxLayout() )[dit];
        int currID;

        for(BoxIterator bit(region);bit.ok();++bit)
        {
            iv = bit();
            
                for(int i=0;i<2;i++)
                {
                    currID = cellIDFab(iv,i);
                    if(currID != -1)
                    {
                        uFab(iv,i) = avec[currID];
                    }
                    
                }
        }
        int asdg;
        asdg=0;
    }

    VecRestoreArrayRead(u,&avec);
}

int pseudoInvert(LAPACKMatrix& A)
{
    CH_TIME("LAPACKMatrix::pseudoInvertUsingSVD");
  
  //int retval = invertUsingSVD(a_maxiter, a_tol);
  int M = A.dims().first;
  int N = A.dims().second;


  LAPACKMatrix B(M, M);
  B.setToIdentity();

  int NRHS = B.dims().second;
  int LDA = M;
  int LDB = Max(M,N);

  const int maxMxN = 40000;
  const int minMN = 100;
  int LWORK[2] = {1,1};
  LWORK[0] = maxMxN;
  LAPACKMatrix WORK(1, maxMxN);
  WORK.setVal(0);
  Real S[minMN];
  int IWORK[250*minMN];

  Real RCOND = -1;
  int INFO;
  int RANK;
  
  
  LAPACKMatrix Bcopy = B;

  LAPACK(GELSD,gelsd)(&M, &N, &NRHS, A.dataPtr(), &LDA, 
            B.dataPtr(), &LDB, S, &RCOND, &RANK,
            WORK.dataPtr(), LWORK, IWORK, &INFO);

  if(INFO != 0)
  {
    CH_assert(false);
  }


  A.define(N,M);

  for(int i=0;i<N;i++)
  {
    for(int j=0;j<M;j++)
    {
      A(i,j) = B(i,j);
    }
  }
  

  return INFO;
}


#include "NamespaceFooter.H"
