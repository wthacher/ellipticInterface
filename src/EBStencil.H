#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef __EBSTENCIL_H_
#define __EBSTENCIL_H_

#include "IntersectionUtils.H"

#include "NamespaceHeader.H" 

/*This class wraps an ivsfab<real> which represents a stencil. It is used in 2 ways:
1) as a stencil for gradients on full and partial faces
2) as a stencil for flux divergence on EB and tagged cells

The 'add' functions modify the ivs, but not the data
The 'add_' functions where _ is some other data structure adds that type of data to the stencil
The 'add_ToFlux' functions adds VTO components using gradients to the flux divergence
 */
class EBStencil
{
public:
    IVSFAB<Real> stencilIVSF; //all the stencil data is in an underlying IVSFAB<Real>
    Real jumpConstant; //this keeps track of the excess from the JCs
    Real rhsConstant; //keeps track of excess from rhs
    //need default constructor so we can put in BoxLayoutData
    EBStencil()
    {
        IntVectSet empt;
        stencilIVSF.define(empt,2);
        jumpConstant = 0;
        rhsConstant = 0;
    }

    void define(int nComp)
    {
        IntVectSet empt;
        stencilIVSF.define(empt,nComp);
        jumpConstant = 0;
    }

    void defineBox(int nComp, const Box& b)
    {
        IntVectSet bIvs(b);
        stencilIVSF.define(bIvs, nComp);
        jumpConstant = 0;
    }

    //add the new box to IVSFAB, zero out new values
    void add(const Box& a_box)
    {
        CH_TIME("EBStencil::add(Box)");
        IVSFAB<Real> stencilIVSFCopy(stencilIVSF.getIVS(), stencilIVSF.nComp());
        Interval i(0,stencilIVSF.nComp()-1);
        Box b = stencilIVSF.getIVS().minBox();
        stencilIVSFCopy.copy(b,i,b,stencilIVSF,i);

        IntVectSet addCells;
        addCells.define(a_box);
        addCells -= stencilIVSF.getIVS(); //these are the new cells

        IntVectSet allCells = addCells;
        allCells |= stencilIVSF.getIVS();
        //redefine on new IVS
        stencilIVSF.define(allCells, stencilIVSF.nComp());

        //copy old data back
        stencilIVSF.copy(b,i,b,stencilIVSFCopy,i);

        //zero out new comps
        for(IVSIterator ivIt(addCells);ivIt.ok();++ivIt)
        {
            for(int n=0;n<stencilIVSF.nComp();n++)
            {
                stencilIVSF(ivIt(),n) = 0;
            }
        }

    };

    //add new IVS to m_ivs, zero new values
    void add(const IVSFAB<Real>& a_IVSF)
    {
        CH_TIME("EBStencil::add(IVSFAB)");
        IVSFAB<Real> stencilIVSFCopy(stencilIVSF.getIVS(), stencilIVSF.nComp());
        Interval i(0,stencilIVSF.nComp()-1);
        Box b = stencilIVSF.getIVS().minBox();
        stencilIVSFCopy.copy(b,i,b,stencilIVSF,i);

        IntVectSet addCells = a_IVSF.getIVS();
        addCells -= stencilIVSF.getIVS(); //these are the new cells

        IntVectSet allCells = addCells;
        allCells |= stencilIVSF.getIVS();
        //redefine on new IVS
        stencilIVSF.define(allCells, stencilIVSF.nComp());

        //copy old data back
        stencilIVSF.copy(b,i,b,stencilIVSFCopy,i);

        //zero out new comps
        for(IVSIterator ivIt(addCells);ivIt.ok();++ivIt)
        {
            for(int n=0;n<stencilIVSF.nComp();n++)
            {
                stencilIVSF(ivIt(),n) = 0;
            }
        }

    };

    //add a single iv, zero out if new value
    void add(const IntVect& a_iv)
    {
        CH_TIME("EBStencil::add(IntVect)");
        if( !(stencilIVSF.getIVS().contains(a_iv) )) //if this is a new member, zero it out
        {
            IVSFAB<Real> stencilIVSFCopy(stencilIVSF.getIVS(), stencilIVSF.nComp());
            Interval i(0,stencilIVSF.nComp()-1);
            Box b = stencilIVSF.getIVS().minBox();
            stencilIVSFCopy.copy(b,i,b,stencilIVSF,i);

            IntVectSet allCells;
            allCells.define( a_iv);
            allCells |= stencilIVSF.getIVS();
            //redefine on new IVS
            stencilIVSF.define(allCells, stencilIVSF.nComp());

            //copy old data back
            stencilIVSF.copy(b,i,b,stencilIVSFCopy,i);
            //zero out new member
            for(int n=0;n<stencilIVSF.nComp();n++)
            {
                stencilIVSF(a_iv,n) = 0;
            }
        }
    };

    //just for convenience
    Real& operator()(const IntVect& a_ndin,
                       const int& a_comp)
    {
        return stencilIVSF(a_ndin, a_comp);
    };

    //add this IVSFAB in, scaled by a_scale
    void addIVSFAB(const IVSFAB<Real>& a_IVSF, const Vector<Real>& a_scale)
    {
        CH_TIME("EBStencil::addIVSFAB");
        CH_assert(a_IVSF.nComp() == stencilIVSF.nComp());
        add(a_IVSF);
        for(IVSIterator ivIt(a_IVSF.getIVS()); ivIt.ok();++ivIt)
        {
            for(int n=0;n<a_IVSF.nComp();n++)
            {
                stencilIVSF(ivIt(),n) += a_scale[n]*a_IVSF(ivIt(),n);
            }
        }
    };
    
    //add an 18 component a_vect, which represents du, by iterating through box how we are supposed to
    void addVect(const Box& a_box, const Vector<Real>& a_vect, Real a_scale, bool init)
    {
        CH_TIME("EBStencil::addVect");
        //a_vect has 2*nStenCells comps
        int nStenCells = a_vect.size()/stencilIVSF.nComp();
        int r = (int( pow(nStenCells,.5)) -1 )/2;
        //CH_assert(a_vect.size() == 18);
        // if(init)
        // {
        //     add(a_box);
        // }
        
        IntVect center = a_box.smallEnd() + r*IntVect::Unit;
        IntVect iv;
        int count = 0;
        for(int l = -r;l<r+1;l++)
        {
            for(int c = -r;c<r+1;c++) 
            {
                iv[0] = center[0] + c;
                iv[1] = center[1] + l;
                for(int i=0;i<stencilIVSF.nComp();i++)
                {
                    stencilIVSF(iv,i) += a_vect[count + i*nStenCells]*a_scale;
                }
                
                count+=1;
            }
        }
        CH_assert(count ==nStenCells);
    };
    
    //add a single value at this intvect
    void addIV(const IntVect& a_iv, Real a_val, int comp,bool init)
    {
        CH_TIME("EBStencil::addIV");
        if(init)
        {
            add(a_iv);
        }
        
        stencilIVSF(a_iv,comp) += a_val;
    };

    void clear()
    {
        stencilIVSF.setVal(0);
    }

    
    EBStencil& operator=(const EBStencil& other){}
    
};

#include "NamespaceFooter.H"

#endif