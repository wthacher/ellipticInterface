#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#include "EBReconstruction.H"

#include <map>


#include "NamespaceHeader.H"


void initializePsiNodes(ImplicitFunction& psi,
                    Real xl, Real yl, Real dx,
                    LevelData<NodeFArrayBox>& psiNodes,
                    LayoutData<IntVectSet>& EBCells,
                    Vector<LayoutData<IntVectSet>* > faceIntersections)
{
    CH_TIME("EBReconstruction::initializeBField");
    for(DataIterator dit = psiNodes.dataIterator(); dit.ok(); ++dit)
    {
        FArrayBox& psiNodesFab = psiNodes[dit()].getFab();
        IntVect iv;
        Real x,y;
        
        //PASS 1: tag the potential EB cells
        for(BoxIterator bit(psiNodesFab.box() );bit.ok();++bit)
        {
            iv = bit();
            x = iv[0] * dx + xl;
            y = iv[1] * dx + yl;

            psiNodesFab(iv) = psi.val(x,y);
        }
    }

    tagEBCells(psiNodes,EBCells,faceIntersections);

}

void tagEBCells(const LevelData<NodeFArrayBox>& psiNodes,
                LayoutData<IntVectSet>& EBCells,
                Vector<LayoutData<IntVectSet>* > faceIntersections)
{
    CH_TIME("EBReconstruction::tagEBCells");
    for(DataIterator dit = psiNodes.dataIterator();dit.ok();++dit)
    {
        const FArrayBox& bNodesFab = psiNodes[dit()].getFab();
        IntVectSet& EBCellsFab = EBCells[dit()];
        IntVect iv;
        //tag EB cells 
        for (int d = 0;d<2;d++)
        {
            LayoutData<IntVectSet>& faceLD = *faceIntersections[d];
            IntVectSet& taggedFaces = faceLD[dit()];
            IntVect tangShift = IntVect::Unit;
            IntVect normShift = IntVect::Zero;
            tangShift[d] = 0; normShift[d]=1;
            Box validNodes = bNodesFab.box();
            validNodes.growHi((d+1)%2,-1);
            for(BoxIterator bit(validNodes);bit.ok();++bit)
            {
                iv = bit();
                if((bNodesFab(iv))* (bNodesFab(iv+tangShift)) <0)
                {
                    taggedFaces |= iv;
                    EBCellsFab |= iv;
                    EBCellsFab |= iv - normShift;
                }
            }
        }
        Box validCells = bNodesFab.box(); validCells.enclosedCells();
        EBCellsFab &= validCells;
        int asd;
        asd+=1;

    }
 
}

void findIntersections(ImplicitFunction& psi,
                Real xl, Real yl, Real dx,
                const LevelData<NodeFArrayBox>& psiNodes,
                LayoutData<IntVectSet>& EBCells,
                Vector<LayoutData<IntVectSet>* > faceIntersections,
                LevelData<FluxBox>& faceIntVals)
{
    CH_TIME("EBReconstruction::findIntersections");
    for(DataIterator dit = psiNodes.dataIterator();dit.ok();++dit)
    {
        FluxBox& faceIntValsFab = faceIntVals[dit()];
        faceIntValsFab.setVal(-1);
 
        IntVectSet& EBCellsFab = EBCells[dit()];
        IntVect iv;

        Real xint,x1,y1,x2,y2;

        Box region = (psiNodes.disjointBoxLayout() )[dit()];
        for (int d = 0;d<2;d++)
        {
            LayoutData<IntVectSet>& faceLD = *faceIntersections[d];
            IntVectSet& taggedFaces = faceLD[dit()];
            FArrayBox& faceIntValsD = faceIntValsFab[d];
            IntVect tangShift = IntVect::Unit;
            IntVect normShift = IntVect::Zero;
            tangShift[d] = 0; normShift[d]=1;

            Box validNodes = region; validNodes.grow(d,1);
            IntVectSet validTaggedFaces = taggedFaces;
            validTaggedFaces &= validNodes;
            for(IVSIterator ivIt(validTaggedFaces);ivIt.ok();++ivIt)
            {
                iv = ivIt();

                x1 = iv[0]*dx + xl;
                y1 = iv[1]*dx + yl;

                if(d==0)
                {
                    x2 = x1;
                    y2 = y1+dx;
                }
                else
                {
                    x2 = x1+dx;
                    y2=y1;
                }

                psi.secantMethod(x1,y1,x2,y2);
 
                if(d==0)
                {
                    xint = (y2-y1)/dx;
                }
                else
                {
                    xint = (x2-x1)/dx;
                }

                if(xint <0 || xint >1 || isnan(xint))
                {
                    printf("whiff at: %d, %d \n" ,iv[0], iv[1]);//CH_assert(false);
                }
                else
                {
                    faceIntValsD(iv) = xint;
                }

            }
            int asd;
            asd+=1;

        }        

    }

     faceIntVals.exchange();   
}


void getMoments(ImplicitFunction& psi,
                Real xl, Real yl, Real dx,
                const LevelData<NodeFArrayBox>& psiNodes,
                const LevelData<FluxBox>& faceIntVals,
                LevelData<IVSFAB<Real>>& moments,
                int Q)
{
    CH_TIME("EBReconstruction::getMoments");
    for(DataIterator dit(faceIntVals.dataIterator());dit.ok();++dit)
    {
        IFPolygon IFP;
        romberg<IFPolygon> rmp;
        Vector<Real> xB, yB, xc, yc;
        const FluxBox& faceIntValsFab = faceIntVals[dit()];
        IVSFAB<Real>& momentsFab = moments[dit()];
        IntVectSet validEBCells = momentsFab.getIVS();
        validEBCells &= faceIntVals.disjointBoxLayout()[dit()];
        const FArrayBox& psiNodesFab = psiNodes[dit()].getFab();
        IntVect iv,nIv;
        int comp;
        //int Q = 2;
        int n = ((Q+1)*(Q+2))/2;
        comp = 0;
        Vector<Real> fullMoments(n,0);
        Real xCtd, yCtd;

        for(int qx = 0; qx<=Q; qx++)
        {
            for(int qy = 0; qy<= Q- qx; qy++)
            {
                fullMoments[comp] = ( (pow(.5,qx+1) - pow(-.5,qx+1) ) / (qx+1) )*((pow(.5,qy+1) - pow(-.5,qy+1) ) / (qy+1));
                comp +=1;
            }
        }
        
        for(IVSIterator ivIt(validEBCells);ivIt.ok();++ivIt)
        {
            iv = ivIt();
            xB.clear();
            yB.clear();
            xc.clear(); xc.push_back(-.5);
            yc.clear(); yc.push_back(-.5);
            //find intersections on this cell
            for(int d=0;d<2;d++)
            {
                for(int s=0;s<2;s++)
                {
                    nIv = iv;
                    nIv[d] = iv[d] +s;
                    if(faceIntValsFab[d](nIv) != -1)
                    {
                        if(d==0)
                        {
                            xB.push_back(.5*(s==0 ? -1:1));
                            yB.push_back(faceIntValsFab[d](nIv) - .5);
                        }
                        else
                        {
                            yB.push_back(.5*(s==0 ? -1:1));
                            xB.push_back(faceIntValsFab[d](nIv) - .5);
                        }
                        
                    }
                }
            }

            xCtd = iv[0]*dx + xl + .5*dx;
            yCtd = iv[1]*dx + yl + .5*dx;


            get_points(xc,yc,xB,yB,1);
            IFP.define(xB,yB,xc,yc, xCtd, yCtd, dx, psi);
            
            comp=0;
            for(int qx = 0; qx<=Q; qx++)
            {
                for(int qy = 0; qy<= Q- qx; qy++)
                {
                    momentsFab(iv,comp) = rmp.integrate(IFP,qx,qy,1e-15,12);
                    if(comp == 0)
                    {
                        if(momentsFab(iv,comp) >1 || momentsFab(iv,comp)<0)
                        {
                            pout() <<"bad";
                        }
                    }
                    momentsFab(iv,n+comp) = fullMoments[comp] - momentsFab(iv,comp); //subtract lo part
                    momentsFab(iv,2*n+comp) = rmp.integrateArea(IFP,qx,qy,2,1e-15,12);
                    momentsFab(iv,3*n+comp) = rmp.integrateArea(IFP,qx,qy,0,1e-15,12);
                    momentsFab(iv,4*n+comp) = rmp.integrateArea(IFP,qx,qy,1,1e-15,12);
                    comp +=1;
                }
            }
            

        }
        
    }
    moments.exchange();
}


//output functions
void tagCellsAndFaces(const LevelData<NodeFArrayBox>& psiNodes,
                const LayoutData<IntVectSet>& EBCells,
                RefCountedPtr<LayoutData<IntVectSet> >& taggedFluxCellsLD,
                Vector< RefCountedPtr< LayoutData<IntVectSet> > >& taggedFacesLD,
                int order,
                int velOrder)
{
    CH_TIME("EBReconstruction::tagCellsAndFaces");
    for(DataIterator dit = psiNodes.dataIterator(); dit.ok(); ++dit)
    {
        const IntVectSet& EB = (EBCells[dit()]);
        const FArrayBox& psiNodesFab = psiNodes[dit()].getFab();
        IntVectSet& taggedFluxCells = (*taggedFluxCellsLD)[dit()];

        Box currBox = psiNodesFab.box();
        currBox.grow(-1*psiNodes.ghostVect()); //you only want faces on full cells
        Box dirBox,stencil; IntVect iv, normalShift, tangShift;
        BoxIterator bit, stencilIt;
        int Q = order;

        int regRadius = (Q+1)/2;
        int regRadiusVel = (velOrder+1)/2;
        
        //tag faces and cells
        for(int dir = 0; dir < CH_SPACEDIM; dir++)
        {
            IntVectSet& taggedFaces = (*(taggedFacesLD[dir]) )[dit()];
            tangShift[dir] = 0; tangShift[(dir+1)%2] = 1;
            normalShift[dir] = 1; normalShift[(dir+1)%2] = 0;
            dirBox = currBox;
            dirBox.growHi((dir+1)%2,-1); //grow the top corner so we dont step off
            bit.define(dirBox);
            for(bit.begin(); bit.ok(); ++bit)
            {
                iv = bit();
                stencil.define(iv-regRadiusVel*(normalShift+tangShift), iv + (regRadiusVel-1)*normalShift + regRadiusVel*tangShift);
                stencilIt.define(stencil);
                for(stencilIt.begin(); stencilIt.ok(); ++stencilIt)
                {
                    if(EB.contains(stencilIt()) ) //if anything in your stencil is cut then dip out
                    {
                        taggedFaces |=iv;
                        taggedFluxCells |= iv;
                        taggedFluxCells |= iv-normalShift;
                        break; //you can just leave, we only need one
                    }
                }  
            }
        
        }
        //iterate through EBCells and tag adjacent ones for rhs reconstruction
        IntVect right(1,0); IntVect up(0,1);
        
        Box validCells(psiNodesFab.box()); validCells.enclosedCells();
        taggedFluxCells&=validCells;
        taggedFluxCells-=EB;//we dont want any EB cells in here

    }
}

#include "NamespaceFooter.H"