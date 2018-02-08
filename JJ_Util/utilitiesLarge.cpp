#include <jem/base/Array.h>
#include <jem/base/array/select.h>
#include <jem/base/array/utilities.h>
#include <jem/base/PrecheckException.h>
#include <jem/base/Error.h>

#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/LUSolver.h>

#include <jem/numeric/utilities.h>

#include <jem/base/System.h>


#include "utilitiesLarge.h"
#include "Array.h"

extern "C"
{
  #include  <math.h>
}

using jem::ALL;
//using jem::numeric::abs;
using jem::numeric::inverse;
using jem::numeric::matmul;

using namespace jem;


//-----------------------------------------------------------------------
//   constants
//-----------------------------------------------------------------------


void               getDeformationGradient
  
  ( const Cubix&   defgrad,
    const Matrix&  coords,
    const Vector&  disp,
    const Cubix&   grads )
{
  // Compute the deformation gradient
  //
  // Input:
  // coord          - nodal coordinates in initial configuration
  // disp           - current nodal displacement values
  // grads          - derivatives of shape function in initial configuration
  //
  // Output:
  // defgrad        - deformation gradient matrix for every integration point
  
  int ipcount=grads.size(2);
  int nodecount=grads.size(1);
  defgrad=0.0;
  defgrad(2,2,ALL)=1.0;
  for (int ip=0;ip<ipcount;ip++)
     {
      for (int nd=0;nd<nodecount;nd++)
         { 
           defgrad(0,0,ip)+=grads(0,nd,ip)*(coords(0,nd)+disp[2*nd]);  //F11
           defgrad(0,1,ip)+=grads(1,nd,ip)*(coords(0,nd)+disp[2*nd]);  //F12
           defgrad(1,0,ip)+=grads(0,nd,ip)*(coords(1,nd)+disp[2*nd+1]);//F21
           defgrad(1,1,ip)+=grads(1,nd,ip)*(coords(1,nd)+disp[2*nd+1]);//F22
         }
      }  
}

// ---------------------------------------------------------------
//    getCurrentgrads  // gives derivations w.r.t the current configuration
// ---------------------------------------------------------------

void             getCurrentgrads

  ( const Cubix&  deformedgrads,
    const Cubix&  grads,
    const Cubix&  defgrad )
    
{
  // Compute derivatives of shape function in current configuration
  //
  // Input:
  // F               - deformation gradient
  // Undeformedgrads - derivatives of shape function in initial configuration
  //
  // Output:
  // deformedgrads   - current nodal displacement values
  
  int ipcount=grads.size(2);
  
  Matrix Finv(3,3),Q(2,2);
  
  for (int ip=0;ip<ipcount;ip++)
  {
    Finv=inverse(defgrad(ALL,ALL,ip));
    Q(0,0)=Finv(0,0);Q(0,1)=Finv(1,0);
    Q(1,0)=Finv(0,1);Q(1,1)=Finv(1,1);
  
    deformedgrads(ALL,ALL,ip) = matmul(Q,grads(ALL,ALL,ip));
  }

}

//-----------------------------------------------------------------------
//   getNLShapeGradsFunc
//-----------------------------------------------------------------------


NLShapeGradsFunc getNLShapeGradsFunc ( idx_t rank )
{
  JEM_ASSERT ( rank >= 1 && rank <= 3 );


  if ( rank == 2 )
  {
    return & get2DNLShapeGrads;
  }
  else
  {
    using namespace jem;
    throw Error (
    JEM_FUNC,
    " Updated Lagrangian shape function gradients only defined for 2D!!! "
    );  
  }
}

//-----------------------------------------------------------------------
//   get2DNLShapeGrads
//-----------------------------------------------------------------------


void              get2DNLShapeGrads

  ( const Matrix&   bNL,
    const Matrix&   g )

{
  JEM_ASSERT ( bNL.size(0) == 4 &&
               g.size(0) == 2 &&
               bNL.size(1) == 2 * g.size(1) );
                 
  const idx_t  nodeCount = g.size (1);

  // construct linear part of shape function
  bNL = 0.0;

  for ( idx_t inode = 0; inode < nodeCount; inode++ )
  {
    idx_t  i = 2 * inode;
   
    // non-linear part
    bNL(0,i + 0) = g(0,inode);
    bNL(1,i + 0) = g(1,inode);
    bNL(2,i + 1) = g(0,inode);
    bNL(3,i + 1) = g(1,inode);
  }

  
}

//-----------------------------------------------------------------------
//   getStressTensor
//-----------------------------------------------------------------------

void                getStressTensor

  ( const Matrix&       S,
    const Vector&  stress )

{
  JEM_ASSERT ( S.size(0) == 3 &&
               S.size(1) == 3 &&
               stress.size() == 3 );
   
   
   S=0.0;
   S(0,0)=stress[0];
   S(1,1)=stress[1];
   S(0,1)=S(1,0)=stress[2];
}

//-----------------------------------------------------------------------
//   getStressMatrix
//-----------------------------------------------------------------------

void                getStressMatrix

  ( const Matrix&       S,
    const Vector&  stress )

{
  JEM_ASSERT ( S.size(0) == 4 &&
               S.size(1) == 4 &&
               stress.size() == 3 );
   
   
   S=0.0;
   S(0,0)=S(2,2)=stress[0];
   S(1,1)=S(3,3)=stress[1];
   S(0,1)=S(1,0)=S(2,3)=S(3,2)=stress[2];
}

//-----------------------------------------------------------------------
//   getStressVector
//-----------------------------------------------------------------------

void                getStressVector

  ( const Vector&  stress,
    const Matrix&       S )

{
  JEM_ASSERT ( S.size(0) == 3 &&
               S.size(1) == 3 &&
               stress.size() == 3 );
   
   stress[0] = S(0,0);
   stress[1] = S(1,1);
   stress[2] = S(1,0);
}




