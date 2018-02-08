#include <jem/base/Array.h>
#include <jem/base/array/select.h>
#include <jem/base/array/utilities.h>
#include <jem/base/PrecheckException.h>
#include <jem/base/Error.h>
#include <jem/numeric/utilities.h>

#include <jem/base/System.h>


#include "utilities.h"
#include "Array.h"

extern "C"
{
  #include  <math.h>
}

using jem::ALL;
//using jem::numeric::abs;

using namespace jem;


//-----------------------------------------------------------------------
//   constants
//-----------------------------------------------------------------------


const idx_t  STRAIN_COUNTS[4] = { 0, 1, 3, 6 };

const double PI               = 3.1415926535897932;
const double RAD_120          = 0.6666666666666667 * PI;
const double RAD_240          = 2.0 * RAD_120;

const double ONE_THIRD        = 0.3333333333333333;

const double EPS              = 1.e-16;

// --------------------------------------------------------------------
//  getShapeFuncsFunc
// --------------------------------------------------------------------
// A general function to compute the shape function values. This function
// will call a dimensional specific function.

ShapeFuncsFunc   getShapeFuncsFunc

  ( idx_t               rank )

{
  JEM_ASSERT ( rank >= 1 && rank <= 3 );


  if      ( rank == 1 )
  {
    return & get1DShapeFuncs;
  }
  else if ( rank == 2 )
  {
    return & get2DShapeFuncs;
  }
  else
  {
    return & get3DShapeFuncs;
  }
}


// get1DShapeFuncs -------------------------------------------------------------

void                  get1DShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n )
{
  sfuncs(0,0) = n[0];
}


// get2DShapeFuncs -------------------------------------------------------------

void                  get2DShapeFuncs

  ( const Matrix&       s,
    const Vector&       n )
{
  JEM_ASSERT ( s.size(0) == 2 &&
               s.size(1) == 2 * n.size() );

  const idx_t  nodeCount = n.size ();

  s = 0.0;

  for ( idx_t inode = 0; inode < nodeCount; inode++ )
  {
    idx_t  i = 2 * inode;

    s(0,i + 0) = n[inode];
    s(1,i + 1) = n[inode];
  }
}


// get3DShapeFuncs -------------------------------------------------------------

void                  get3DShapeFuncs

  ( const Matrix&       s,
    const Vector&       n )
{
  JEM_ASSERT ( s.size(0) == 3 &&
                 s.size(1) == 3 * n.size() );

  const idx_t  nodeCount = n.size ();

  s = 0.0;

  for ( idx_t inode = 0; inode < nodeCount; inode++ )
  {
    idx_t  i = 3 * inode;

    s(0,i + 0) = n[inode];
    s(1,i + 1) = n[inode];
    s(2,i + 2) = n[inode];
  }
}


//-----------------------------------------------------------------------
//   getShapeGradsFunc
//-----------------------------------------------------------------------
// A general function to compute the shape function derivative values. 
// This function will call a dimensional specific function.

ShapeGradsFunc getShapeGradsFunc ( idx_t rank )
{
  JEM_ASSERT ( rank >= 1 && rank <= 3 );


  if      ( rank == 1 )
  {
    return & get1DShapeGrads;
  }
  else if ( rank == 2 )
  {
    return & get2DShapeGrads;
  }
  else
  {
    return & get3DShapeGrads;
  }
}


// get1DShapeGrads -------------------------------------------------------------

void              get1DShapeGrads

  ( const Matrix&   b,
    const Matrix&   g )

{
  JEM_ASSERT ( b.size(0) == 1 &&
                 g.size(0) == 1 &&
                 b.size(1) == g.size(1) );

  b = g;
}


// get2DShapeGrads -------------------------------------------------------------

void              get2DShapeGrads

  ( const Matrix&   b,
    const Matrix&   g )

{
  JEM_ASSERT ( (b.size(0) == 3 || b.size(0) == 4) &&
                g.size(0) == 2 &&
                b.size(1) == 2 * g.size(1) );

  const idx_t  nodeCount = g.size (1);


  b = 0.0;

  for ( idx_t inode = 0; inode < nodeCount; inode++ )
  {
    idx_t  i = 2 * inode;

    b(0,i + 0) = g(0,inode);
    b(1,i + 1) = g(1,inode);

    b(2,i + 0) = g(1,inode);
    b(2,i + 1) = g(0,inode);
  }
}


// get3DShapeGrads -------------------------------------------------------------

void              get3DShapeGrads

  ( const Matrix&   b,
    const Matrix&   g )

{
  JEM_ASSERT ( b.size(0) == 6 &&
               g.size(0) == 3 &&
               b.size(1) == 3 * g.size(1) );

  const idx_t  nodeCount = g.size (1);


  b = 0.0;

  for ( idx_t inode = 0; inode < nodeCount; inode++ )
  {
    idx_t  i = 3 * inode;

    b(0,i + 0) = g(0,inode);
    b(1,i + 1) = g(1,inode);
    b(2,i + 2) = g(2,inode);

    b(3,i + 0) = g(1,inode);
    b(3,i + 1) = g(0,inode);

    b(4,i + 1) = g(2,inode);
    b(4,i + 2) = g(1,inode);

    b(5,i + 2) = g(0,inode);
    b(5,i + 0) = g(2,inode);
  }
}



//-----------------------------------------------------------------------
//   Other Functions
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
//   addAxiShapeGrads
//-----------------------------------------------------------------------
void              addAxiShapeGrads

  ( const Matrix&   b,
    const Vector&   n,
    const double&   r )

{
  JEM_ASSERT ( (b.size(0) == 4) );

  const idx_t  nodeCount = b.size(1) / 2;

  for ( idx_t inode = 0; inode < nodeCount; inode++ )
  {
    idx_t  i = 2 * inode;
    
    b(3,i + 0) = n[inode] / r;
  }
}

//-----------------------------------------------------------------------
//   Cauchy2DVec2Mat
//-----------------------------------------------------------------------

void                  Cauchy2DVec2Mat

  ( const Matrix&       CMat,
    const Vector&       CVec )
    
{
  JEM_ASSERT ( CMat.size(0) == 4 &&
               CMat.size(1) == 4 &&
               CVec.size()  == 3 );
  
  CMat = 0.0;
  CMat(0,0) = CMat(2,2) = CVec[0];
  CMat(1,1) = CMat(3,3) = CVec[1];
  CMat(1,0) = CMat(0,1) = CVec[2];
  CMat(2,3) = CMat(3,2) = CVec[2];
}

//-----------------------------------------------------------------------
//   solveLinearEqua
//-----------------------------------------------------------------------

void                  solveLinearEqua

  ( double&       ans,
    const double  a,
    const double  b )
{
  ans = -b/a;
}

//-----------------------------------------------------------------------
//   solveQuadEqua
//-----------------------------------------------------------------------

void                  solveQuadEqua

  ( Vector&       ans,
    const double  a,
    const double  b,
    const double  c )
{
  if ( std::abs(a) < EPS )
  {
    ans.resize( 1 );
    ans[0] = - c / b;
  }
  else
  {
    double d = b * b - 4.0 * a * c;
    
    if ( d < 0 )
    {
      using namespace jem;
      throw Error (
      JEM_FUNC,
      " imaginary principal values !!! "
      );  
    }
    else
    {
      ans.resize( 2 );
      ans[0] = 0.5 * (-b + sqrt( d )) / a;
      ans[1] = 0.5 * (-b - sqrt( d )) / a;
    }
  }
}

//-----------------------------------------------------------------------
//   solveCubicEqua
//-----------------------------------------------------------------------

idx_t                 solveCubicEqua

  ( Tuple<double,3>& ans, 
    const double     a,
    const double     b,
    const double     c,
    const double     d)

{  
  ans = NAN;

  if ( std::abs(a) < EPS )
  {
    if ( std::abs(b) < EPS )
    {
      ans[0] = - d / c;

      return 1;
    }
    else
    {
      double D = c * c - 4.0 * b * d;
      
      if ( D < 0.0 )
      {        
        return 0;
      }
      else
      {
        D      = sqrt(D);

        ans[0] = 0.5 * (-c + D) / b;
        ans[1] = 0.5 * (-c - D) / b;
       
        return 2;
      }
    }
  }
  else
  {
    double kk  = 1.0 / a;
    
    double aa  = b * kk;
    double bb  = c * kk;
    double cc  = d * kk;

    double aa3 = aa / 3.0;

    double p, q, r;
    double phi, help;
   
    q  = ( aa * aa - 3.0 * bb ) / 9.0;
    r  = ( 2.0 * aa * aa * aa - 9.0 * aa * bb + 27.0 * cc ) / 54.0;
 
    help = r / sqrt( q * q * q );
   
    if ( std::abs (help) > 1.0 )
    {
      help = ( help < 0 ) ? -1.0 : 1.0; // prevent rounding errors
    }
   
    phi = acos ( help );
    p   = sqrt ( q    );

    ans[0] = -2.0 * p * cos ( phi / 3.0 )           - aa3;
    ans[1] = -2.0 * p * cos ( phi / 3.0 + RAD_120 ) - aa3;
    ans[2] = -2.0 * p * cos ( phi / 3.0 - RAD_120 ) - aa3;

    return 3;
  }
}
  
//-----------------------------------------------------------------------
//   invert2x2
//-----------------------------------------------------------------------

void  invert2x2

  ( const Matrix& m )

{
  // no check for size and singularity!

  double d = m(0,0) * m(1,1) - m(0,1) * m(1,0);
  double s = 1. / d;

  double t = m(0,0);

  m(0,0) =  m(1,1) * s;
  m(0,1) = -m(0,1) * s;
  m(1,0) = -m(1,0) * s;
  m(1,1) =  t      * s;
}


//-----------------------------------------------------------------------
//   evalMcAuley
//-----------------------------------------------------------------------

double                evalMcAuley

  ( const double x )
{
  return (x > 0 ? x : 0.0);
}

//-----------------------------------------------------------------------
//    evalHeaviside
//-----------------------------------------------------------------------

double               evalHeaviside

  ( const double x )
{
  return (x < 0 ? -1.0 : 1.0);
} 

//-----------------------------------------------------------------------
//    getUnique
//-----------------------------------------------------------------------

void               getUnique

  (       Vector& U,
    const Vector& V )
{
  IdxVector V_double (V.size());
  V_double = 1;
  U.resize(0);

  for ( int i=0; i<V.size(); i++ )
  {
    if (V_double[i] == 0)
    {
      continue;
    }
    
    for ( int j=i+1; j<V.size(); j++ )
    {
      (V[i] == V[j]) ? V_double[j] = 0: 0;
    }
    
    if ( V_double[i] == 1 )
    {
      U.reshape(U.size()+1);
      U[U.size()-1] = V[i];
    }
  }

} 

void            getUnique

  (       IdxVector& U,
    const IdxVector& V )
{
  IdxVector V_double (V.size());
  V_double = 1;
  U.resize(0);

  for ( int i=0; i<V.size(); i++ )
  {
    if (V_double[i] == 0)
    {
      continue;
    }
    
    for ( int j=i+1; j<V.size(); j++ )
    {
      (V[i] == V[j]) ? V_double[j] = 0: 0;
    }
    
    if ( V_double[i] == 1 )
    {
      U.reshape(U.size()+1);
      U[U.size()-1] = V[i];
    }
  }

} 

//-----------------------------------------------------------------------
//    getThisValue
//----------------------------------------------------------------------

void            getThisValue

  (       IdxVector& pos,
    const Vector&    V,
    const double    d )
{
  pos.resize(0);
  for ( idx_t i=0; i<V.size(); i++ )
  {
    if ( V[i] == d )
    {
      pos.reshape(pos.size()+1);
      pos[pos.size()-1]=i;
    }
  }

}

void            getThisValue

  (       IdxVector& pos,
    const IdxVector& V,
    const double    d )
{
  pos.resize(0);
  for ( idx_t i=0; i<V.size(); i++ )
  {
    if ( V[i] == d )
    {
      pos.reshape(pos.size()+1);
      pos[pos.size()-1]=i;
    }
  }

}

void            getThisValue

  (       int&       pos,
    const Vector&    V,
    const double     d )
{
  
  for ( int i=0; i<V.size(); i++ )
  {
    if ( V[i] == d )
    {
      pos=i;
    }
  }

}

void            getThisValue

  (       int&       pos,
    const IdxVector& V,
    const double     d )
{
  
  for ( int i=0; i<V.size(); i++ )
  {
    if ( V[i] == d )
    {
      pos=i;
    }
  }

}

//-----------------------------------------------------------------------
//   makeNorm
//-----------------------------------------------------------------------

void            makeNorm

    (       Vector&             norm   )
    
{
  norm /= getLength(norm);
}

//-----------------------------------------------------------------------
//   LVec_
//-----------------------------------------------------------------------

double          getLength

    ( const Vector&             norm   )
    
{
  double l = 0.0;
  for ( int i=0; i<norm.size(); i++ )
  {
    l += norm[i]*norm[i];
  }
  
  return sqrt(l);  
}

//-----------------------------------------------------------------------
//   signERIK
//-----------------------------------------------------------------------

double                  signERIK

  ( const double     d )
  
{
  return (d>=0) ? 1.0: -1.0;
}

