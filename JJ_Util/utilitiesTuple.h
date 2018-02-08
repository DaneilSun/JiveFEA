#ifndef UTILITIESTUPLE_H
#define UTILITIESTUPLE_H

#include <jive/Array.h>
#include <jem/base/Tuple.h>

#include "Invariants.h"

using jive::Vector;
using jive::Matrix;
using jive::IntVector;
using jem::Tuple;
using jem::idx_t;

//-----------------------------------------------------------------------
//   typedefs
//-----------------------------------------------------------------------

// A pointer to a function that computes the spatial derivatives of
// the interpolation matrix. This is the so-called B-matrix.

void            m6_to_tt6

    (       Mat6& tt6,
      const Matrix&            m6  );
      
void            tt6_to_m6

    ( const Matrix&            m6,
      const Mat6& tt6 );
      
void v6_to_t6

  (       Vec6&   t6,
    const Vector& v6 );
    
Vec6 v6_to_t6

  ( const Vector& v6 );
  
void t6_to_v6

  ( const Vector& v6,
    const Vec6&   t6 ); 


   
Vec6 tmatmul

  ( const Matrix& mat,
    const Vec6&   t6 );
    
Vec6 tmatmul

  ( const Vec6&   t6,
    const Matrix& mat );
    
Vec6 tmatmul

  ( const Vec6&   t6,
    const Mat6&   tt6 );
    
Vec6 tmatmul

  ( const Mat6&   tt6,
    const Vec6&   t6 );
    
Mat6 tmatmul

  ( const Vec6&   t0,
    const Vec6&   t1  );
    
Mat6 tmatmul

  ( const Mat6&     tt6,
    const Matrix&   mat  );    
    
Mat6 tmatmul

  ( const Matrix&   mat,
    const Mat6&     tt6  ); 
    
Mat6 tmatmul

  ( const Mat6&     tt6_1,
    const Mat6&     tt6_2  ); 
    
    
double dot

  ( const Vector&   v6, 
    const Vec6&     t6 );    

double dot

  ( const Vec6&     t6, 
    const Vector&   v6 );

#endif





