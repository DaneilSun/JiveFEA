#ifndef UTILITIESL_H
#define UTILITIESL_H

#include <jive/Array.h>
#include <jem/base/Tuple.h>

#include "Array.h"

using jive::Vector;
using jive::Matrix;
using jive::IntVector;
using jem::Tuple;
using jem::idx_t;

//-----------------------------------------------------------------------
//   typedefs
//-----------------------------------------------------------------------

// Similar stuff, but now for Updated Lagrangian large strain Formulation
typedef void        (*NLShapeGradsFunc)

  ( const Matrix&       bNL,
    const Matrix&       g );

//-----------------------------------------------------------------------
//   constants
//-----------------------------------------------------------------------

// An integer array that maps the number of spatial dimensions (1, 2,
// or 3) to the number of strain/stress components.

extern const idx_t    STRAIN_COUNTS[4];


//-----------------------------------------------------------------------
//   public functions
//-----------------------------------------------------------------------
    
void                  getDeformationGradient
  
  ( const Cubix&       defgrad,
    const Matrix&      coords,
    const Vector&      disp,
    const Cubix&       grads );
    
void                 getCurrentgrads

  ( const Cubix&       deformedgrads,
    const Cubix&       grads,
    const Cubix&       defgrad );


// These functions compute the B-matrix given an interpolation matrix.
    
void                  get2DNLShapeGrads

  ( const Matrix&       bNL,
    const Matrix&       g );    

// A function that returns a pointer to a function that computes the
// B-matrix given the number of spatial dimensions.
  
NLShapeGradsFunc        getNLShapeGradsFunc

  ( idx_t               rank );

// -----------------------------------------------------------------------
//   other utility functions
// -----------------------------------------------------------------------

void                getStressTensor

  ( const Matrix&       S,
    const Vector&  stress );
    
void                getStressMatrix

  ( const Matrix&       S,
    const Vector&  stress );

void                getStressVector

  ( const Vector&  stress,
    const Matrix&       S );
    
#endif

