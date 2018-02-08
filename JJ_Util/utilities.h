#ifndef UTILITIES_H
#define UTILITIES_H

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

// A pointer to a function that computes the spatial derivatives of
// the interpolation matrix. This is the so-called B-matrix.

typedef void        (*ShapeGradsFunc)

  ( const Matrix&       b,
    const Matrix&       g );

// Similar stuff, from the shape functions N_i, compute the shape matrix N

typedef void        (*ShapeFuncsFunc)

  ( const Matrix&       sfuncs,
    const Vector&       n );

// Similar stuff, from computing grad(U)

typedef void        (*GradUFunc)

  ( const Matrix&       gradU,
    const Matrix&       gradN,
    const Vector&       elemDisp );

//-----------------------------------------------------------------------
//   constants
//-----------------------------------------------------------------------

// An integer array that maps the number of spatial dimensions (1, 2,
// or 3) to the number of strain/stress components.

extern const idx_t    STRAIN_COUNTS[4];
extern const double   PI;

//-----------------------------------------------------------------------
//   public functions
//-----------------------------------------------------------------------

// -----------------------------------------------------------------------
//   matrix of shape functions N
// -----------------------------------------------------------------------

// A function that returns a pointer to a function that computes the
// N-matrix given the number of spatial dimensions.

ShapeFuncsFunc        getShapeFuncsFunc

  ( idx_t               rank );

void                  get1DShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n );

void                  get2DShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n );

void                  get3DShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n );
    


// -----------------------------------------------------------------------
//   matrix of shape function derivatives B
// -----------------------------------------------------------------------

// A function that returns a pointer to a function that computes the
// B-matrix given the number of spatial dimensions.

ShapeGradsFunc        getShapeGradsFunc

  ( idx_t               rank );

void                  get1DShapeGrads

  ( const Matrix&       b,
    const Matrix&       g );

void                  get2DShapeGrads

  ( const Matrix&       b,
    const Matrix&       g );

void                  get3DShapeGrads

  ( const Matrix&       b,
    const Matrix&       g );


// -----------------------------------------------------------------------
//   other utility functions
// -----------------------------------------------------------------------

// 

void                  addAxiShapeGrads

  ( const Matrix&       b,
    const Vector&       n,
    const double&       r  );

// 2d Cauchy stress vector to 2d Cauchy stress matrix.

void                  Cauchy2DVec2Mat

  ( const Matrix&       CMat,
    const Vector&       CVec );

// Solving linear equation of one variable
// a*x + b = 0

void                  solveLinearEqua

  ( double&      ans,
    const double a,
    const double b );

// Solving quadratic equation
// a*x^2 + b*x + c = 0

void                  solveQuadEqua

  ( Vector&      ans,
    const double a,
    const double b,
    const double c );

// Solving cubic equation
// a*x^3 + b*x^2 + c*x + d = 0
// only for functions with real roots
// returns the number of roots

idx_t                 solveCubicEqua

  ( Tuple<double,3>& ans,
    const double     a,
    const double     b,
    const double     c,
    const double     d);

// Invert 2x2 matrix (no size check performed!)

void                  invert2x2

  ( const Matrix& mat );

// Compute the ramp function <x> = 1/2(x+x)

double                  evalMcAuley

  ( const double );

double                  evalHeaviside

  ( const double );

// Find unique values 'U' from a Vector 'V'
void                    getUnique

  (       Vector& U,
    const Vector& V );

void                    getUnique

  (       IdxVector& U,
    const IdxVector& V );
  
// Find positions 'pos' of value 'd' in vector 'V'.
void                    getThisValue

  (       IdxVector& pos,
    const Vector&    V,
    const double     d );

void                    getThisValue

  (       IdxVector& pos,
    const IdxVector& V,
    const double     d );
    
// Find position 'pos' of value 'd' in vector 'V', returns one value only! For 
// vectors with multiple occurences of 'd' the last occurence is stored in pos.
void                    getThisValue

  (       int&       pos,
    const Vector&    V,
    const double     d );
    
void                    getThisValue

  (       int&       pos,
    const IdxVector& V,
    const double     d );
    
// Normalize a vector with repect to its length (L_2 norm)
void                    makeNorm

  (       Vector&    norm );
  
// Returns the length (norm) of a vector
double                  getLength

  ( const Vector&    norm );
  
// Returns the sign of a double, may be obsolete.
double                  signERIK

  ( const double     d );
    
#endif

