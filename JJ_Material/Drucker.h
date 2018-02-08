/*
 * 
 *  Copyright (C) 2016 TU Delft. All rights reserved.
 *
 *  This model implements a Drucker-Prager material model. Regular return to the
 *  yield surface is done with an Euler Backward scheme as described in (1). The
 *  apex return requires special attention and uses an algorithm described in (2).
 *
 *    (1) "Non-Linear Finite Element Analysis of Solids and Structures",
 *        Volume 1, M.A. Crisfield 1997. 
 *    (2) "Computational Methods for Plasticity", E.A. de Souza Neto et.al. 2008
 *
 *
 *  Author  : E.C. Simons, e.c.simons@tudelft.nl
 *  Version : December 2016
 *
 *
 * Input parameters:
 *
 *  c   - initial cohesion / shear strength of the material
 *  h   - hardening parameter, will change the cohesion with plastic flow.
 *  phi - angle of internal friction (degrees)
 *  psi - angle of dilatancy (degrees)
 *
 */

#ifndef DRUCKER_H
#define DRUCKER_H

#include "jem/numeric/func/Function.h"
#include "jem/util/Flex.h"

#include "Plasticity.h"
#include "HookeMaterial.h"
#include "Invariants.h"

using jem::numeric::Function;
using jem::util::Flex;

// =======================================================
//  class Drucker
// =======================================================

class Drucker : public HookeMaterial,
                public Plasticity
{

  //********** GENERAL PROGRAM FUNCTIONS **********//
 public:

  typedef Drucker         Self;
  typedef HookeMaterial   Super;

  explicit                Drucker

    ( idx_t                 rank,
      const Properties&     globdat );

  virtual void            configure

    ( const Properties&     props );

  virtual void            getConfig

    ( const Properties&     conf )         const;

  virtual void            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      idx_t                 ipoint );

  void                    commit ();
  
  
  //********** HISTORY RELATED FUNCTIONS **********//

  void                    setHistory

    ( const Vec6&           epsp,
      const double          epspeq,
      const idx_t           ipoint );

  virtual void            getHistory

    ( const Vector&         hvals,
      const idx_t           mpoint ) const;
      
  virtual void            getStress

    ( const Vector&         stress,
      const idx_t           mpoint ) const;
      
  virtual void            getStrain

    ( const Vector&         strain,
      const idx_t           mpoint ) const;

  virtual void            getDissipationStress
    
    ( const Vector&         sstar,
      const Vector&         strain,
      const idx_t           ipoint ) const;

  virtual double          giveHistory     ( const idx_t ipoint ) const;

  virtual double          giveDissipation ( const idx_t point  ) const;

  virtual void            allocPoints

    ( const idx_t           count );

  inline virtual idx_t    pointCount () const;

  inline virtual idx_t    isLoading         ( idx_t point  ) const;

  inline virtual idx_t    wasLoading        ( idx_t point  ) const;
    
    
 //********** MATERIAL RELATED FUNCTIONS **********//
    
 void                    Yield
 
   (  const double&         epspeq,
            double&         sigY,
            double&         dsigY );
    
 void                    YieldSurf

    (       double&        sigeq,
            Vec6&          n,
            Vec6&          m,
      Tuple<double,6,6>&   dnormdsig,
      const Vec6&          sig );
      
 protected:

  virtual                ~Drucker   ();


 //********** MATERIAL RELATED VARIABLES **********//
 
 protected: 

  // history class, acts as a contained for all history variables
  class                  Hist_
  {
   public:
    Hist_();
    void                   toVector ( const Vector& vec ) const;
    void                   print () const;
    
    Mat6                   Dep;
    Vec6                   sig;
    Vec6                   eps;
    Vec6                   epsp;    // plastic strain
    double                 epspeq;  // equivalent plastic strain
    double                 damage;
    double                 dissipation;
    double                 sigeq;
    double                 p;
    int                    loading;

  };
  
  // preallocated history arrays
  Flex<Hist_>              preHist_;
  Flex<Hist_>              newHist_;
  Flex<Hist_>*             latestHist_;
  
  //********** MATERIAL RELATED VARIABLES **********//
  
  // need own rank, because base rank is always 3, for 3D stiffMat
  idx_t                    DruckerRank_;

  // input properties
  double                   c_;
  double                   phi_;
  double                   psi_;
  double                   h_;

  // dependent constant
  double                   G_;
  double                   K_;
  double                   alpha0_; // associated
  double                   alpha1_; // non-associated
  double                   tau_;
  

  

  // other parameters
  double                   rmTolerance_;
  idx_t                    rmMaxIter_;
  Properties               globdat_;
  double                   limsigfact_;
  double                   limsig_;

  // preallocated arrays
  Matrix                   P_;
};



//********** INLINE FUNCTIONS **********//

inline idx_t Drucker::pointCount () const
{
  return newHist_.size();
}

inline idx_t Drucker::isLoading ( idx_t ipoint ) const
{
  return newHist_[ipoint].loading;
}

inline idx_t Drucker::wasLoading ( idx_t ipoint ) const
{
  return preHist_[ipoint].loading;
}


#endif 
