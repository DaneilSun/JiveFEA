/*
 * 
 *  Copyright (C) 2014 TU Delft. All rights reserved.
 *  
 *  This class implements the material model for polymers
 *  from Melro et al. (2013)
 *  
 *  Author: F.P. van der Meer, f.p.vandermeer@tudelft.nl
 *  Date: October 2014
 *
 */

#ifndef LINHARDPLAST_H
#define LINHARDPLAST_H

#include "jem/numeric/func/Function.h"
#include "jem/util/Flex.h"

#include "Plasticity.h"
#include "HookeMaterial.h"
#include "Invariants.h"

using jem::numeric::Function;
using jem::util::Flex;

// =======================================================
//  class LinHardPlast
// =======================================================

class LinHardPlast : public HookeMaterial,
                      public Plasticity
{
 public:

  typedef LinHardPlast   Self;
  typedef HookeMaterial   Super;

  explicit                LinHardPlast

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
  
  // History related functions

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
    
  // Simple output requests

  inline virtual idx_t    pointCount () const;

  inline virtual idx_t    isLoading         ( idx_t point  ) const;

  inline virtual idx_t    wasLoading        ( idx_t point  ) const;
  
  // Stress related operations

  Vec6                    deviatoric

    ( const Vec6&           full,
      const double          p ) const;
      
  void                    isodev

    ( const Vec6&           all,
      Vec6&                 iso,
      Vec6&                 dev  ) const;

  void                    setPlaneStress

    (       Vec6&           sig );
    
 // Material model related functions
    
 double                  Yield

    ( const double&         epspeq ); 
    
 double                  dYield

    ( const double&         epspeq );
    
 void                    YieldSurf

    (       double&        sigeq,
            Vec6&          norm,
      Tuple<double,6,6>&   dnormdsig,
      const Vec6&          sig );
      
 protected:

  virtual                ~LinHardPlast   ();

  Ref<Function>           makeFunc_

    ( const Properties&     props,
      const String&         name )         const;

 // Variables ///////////////////////////////////////////////////
 
 protected: 

  class                  Hist_
  {
   public:
    Hist_();
    void                   toVector ( const Vector& vec ) const;
    void                   print () const;
    
    Vec6                   sig;
    Vec6                   eps;
    Vec6                   epsp;    // plastic strain
    double                 epspeq;  // equivalent plastic strain
    double                 dissipation;
    double                 sigeq;
    double                 p;
    bool                   loading;
  };

  // need own rank, because base rank is always 3, for 3D stiffMat

  idx_t                    LinHardPlastRank_;

  // input properties

  double                   rmTolerance_;
  idx_t                    rmMaxIter_;

  // dependent constants

  double                   G_;
  double                   K_;
  double                   alpha_;
  double                   alpha2_;
  double                   nuPFac_;
  double                   Ka_;

  // history 

  Flex<Hist_>              preHist_;
  Flex<Hist_>              newHist_;
  Flex<Hist_>*             latestHist_;

  // hardening functions

  Properties               globdat_;
  Ref<Function>            sigmaC_;
  double                   limsigfact_;
  double                   limsig_;

  // preallocated arrays

  Vector                   v61_;
  Vector                   v62_;
  Matrix                   m6_;
  Vector                   tmpm6_;
  Matrix                   P_;
  Matrix                   Q_;
};

inline idx_t LinHardPlast::pointCount () const
{
  return newHist_.size();
}

inline idx_t LinHardPlast::isLoading ( idx_t ipoint ) const
{
  return newHist_[ipoint].loading;
}

inline idx_t LinHardPlast::wasLoading ( idx_t ipoint ) const
{
  return preHist_[ipoint].loading;
}

#endif 
