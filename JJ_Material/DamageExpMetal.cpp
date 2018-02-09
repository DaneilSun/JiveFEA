/*
 * 
 *  Copyright (C) 2014 TU Delft. All rights reserved.
 *  
 *  The original structure by this code was by:
 *
 *  Author: F.P. van der Meer, f.p.vandermeer@tudelft.nl
 *  Date: October 2014
 *
 *  The original class implemented the material model for polymers
 *  from Melro et al. (2013)
 *  
 *  The Melro model was removed and replaced by simple J2 plasticity by:
 *
 *  Author: E.C. Simons, e.c.simons@tudelft.nl
 *  Date: March 2015
 *
 *  The J2 plasticity was removed and replaced by isotropic elasticity-based damage by:
 *
 *  Author: R. Bharali, ritukesh.bharali@gmail.com
 *  Date: February 2018
 *
 */


#include <jem/base/limits.h>
#include <jem/base/Float.h>
#include <jem/base/Error.h>
#include <jem/base/PrecheckException.h>
#include <jem/io/PrintWriter.h>
#include <jem/numeric/func/UserFunc.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/EigenUtils.h>
#include <jem/numeric/algebra/LUSolver.h>
#include <jem/numeric/utilities.h>
#include <jem/base/System.h>
#include <jive/util/FuncUtils.h>

#include "utilities.h"
#include "utilitiesTuple.h"
#include "DamageExpMetal.h"

using namespace jem;
using jem::numeric::dotProduct;
using jem::numeric::inverse;
using jem::numeric::matmul;
using jem::numeric::EigenUtils;
using jem::numeric::UserFunc;
using jem::io::PrintWriter;
using jive::util::FuncUtils;

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


DamageExpMetal::DamageExpMetal 

  ( idx_t rank, const Properties& globdat )
    : Super ( 3, globdat ), IsoDamageRank_ ( rank ), globdat_ ( globdat )

{
  rmTolerance_ = 1.e-10;
  //rmMaxIter_   = 100;
  //limsigfact_  = 1e-6;
  kappa1_      = 1.e-10;
  kappa2_      = 1.e-01; 

  v61_.resize ( 6 );
  v62_.resize ( 6 );
  m6_.resize ( 6, 6 );

}


DamageExpMetal::~DamageExpMetal ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void DamageExpMetal::configure ( const Properties& props )
{
  // configure DamageExpMetal
  Super::configure ( props );

  // configure DamageExpMetal
  props.find ( rmTolerance_, "rmTolerance" );
   
  // Added for Damage Model
  props.find ( kappa1_, "kappa1" );
  props.find ( kappa2_, "kappa2" );


  G_ = young_ / 2. / ( 1. + poisson_ );
  K_ = young_ / 3. / ( 1. - 2. * poisson_ );


  // need to read state, because Super is 3D
  if ( IsoDamageRank_ == 2  )
  {
    props.get( stateString_, STATE_PROP );

    if      ( stateString_ == "PLANE_STRAIN" )
    {
      state_ = PlaneStrain;
    }
    else if ( stateString_ == "PLANE_STRESS" )
    {
      state_ = PlaneStress;
      
      throw Error( JEM_FUNC, "DamageExpMetal does not have a consistent implementation of plane stress damage." );
    }
    else if ( stateString_ == "AXISYMMETRIC" )
    {
      state_ = AxiSymmetric;
    }
  }

  historyNames_.resize ( 12 );
  historyNames_[0] = "epsp_xx";
  historyNames_[1] = "epsp_yy";
  historyNames_[2] = "epsp_zz";
  historyNames_[3] = "epsp_xy";
  historyNames_[4] = "epsp_yz";
  historyNames_[5] = "epsp_zy";
  historyNames_[6] = "epspeq";
  historyNames_[7] = "diss";
  historyNames_[8] = "loading";
  historyNames_[9] = "p";
  historyNames_[10] = "history";
  historyNames_[11] = "damage";
  

  // Build P and Q matrix
  P_.resize(6,6);
  Q_.resize(6,6);
  P_ = Q_ = 0.0;
  for ( int i=0; i<3; i++ )
  {
    for ( int j=0; j<3; j++ )
    {
      i==j ? P_(i,j)=2./3.: P_(i,j)=-1./3.;
      i==j ? Q_(i,j)=2./3.: Q_(i,j)=-1./3.;
    }
  }
  P_(3,3) = P_(4,4) = P_(5,5) = 2.0;
  Q_(3,3) = Q_(4,4) = Q_(5,5) = 0.5;

}

//-----------------------------------------------------------------------
//   makeFunc_
//-----------------------------------------------------------------------

Ref<Function> DamageExpMetal::makeFunc_ 

  ( const Properties& props,
    const String&     name ) const

{
  // A function is created to define the hardening/softening behaviour. Input 
  // for this function is given as 'x' in the input file.

  String args = "x";
  props.find ( args, "args" );

  Ref<Function> func = FuncUtils::newFunc ( args, name, props, globdat_ );

  FuncUtils::resolve ( *func, globdat_ );

  return func;
}

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void DamageExpMetal::getConfig ( const Properties& conf ) const
{
  Super::getConfig ( conf );

  conf.set ( "rmTolerance", rmTolerance_ );
  conf.set ( "kappa1"     , kappa1_  );
  conf.set ( "kappa2"     , kappa2_  );

  if ( IsoDamageRank_ == 2 )
  {
    conf.set ( STATE_PROP, stateString_ );
  }
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void DamageExpMetal::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      idx_t                 ipoint )

{ 
  // The update routine is called when building the stiffness matrix and the
  // force vector in the SolidModel.
  // 
  // Input: 
  //   strain - vector containing strain values
  //   ipoint - integration point number, used to read and write history values
  //
  // Output:
  //   stress - vector contraining stress values, same size as strain
  //   stiff  - consistent tangent, square array
  //
  // The algorithm presented in this function is taken from:
  // Book by R. de Borst et. al., "Non-Linear Finite Element Analysis of Solids and 
  // Structures", Volume 1, 2012. 
  
  // initialize doubles
  double epspeq, sigeq, histvar , damage , hist0 , damage0 , epseq , damexp , deriv;;
  histvar = damage = hist0 = damage0 = epseq = damexp = deriv = epspeq = sigeq = 0.0;
  
  // initialize vectors
  Vec6 eps, epsp, sig;
  eps = epsp = sig = 0.0;
  
  // initialize tuples
  Tuple<double,6,6>  C, dmat;
  C = dmat = 0.0;
  m6_to_tt6(C,stiffMat_);
  
  
  // convert input strain vector of variable length to 
  // six component strain vector
  eps     = fill3DStrain ( strain );
  
  // read history values
  hist0   = preHist_[ipoint].histvar;

  // compute the equivalent strain 
  epseq  =  std::abs(1./2.*dot(eps,tmatmul(stiffMat_,eps))) ;

  // update the history parameter
  
  histvar = jem::max(hist0, epseq);

  if ( histvar <= kappa1_ )
      {
        // No damage
		    damage = 0.0;
	    }
  else
      {
        // Exponential damage formulation.

        damexp = std::exp ( -(histvar - kappa1_) / (kappa2_ - kappa1_) );

        damage = 1.0 - kappa1_ / histvar *  damexp;

        damage = jem::max ( damage, 0.0 );

        damage = jem::min ( damage, 1.0 );		
	    }

	// Update history variables
    
	newHist_[ipoint].histvar = histvar;
	newHist_[ipoint].damage = damage;
  newHist_[ipoint].loading = false;
    
	// Compute stress 
		
	dmat = (1 - damage) * C;

	sig = matmul(dmat,eps);
        
	if (histvar > hist0 && histvar > kappa1_)
      {
			
        // Compute derivate of the damage function
		    
			  damexp    = std::exp ( -(histvar - kappa1_) / (kappa2_ - kappa1_) );
			  deriv     = -(1.0 / histvar + 1.0 / (kappa2_ - kappa1_));
			  deriv    *= kappa1_ * damexp / histvar;
			
        // Build consistent tangent stiffness matrix 

        dmat     -= (1/histvar) * deriv * tmatmul(sig,sig);
			  newHist_[ipoint].loading = true;
		  }
		
  // Update history variables
    
  newHist_[ipoint].sig = sig;
  newHist_[ipoint].eps = eps;
  newHist_[ipoint].epsp = 0;
  newHist_[ipoint].epspeq = 0;
  newHist_[ipoint].dissipation = 0;
  newHist_[ipoint].sigeq = Mises(sig);
  newHist_[ipoint].p     = pressure(sig);
  
  // Put stress and stiffness is correct size (1d / 2d / 3d)
  reduce3DVector ( stress, sig );
  reduce3DMatrix ( stiff, dmat );

  // set latest history to current one
  latestHist_ = &newHist_;   
}


//-----------------------------------------------------------------------
//   getDissipationStress
//-----------------------------------------------------------------------

void DamageExpMetal::getDissipationStress

    ( const Vector&   sstar,
      const Vector&   strain,
      const idx_t     ipoint ) const

{
  System::out() << "No dissipation-stress algorithm implemented!\n";
}

//-----------------------------------------------------------------------
//   getHistory
//-----------------------------------------------------------------------

void DamageExpMetal::getHistory

  ( const Vector&  hvals,
    const idx_t    mpoint ) const

{
  (*latestHist_)[mpoint].toVector ( hvals );
}

//-----------------------------------------------------------------------
//   getStress
//-----------------------------------------------------------------------

void DamageExpMetal::getStress

  ( const Vector&  stress,
    const idx_t    mpoint ) const

{
  reduce3DVector(stress, (*latestHist_)[mpoint].sig );
}

//-----------------------------------------------------------------------
//   getStrain
//-----------------------------------------------------------------------

void DamageExpMetal::getStrain

  ( const Vector&  strain,
    const idx_t    mpoint ) const

{
  reduce3DVector(strain, (*latestHist_)[mpoint].eps );
}

//-----------------------------------------------------------------------
//   setHistory
//-----------------------------------------------------------------------

void DamageExpMetal::setHistory

  ( const Vec6&    epsp,
    const double   epspeq,
    const idx_t    ipoint )

{
  preHist_[ipoint].epspeq = epspeq;
  preHist_[ipoint].epsp   = epsp;
}



//-----------------------------------------------------------------------
//   commit
//-----------------------------------------------------------------------

void DamageExpMetal::commit () 

{
  newHist_.swap ( preHist_ );

  latestHist_ = &preHist_;
}

//-----------------------------------------------------------------------
//   deviatoric
//-----------------------------------------------------------------------

Vec6 DamageExpMetal::deviatoric

  ( const Vec6&   full,
    const double  p ) const

{
  Vec6 ret = full;

  ret[0] -= p;
  ret[1] -= p;
  ret[2] -= p;

  return ret;
}

//-----------------------------------------------------------------------
//   split into pressure and deviatoric
//-----------------------------------------------------------------------

void DamageExpMetal::isodev

  ( const Vec6&   all,
    Vec6&   iso,
    Vec6&   dev  ) const

{
  iso = dev = 0.0;
  
  double p = (all[0] + all[1] + all[2])*1./3.;

  iso[0] = iso[1] = iso[2] = -p;
  
  dev = all + iso;
}

//-----------------------------------------------------------------------
//   set stress vector to plane stress vector
//-----------------------------------------------------------------------

void DamageExpMetal::setPlaneStress

    (      Vec6&           sig )
{
  if ( state_ == 1 ) // Plane stress state
  {
    sig[2] = 0.0;
    sig[4] = 0.0;
    sig[5] = 0.0;
  }
  else //if ( state_ == 2 )
  {
    sig[4] = 0.0;
    sig[5] = 0.0;
  }
}

//-----------------------------------------------------------------------
//   Yield
//-----------------------------------------------------------------------

double DamageExpMetal::Yield

    ( const double&           epspeq )
{
  double sigC  = sigmaC_->eval  ( epspeq );
  
  (sigC  < limsig_) ? sigC = limsig_ : 0.0;
  
  return sigC;
}

//-----------------------------------------------------------------------
//   dYield
//-----------------------------------------------------------------------

double DamageExpMetal::dYield

    ( const double&           epspeq )
{
  double sigC  = sigmaC_->eval   ( epspeq  );
  double dsigC = sigmaC_->deriv  ( epspeq  );
  
  (sigC  < limsig_) ? dsigC = 0.0 : 0.0;
  
  return dsigC;
}

//-----------------------------------------------------------------------
//   YieldSurf
//-----------------------------------------------------------------------

void DamageExpMetal::YieldSurf

    ( double&              sigeq,
      Vec6&                norm,
      Tuple<double,6,6>&   dnormdsig,
      const Vec6&          sig )
      
{
  sigeq = sqrt( std::abs(3./2.*dot(sig,tmatmul(P_,sig))) );
  
  norm = 3./2.*1./sigeq * tmatmul(P_,sig);

  for ( int i=0; i<6; i++ )
  {
    for ( int j=0; j<6; j++ )
    {
      dnormdsig(i,j) = (3./(2.*sigeq)) * P_(i,j) - 1./sigeq * norm[i] * norm[j];
    }
  }
}

//-----------------------------------------------------------------------
//   giveHistory
//-----------------------------------------------------------------------

double DamageExpMetal::giveHistory ( const idx_t ip ) const
{
  return (*latestHist_)[ip].epspeq;
}

//-----------------------------------------------------------------------
//   allocPoints
//-----------------------------------------------------------------------

void DamageExpMetal::allocPoints

  ( idx_t count )

{
  for ( idx_t i = 0; i < count; ++i )
  {
    preHist_.pushBack ( Hist_() );
    newHist_.pushBack ( Hist_() );
  }
}

//-----------------------------------------------------------------------
//   giveDissipation
//-----------------------------------------------------------------------

double DamageExpMetal::giveDissipation ( const idx_t ipoint ) const
{   
  return (*latestHist_)[ipoint].dissipation;
}

//-----------------------------------------------------------------------
//   Hist_ constructor
//-----------------------------------------------------------------------

DamageExpMetal::Hist_::Hist_ () : epspeq ( 0. ), dissipation ( 0. ), sigeq (0.), p (0.), histvar(0.), damage(0.)
{
  epsp = 0.0;
  eps  = 0.0;
  sig  = 0.0;
  loading = false;
  histvar = 0.0;
  damage  = 0.0;
}

//-----------------------------------------------------------------------
//   Hist_ print function
//-----------------------------------------------------------------------

void DamageExpMetal::Hist_::print () const

{
  System::out() << "epsp " << epsp << ", epspeq " << epspeq << endl;
}

// -------------------------------------------------------------------
//  Hist_::toVector
// -------------------------------------------------------------------

inline void DamageExpMetal::Hist_::toVector

 ( const Vector&  vec ) const

{
  vec[0] = sig[0];
  vec[1] = sig[1];
  vec[2] = sig[2];
  vec[3] = sig[3];
  vec[4] = sig[4];
  vec[5] = sig[5];
  vec[6] = 0;  // epspeq
  vec[7] = 0;  // dissipation
  vec[8] = sigeq;
  vec[9] = p;
  vec[10] = histvar;
  vec[11] = damage;
}

