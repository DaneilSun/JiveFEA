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
#include <jem/numeric/algebra/LUSolver.h>
#include <jem/numeric/utilities.h>
#include <jem/base/System.h>
#include <jive/util/FuncUtils.h>

#include "utilities.h"
#include "utilitiesTuple.h"
#include "LinHardPlast.h"

using namespace jem;
using jem::numeric::dotProduct;
using jem::numeric::inverse;
using jem::numeric::matmul;
using jem::numeric::UserFunc;
using jem::io::PrintWriter;
using jive::util::FuncUtils;

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


LinHardPlast::LinHardPlast 

  ( idx_t rank, const Properties& globdat )
    : Super ( 3, globdat ), LinHardPlastRank_ ( rank ), globdat_ ( globdat )

{
  rmTolerance_ = 1.e-10;
  rmMaxIter_   = 100;
  limsigfact_  = 1e-6;

  v61_.resize ( 6 );
  v62_.resize ( 6 );
  m6_.resize ( 6, 6 );

}


LinHardPlast::~LinHardPlast ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void LinHardPlast::configure ( const Properties& props )
{
  // configure LinHardPlast
  Super::configure ( props );

  // configure LinHardPlast
  props.find ( rmTolerance_, "rmTolerance" );
  props.find ( rmMaxIter_, "rmMaxIter" );
  props.find ( limsigfact_, "limsig" );
  
  sigmaC_ = makeFunc_ ( props, "sigmaC" );
  limsig_ = limsigfact_ * (sigmaC_->eval  ( 0.0 ));


  G_ = young_ / 2. / ( 1. + poisson_ );
  K_ = young_ / 3. / ( 1. - 2. * poisson_ );


  // need to read state, because Super is 3D
  if ( LinHardPlastRank_ == 2  )
  {
    props.get( stateString_, STATE_PROP );

    if      ( stateString_ == "PLANE_STRAIN" )
    {
      state_ = PlaneStrain;
    }
    else if ( stateString_ == "PLANE_STRESS" )
    {
      state_ = PlaneStress;
      
      throw Error( JEM_FUNC, "LinHardPlast does not have a consistent implementation of plane stress plasticity." );
    }
    else if ( stateString_ == "AXISYMMETRIC" )
    {
      state_ = AxiSymmetric;
    }
  }

  historyNames_.resize ( 9 );
  historyNames_[0] = "epsp_xx";
  historyNames_[1] = "epsp_yy";
  historyNames_[2] = "epsp_zz";
  historyNames_[3] = "epsp_xy";
  historyNames_[4] = "epsp_yz";
  historyNames_[5] = "epsp_zy";
  historyNames_[6] = "epspeq";
  historyNames_[7] = "diss";
  historyNames_[8] = "loading";

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

Ref<Function> LinHardPlast::makeFunc_ 

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


void LinHardPlast::getConfig ( const Properties& conf ) const
{
  Super::getConfig ( conf );

  conf.set ( "rmMaxIter"  , rmMaxIter_   );
  conf.set ( "rmTolerance", rmTolerance_ );
  conf.set ( "limsig"     , limsigfact_  );

  FuncUtils::getConfig ( conf, sigmaC_, "sigmaC" );

  if ( LinHardPlastRank_ == 2 )
  {
    conf.set ( STATE_PROP, stateString_ );
  }
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void LinHardPlast::update

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
  // CT5142 lecture notes, ``Computational Methods in Non-linear Solid 
  // Mechanics'', May 2010.
  // Book by M.A. Crisfield, "Non-Linear Finite Element Analysis of Solids and 
  // Structures", Volume 1, 1997. 
  
  // initialize doubles
  double epspeq, epspeq0, depspeq, f, G0, sigeq;
  epspeq = epspeq0 = depspeq = f = G0 = sigeq = 0.0;
  
  // initialize vectors
  Vec6 eps, eps0, deps, epsp0, depsp, epsp, sig, sig0, dsig, sig_b, sig_c, dsig_c, r_c, r_0, norm_b, norm_c;
  eps = eps0 = deps = epsp0 = depsp = epsp = sig = sig0 = dsig = sig_b = sig_c = dsig_c = r_c = norm_b = norm_c = 0.0;
  
  // initialize tuples
  Tuple<double,6,6>  C, dmat, dnormdsig, ident, Q1, H0;
  C = dmat = dnormdsig = ident = Q1 = H0 = 0.0;
  ident(0,0) = ident(1,1) = ident(2,2) = ident(3,3) = ident(4,4) = ident(5,5) = 1.0;
  Matrix  Q (6,6), H (6,6);
  Q = H = 0.0;
  m6_to_tt6(C,stiffMat_);
  
  
  // convert input strain vector of variable length to 
  // six component strain vector
  eps     = fill3DStrain ( strain );
  
  // read history values
  sig0    = preHist_[ipoint].sig;
  eps0    = preHist_[ipoint].eps;
  epsp0   = preHist_[ipoint].epsp;
  epspeq0 = preHist_[ipoint].epspeq;
  G0      = preHist_[ipoint].dissipation;
  
  // compute step in strain and stress
  deps = eps - eps0;
  dsig = matmul ( C, deps );
  
  // trial stress
  sig_b = sig0 + dsig;
  
  // find current yield stress and slope of the softening curve
  double sigC = Yield ( epspeq0 );
  double HC   = dYield( epspeq0 );
  YieldSurf(sigeq,norm_b,dnormdsig,sig_b);
  f = sigeq - sigC;
  
  if ( f >= -rmTolerance_ ) // Plastic step!
  {
    // simple radial return step to compute initial plastic multiplier
    double Delta_lambda = f / ( dot(norm_b,matmul(C,norm_b)) + HC );
    
    // compute new stress state, based on trial stress sig_b
    sig_c = sig_b - Delta_lambda*matmul(C,norm_b);
    
    // compute step in plastic strain and equivalent plastic strain
    depsp = Delta_lambda * norm_b;
    depspeq = MisesE(depsp);
    epspeq = epspeq0 + depspeq;
    
    // compute new yieldstress and slope of softening curve
    sigC = Yield ( epspeq );
    HC   = dYield( epspeq );
    
    // compute equivalent stress, yield function and its derivatives
    YieldSurf(sigeq,norm_c,dnormdsig,sig_c);
    f = sigeq - sigC;
    
    // initialize some more values
    int conv = 0;
    int iiter=0;

    // local NR scheme to find converged stress state on the yield function
    while ( conv == 0 && iiter < 25)
    {
      iiter++;

      // compute stress residual
      // Crisfield (6.79)
      r_c = sig_c-(sig_b-Delta_lambda*matmul(C,norm_c));
      
      // compute 'Q' matrix and its inverse
      // Crisfield (6.81)
      Q1 = ident + Delta_lambda * matmul(C,dnormdsig);
      tt6_to_m6(Q,Q1);
      Matrix Qinv (6,6);
      Qinv = inverse(Q);
      
      // compute change of plastic multiplier increment
      // Crisfield (6.83)
      double lambda_dot;
      lambda_dot = (f - dot(norm_c,tmatmul(Qinv,r_c))) / ( dot(norm_c,tmatmul(Qinv,tmatmul(stiffMat_,norm_c))) + HC );

      // update plastic multiplier increment
      Delta_lambda += lambda_dot;
      
      // compute change of stress increment
      // Crisfield (6.81)
      dsig_c = -tmatmul(Qinv,r_c) - lambda_dot*tmatmul(Qinv,tmatmul(stiffMat_,norm_c));
      
      // update stress increment
      sig_c += dsig_c;
      
      // compute step in plastic strain and equivalent plastic strain
      depsp   = Delta_lambda * norm_c;
      depspeq = MisesE(depsp);
      epspeq  = epspeq0 + depspeq;
      
      // compute new yieldstress and slope of softening curve
      sigC = Yield ( epspeq );
      HC   = dYield( epspeq );
      
      // compute equivalent stress, yield function and its derivatives
      YieldSurf(sigeq,norm_c,dnormdsig,sig_c);
      f = sigeq - sigC;
      
      // check error
      double err = std::abs(lambda_dot/Delta_lambda);

      // decide to terminal the local NR scheme
      if ( err < 1e-10 || std::abs(lambda_dot) < 1e-14 )
      {
        conv = 1;
      }
      else if ( iiter == 24 && std::abs(lambda_dot) > 1e-10 )
      {
        System::warn() << "Local NR failed to converge?!..... Point = " << ipoint << endl;
        System::out()  << "**************************************************" << endl;
        System::out()  << "Yield function = " << f << endl;
        System::out()  << "dLambda        = " << Delta_lambda << endl;
        System::out()  << "ddLambda       = " << lambda_dot << endl;
        System::out()  << "rel. error     = " << err << endl;
        System::out()  << "-------------------------------------------------" << endl;
        System::out()  << "total   strain      = " << eps     << endl;
        System::out()  << "plastic strain      = " << epsp0   << endl;
        System::out()  << "softening parameter = " << HC      << endl;
        System::out()  << "initial epspeq      = " << epspeq0 << endl;
        System::out()  << "**************************************************" << endl;
        throw Error( JEM_FUNC, "LinHardPlast could not find convergence in local NR scheme..." );
      }
      
    }


    // Build consistent tangent matrix
    // Lecture notes (eq 6.15)
    Q1 = ident + Delta_lambda * matmul(C,dnormdsig);
    tt6_to_m6(Q,Q1);   
    
    // Lecture notes (eq 6.20)
    H = matmul(inverse(Q),stiffMat_);
    m6_to_tt6(H0,H);
    
    // Lecture notes (eq 6.23)
    dmat = H0 - tmatmul(tmatmul(H,norm_c),tmatmul(norm_c,H)) / ( dot(norm_c,tmatmul(H,norm_c)) + HC );

    
    // Update history variables
    Vec6 depsp;
    depsp = Delta_lambda * norm_c;
    double dG = dot ( sig, depsp );
    sig = sig_c;
    newHist_[ipoint].sig = sig_c;
    newHist_[ipoint].eps = eps;
    newHist_[ipoint].epsp = epsp0 + depsp;
    newHist_[ipoint].epspeq = epspeq;
    newHist_[ipoint].dissipation = G0 + dG;
    newHist_[ipoint].loading = true;

  }
  else // Elastic step!
  {
    //System::out() << "Elastic step!\n";
    sig = sig_b;
    
    for ( idx_t i = 0; i < 6; ++i )
    {
      for ( idx_t j = 0; j < 6; ++j )
      {
        dmat(i,j) = stiffMat_(i,j);
      }
    }

    // Update history variables
    newHist_[ipoint].sig = sig;
    newHist_[ipoint].eps = eps;
    newHist_[ipoint].epsp = epsp0;
    newHist_[ipoint].epspeq = epspeq0;
    newHist_[ipoint].dissipation = G0;
    newHist_[ipoint].loading = false;
    
  }

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

void LinHardPlast::getDissipationStress

    ( const Vector&   sstar,
      const Vector&   strain,
      const idx_t     ipoint ) const

{
  System::out() << "No dissipation-stress algorithm for plasticity implemented!\n";
}

//-----------------------------------------------------------------------
//   getHistory
//-----------------------------------------------------------------------

void LinHardPlast::getHistory

  ( const Vector&  hvals,
    const idx_t    mpoint ) const

{
  (*latestHist_)[mpoint].toVector ( hvals );
}

//-----------------------------------------------------------------------
//   getStress
//-----------------------------------------------------------------------

void LinHardPlast::getStress

  ( const Vector&  stress,
    const idx_t    mpoint ) const

{
  reduce3DVector(stress, (*latestHist_)[mpoint].sig );
}

//-----------------------------------------------------------------------
//   getStrain
//-----------------------------------------------------------------------

void LinHardPlast::getStrain

  ( const Vector&  strain,
    const idx_t    mpoint ) const

{
  reduce3DVector(strain, (*latestHist_)[mpoint].eps );
}

//-----------------------------------------------------------------------
//   setHistory
//-----------------------------------------------------------------------

void LinHardPlast::setHistory

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

void LinHardPlast::commit () 

{
  //System::out() << "Erik: Commit\n";
  newHist_.swap ( preHist_ );

  latestHist_ = &preHist_;
}

//-----------------------------------------------------------------------
//   deviatoric
//-----------------------------------------------------------------------

Vec6 LinHardPlast::deviatoric

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

void LinHardPlast::isodev

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

void LinHardPlast::setPlaneStress

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

double LinHardPlast::Yield

    ( const double&           epspeq )
{
  double sigC  = sigmaC_->eval  ( epspeq );
  
  (sigC  < limsig_) ? sigC = limsig_ : 0.0;
  
  return sigC;
}

//-----------------------------------------------------------------------
//   dYield
//-----------------------------------------------------------------------

double LinHardPlast::dYield

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

void LinHardPlast::YieldSurf

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

double LinHardPlast::giveHistory ( const idx_t ip ) const
{
  return (*latestHist_)[ip].epspeq;
}

//-----------------------------------------------------------------------
//   allocPoints
//-----------------------------------------------------------------------

void LinHardPlast::allocPoints

  ( idx_t count )

{
  //System::out() << "Erik: allocPoints\n";

  for ( idx_t i = 0; i < count; ++i )
  {
    preHist_.pushBack ( Hist_() );
    newHist_.pushBack ( Hist_() );
  }
}

//-----------------------------------------------------------------------
//   giveDissipation
//-----------------------------------------------------------------------

double LinHardPlast::giveDissipation ( const idx_t ipoint ) const
{   
  return (*latestHist_)[ipoint].dissipation;
}

//-----------------------------------------------------------------------
//   Hist_ constructor
//-----------------------------------------------------------------------

LinHardPlast::Hist_::Hist_ () : epspeq ( 0. ), dissipation ( 0. ), sigeq (0.), p (0.)
{
  //System::out() << "Erik: contructor\n";
  epsp = 0.0;
  eps  = 0.0;
  sig  = 0.0;
  loading = false;
}

//-----------------------------------------------------------------------
//   Hist_ print function
//-----------------------------------------------------------------------

void LinHardPlast::Hist_::print () const

{
  System::out() << "epsp " << epsp << ", epspeq " << epspeq << endl;
}

// -------------------------------------------------------------------
//  Hist_::toVector
// -------------------------------------------------------------------

inline void LinHardPlast::Hist_::toVector

 ( const Vector&  vec ) const

{
  vec[0] = sig[0];
  vec[1] = sig[1];
  vec[2] = sig[2];
  vec[3] = sig[3];
  vec[4] = sig[4];
  vec[5] = sig[5];
  vec[6] = epspeq;
  vec[7] = dissipation;
  vec[8] = sigeq;
  vec[9] = p;
}

