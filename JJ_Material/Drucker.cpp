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
#include "Drucker.h"

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


Drucker::Drucker 

  ( idx_t rank, const Properties& globdat )
    : Super ( 3, globdat ), DruckerRank_ ( rank ), globdat_ ( globdat )

{
  rmTolerance_ = 1.e-10;
  rmMaxIter_   = 25;
  limsigfact_  = 1e-6;
}


Drucker::~Drucker ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void Drucker::configure ( const Properties& props )
{
  // configure Drucker-Prager
  Super::configure ( props );

  // configure Drucker-Prager
  props.find ( rmTolerance_, "rmTolerance" );
  props.find ( rmMaxIter_, "rmMaxIter" );
  props.find ( limsigfact_, "limsig" );
  props.find ( c_, "c" );
  props.find ( phi_, "phi" );
  props.find ( psi_, "psi" );
  props.find ( h_, "h" );

  G_ = young_ / 2. / ( 1. + poisson_ );
  K_ = young_ / 3. / ( 1. - 2. * poisson_ );
  
  alpha0_ = (6 * sin( 2.0*PI * phi_/360 ) ) / ( 3 - sin( 2.0*PI * phi_/360 ) );
  alpha1_ = (6 * sin( 2.0*PI * psi_/360 ) ) / ( 3 - sin( 2.0*PI * psi_/360 ) );
  tau_ = (6 * c_ * cos( 2.0*PI * phi_/360 ) ) / ( 3 - sin( 2.0*PI * phi_/360 ) );

  //System::out() << alpha0_ << " " << tau_ << endl;

  limsig_ = limsigfact_ * tau_;
  

  // need to read state, because Super is 3D
  if ( DruckerRank_ == 2  )
  {
    props.get( stateString_, STATE_PROP );

    if      ( stateString_ == "PLANE_STRAIN" )
    {
      state_ = PlaneStrain;
    }
    else if ( stateString_ == "PLANE_STRESS" )
    {
      state_ = PlaneStress;
      
      throw Error( JEM_FUNC, "Drucker does not have a consistent implementation of plane stress plasticity." );
    }
    else if ( stateString_ == "AXISYMMETRIC" )
    {
      state_ = AxiSymmetric;
    }
  }


  historyNames_.resize ( 10 );
  
  historyNames_[0] = "sigma_xx";
  historyNames_[1] = "sigma_yy";
  historyNames_[2] = "sigma_zz";
  historyNames_[3] = "sigma_xy";
  historyNames_[4] = "sigma_yz";
  historyNames_[5] = "sigma_zy";
  
  historyNames_[6] = "epspeq";
  historyNames_[7] = "loading";
  historyNames_[8] = "sigeq";
  historyNames_[9] = "p";

  
  // Build P matrix
  P_.resize(6,6);
  P_ = 0.0;
  for ( int i=0; i<3; i++ )
  {
    for ( int j=0; j<3; j++ )
    {
      i==j ? P_(i,j)=2./3.: P_(i,j)=-1./3.;
    }
  }
  P_(3,3) = P_(4,4) = P_(5,5) = 2.0;

}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void Drucker::getConfig ( const Properties& conf ) const
{
  Super::getConfig ( conf );

  conf.set ( "rmMaxIter"  , rmMaxIter_   );
  conf.set ( "rmTolerance", rmTolerance_ );
  conf.set ( "limsig"     , limsigfact_  );
  
  conf.set ( "G"      , G_  );
  conf.set ( "K"      , K_  );
  
  conf.set ( "c"      , c_  );
  conf.set ( "phi"    , phi_  );
  conf.set ( "psi"    , psi_  );
  conf.set ( "h"      , h_  );

  conf.set ( "alpha0" , alpha0_  );
  conf.set ( "alpha1" , alpha1_  );
  conf.set ( "tau"    , tau_  );

  if ( DruckerRank_ == 2 )
  {
    conf.set ( STATE_PROP, stateString_ );
  }
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void Drucker::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      idx_t                 ipoint )

{ //System::out() << "update for point " << ipoint << endl;
  // The current function implements a return mapping algorithm for Drucker 
  // Prager hardening plasticity. Softening is accepted, the hardening 
  // relation is currently only linear. Other types may be added if Yield is 
  // properly modified. 
  // The algorithm is an Euler backward scheme taken from:
  // "Non-Linear Finite Element Analysis of Solids and Structures",
  // Volume 1, M.A. Crisfield 1997. 
  // 
  // Input: 
  //  - strain : Integration point strains, can be length 3, 4 or 6
  //  - ipoint : Integration point number, used to find and store history
  //
  // Output:
  //  - stress : Integration point stress, same size as strain
  //  - stiff  : Consistent tangent
  // 
  
  
  // Initialize doubles
  double epspeq, epspeq0, depspeq, f, G0, sigeq, sigC, HC;
  epspeq = epspeq0 = depspeq = f = G0 = sigeq = sigC = HC = 0.0;
  
  // Initialize vectors
  Vec6 eps, eps0, deps, epsp0, depsp, epsp, sig, sig0, dsig, sig_b, sig_c, dsig_c, r_c, r_0, n_b, n_c, m_b, m_c, dfds;
  eps = eps0 = deps = epsp0 = depsp = epsp = sig = sig0 = dsig = sig_b = sig_c = dsig_c = r_c = n_b = n_c = m_b = m_c = dfds = 0.0;
  
  // Initialize matrices
  Tuple<double,6,6>  C, dmat, dnormdsig, ident, Q1;
  C = dmat = dnormdsig = ident = Q1 = 0.0;
  ident(0,0) = ident(1,1) = ident(2,2) = ident(3,3) = ident(4,4) = ident(5,5) = 1.0;
  Matrix  Q (6,6);
  Q = 0.0;
  m6_to_tt6(C,stiffMat_);
  
  
  // Construct full strain vector and obtain history data
  eps     = fill3DStrain ( strain );
  sig0    = preHist_[ipoint].sig;
  eps0    = preHist_[ipoint].eps;
  epsp0   = preHist_[ipoint].epsp;
  epspeq0 = preHist_[ipoint].epspeq;
  G0      = preHist_[ipoint].dissipation;
  
  // compute step in strain and stress
  deps = eps - eps0;
  dsig = matmul ( C, deps );
  
  // trial stress, based on elastic assumption
  sig_b = sig0 + dsig;
  
  // find current yield stress and slope of the softening curve  
  Yield(epspeq0,sigC,HC);
  YieldSurf(sigeq,n_b,m_b,dnormdsig,sig_b);
  f = sigeq - sigC;
        
  if ( f >= -rmTolerance_ ) // Plastic step!
  {   
    //System::out() << "Plastic step in Drucker!\n";
    newHist_[ipoint].loading = 1;
    
    // Compute initial return map. This step assumes a return to a smooth part
    // of the yield surface, with non-associated flow direction m_b
    double Delta_lambda = f / ( dot(n_b,matmul(C,m_b)) + HC );
    sig_c = sig_b - Delta_lambda*matmul(C,m_b);    
    
    // increment in plastic strain
    depsp   = Delta_lambda * m_b;
    depspeq = MisesE(depsp);
    epspeq  = epspeq0 + depspeq;
    
    // evaluate yield function and the derivatives
    Yield(epspeq,sigC,HC);
    YieldSurf(sigeq,n_c,m_c,dnormdsig,sig_c);
    f = sigeq - sigC;
    
    int conv  = 0;
    int iiter = 0;

    // Check if the initial return point still lies on the smooth part of the
    // yield surface. If not an apex return scheme is applied, otherwise regular
    // return can be used.
    bool chk0 = ( pressure(sig_c) > -sigC/alpha0_ );

    // Check for regular return
    if ( chk0 )
    {
      while ( conv == 0 && iiter < rmMaxIter_ )
      {
        iiter++;
      
        // stress residual
        r_c = sig_c-(sig_b-Delta_lambda*matmul(C,m_c));
        
        Q1 = ident + Delta_lambda * matmul(C,dnormdsig);
        tt6_to_m6(Q,Q1);
        Matrix Qinv (6,6);
        Qinv = inverse(Q);
      
        // Compute increment in plastic multiplier
        double lambda_dot;
        lambda_dot = (f - dot(n_c,tmatmul(Qinv,r_c))) / ( dot(n_c,tmatmul(Qinv,tmatmul(stiffMat_,m_c))) + HC );
        Delta_lambda += lambda_dot;
      
        // Compute increment in stress
        dsig_c = -tmatmul(Qinv,r_c) - lambda_dot*tmatmul(Qinv,tmatmul(stiffMat_,m_c));
        sig_c += dsig_c;
        
        // Update equivalent strain
        depsp   = Delta_lambda * m_c;
        depspeq = MisesE(depsp);
        epspeq  = epspeq0 + depspeq;
        
        // Update yield stress and yieldsurface related properties
        Yield(epspeq,sigC,HC);
        YieldSurf(sigeq,n_c,m_c,dnormdsig,sig_c);
        f = sigeq - sigC;
      
        // check for convergence
        double err = std::abs(lambda_dot/Delta_lambda);
      
        if ( err < 1e-10 || std::abs(lambda_dot) < 1e-14 || ( err < 1e-3 && HC == 0.0 ) ) 
        {
          conv = 1;
          if ( iiter > 5 )
          {
            System::out() << "Converged in " << iiter << " iterations.\n";
          }
        }
        else if ( iiter == rmMaxIter_-1 && std::abs(lambda_dot) > 1e-10 )
        {
          System::warn() << "Local NR failed to converge?!..... Point = " << ipoint << endl;
          System::out()  << "**************************************************" << endl;
          System::out()  << "Yield function = " << f << endl;
          System::out()  << "dLambda        = " << Delta_lambda << endl;
          System::out()  << "ddLambda       = " << lambda_dot << endl;
          System::out()  << "rel. error     = " << err << endl;
          System::out()  << "-------------------------------------------------" << endl;
          System::out()  << "total   strain = " << eps     << endl;
          System::out()  << "plastic strain = " << epsp0   << endl;
          System::out()  << "initial epspeq = " << epspeq0 << endl;
          System::out()  << "**************************************************" << endl;
          throw Error( JEM_FUNC, "Drucker could not find convergence in local NR scheme..." );
        }
        
      }
      
      // Build consistent tangent matrix
      Q1 = ident + Delta_lambda * matmul(C,dnormdsig);
      tt6_to_m6(Q,Q1);   
      
      Matrix H (6,6);
      Tuple<double,6,6> H0;
      H = 0.0; H0 = 0.0;
      H = matmul(inverse(Q),stiffMat_);
      m6_to_tt6(H0,H);
      
      Tuple<double,6,6> Const1;
      Const1 = H0 - tmatmul(tmatmul(H,m_c),tmatmul(n_c,H)) / ( dot(n_c,tmatmul(H,m_c)) + HC );
      
      for ( idx_t i = 0; i < 6; ++i )
      {
        for ( idx_t j = 0; j < 6; ++j )
        {
          dmat(i,j) = Const1(i,j);
        }
      }
    }
    else
    {
      System::debug() << "==============================\n";
      System::debug() << "Apex return for Drucker-Prager! " << ipoint << "\n";
      newHist_[ipoint].loading = 2;
      
      // ECS 2016-11-10
      //depsp = 0.0;
      //
      
      double A, B, D, depv, res, res0, ptr, p;
      A = 1.0 / alpha0_;
      B = 1.0 / alpha1_;
      D = depv = ptr = 0.0;
      
      // trial pressure and residual
      Yield(epspeq0,sigC,HC);
      ptr = pressure(sig_b);
      res = sigC*A + ptr;
      res0 = res;

      conv = 0;
      
      while ( conv == 0  && alpha1_ > 0.0 )
      {
        // alpha1_ should be larger than zero. For alpha1_ = 0.0 the plastic 
        // flow is pure deviatoric. Once the stress reaches an isobaric state
        // no additional plastic strain can be generated. This means no more 
        // softening can occur and the material can never reach the Apex point.
        iiter++;
      
        D      = A*B*HC + K_;
        depv  -= res / D;
        epspeq = epspeq0 + B*depv;
        p      = ptr + K_*depv;
        
        Yield(epspeq,sigC,HC);
        res = sigC*A + p;
        
        sig_c = 0.0;
        sig_c[0] = sig_c[1] = sig_c[2] = -p;
        
        if ( std::abs(res) < 1e-10 )
        {
          conv = 1;
        }
        else if ( iiter > 100 )
        {
          System::debug() << "Apex return does not seem to converge?" << endl;
          
          epspeq = 100.0;
          Yield(epspeq,sigC,HC);
          sig_c[0] = sig_c[1] = sig_c[2] = sigC*A;
          
          conv = 1;
          //throw Error( JEM_FUNC, "Apex return does not seem to converge?" );
        }
        
      }

      dmat = 0.0;
      // build consistent tangent matrix
      for ( idx_t i = 0; i < 3; ++i )
      {
        for ( idx_t j = 0; j < 3; ++j )
        {
          dmat(i,j) = K_ * ( 1 - K_ / ( K_ + A*B*HC ) );
        }
      }
      
      if ( alpha1_ == 0.0 )
      {
        // When apex return is invoked it means the initial returned stress 
        // has reached an isobaric state, whilst still being outside the yield 
        // surface. Plastic deformation should be volumetric to return to the 
        // yield surface. For alpha1_ == 0.0 this is however impossible. 
        // Three cases are considered:
        //
        //  h < 0.0    - Full softening is applied
        //  h = 0.0    - Return to (fixed) apex point
        //  h > 0.0    - Apex expanded to trial pressure state
        //
        
        System::debug() << "Alpha1_ == 0.0! Return theoretically impossible..." << endl;
        dmat   = 0.0;
        epspeq = 100.0;
        sig_c  = 0.0;
        
        if ( h_ == 0.0 )
        {
          System::debug() << "h = 0.0  -  returned to fixed apex" << endl;
          Yield(epspeq,sigC,HC);
          sig_c[0] = sig_c[1] = sig_c[2] = sigC*A;
        }
        else if ( h_ > 0.0 )
        {
          System::debug() << "h > 0.0  -  expanded surface to match trial pressure" << endl;
          double sigC0, sigC1;
          Yield(epspeq0,sigC0,HC);
          sigC1 = -ptr / A;
          depspeq = (sigC1-sigC0)/HC;
          epspeq = epspeq0 + depspeq;
          Yield(epspeq,sigC,HC);
          sig_c[0] = sig_c[1] = sig_c[2] = sigC*A;
        }
        else
        {
          System::debug() << "h < 0.0  -  full softening is applied" << endl;
        }
      }
      else if ( HC <= -(K_/(A*B)) )
      {
        System::debug() << "Critical softening detected, full failure applied" << endl;
        dmat   = 0.0;
        epspeq = std::abs(tau_/h_);//100.0;
        sig_c  = 0.0;
      }

    }
    
    
    // Update history variables
    depsp = Delta_lambda * m_c;
    double dG = dot ( sig, depsp );
    sig = sig_c;
    newHist_[ipoint].sig         = sig_c;
    newHist_[ipoint].eps         = eps;
    newHist_[ipoint].epsp        = epsp0 + depsp;
    newHist_[ipoint].epspeq      = epspeq;
    newHist_[ipoint].dissipation = G0 + dG;
    newHist_[ipoint].Dep         = dmat;

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
    newHist_[ipoint].sig         = sig;
    newHist_[ipoint].eps         = eps;
    newHist_[ipoint].epsp        = epsp0;
    newHist_[ipoint].epspeq      = epspeq0;
    newHist_[ipoint].dissipation = G0;
    newHist_[ipoint].Dep         = dmat;
    newHist_[ipoint].loading     = 0;
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

//-----------------------------------------------------------------------
//   getDissipationStress
//-----------------------------------------------------------------------

void Drucker::getDissipationStress

    ( const Vector&   sstar,
      const Vector&   strain,
      const idx_t     ipoint ) const

{
  // Used for the dissipation based arclength scheme by Verhoosel
  // this function implements the integrand of eq(30) and is called by 
  // the SolidModel to build the full integral.
  // ECS: Please note this function is called after COMMIT. The values stored
  // in preHist_ are therefore the latest values
  Mat6 Dep, DepT, De, DeInv;
  Vec6 sig, sstarT;
  
  Dep  = preHist_[ipoint].Dep;//(*latestHist_)[ipoint].Dep;
  DepT = Dep.transpose();
  
  m6_to_tt6(De,stiffMat_);
  DeInv = inverse(De);  
  
  sig = preHist_[ipoint].sig;//(*latestHist_)[ipoint].sig;
  
  sstarT = tmatmul(DepT,tmatmul(DeInv,sig));
  
  reduce3DVector ( sstar, sstarT );
}

//-----------------------------------------------------------------------
//   getHistory
//-----------------------------------------------------------------------

void Drucker::getHistory

  ( const Vector&  hvals,
    const idx_t    mpoint ) const

{
  (*latestHist_)[mpoint].toVector ( hvals );
}

//-----------------------------------------------------------------------
//   getStress
//-----------------------------------------------------------------------

void Drucker::getStress

  ( const Vector&  stress,
    const idx_t    mpoint ) const

{
  reduce3DVector(stress, (*latestHist_)[mpoint].sig );
}

//-----------------------------------------------------------------------
//   getStrain
//-----------------------------------------------------------------------

void Drucker::getStrain

  ( const Vector&  strain,
    const idx_t    mpoint ) const

{
  reduce3DVector(strain, (*latestHist_)[mpoint].eps );
}


//-----------------------------------------------------------------------
//   setHistory
//-----------------------------------------------------------------------

void Drucker::setHistory

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

void Drucker::commit () 

{

  System::debug() << "JH2SCon2: Commit\n";
  
  int NumApex  = 0;
  int NumPlast = 0;
  for ( int i=0; i<newHist_.size(); i++ )
  {
    if ( newHist_[i].loading == 2 )
    {
      NumApex++;
    }
    else if ( newHist_[i].loading == 1 )
    {
      NumPlast++;
    }
  }
  
  System::out() << "----------------------------------------" << endl;
  System::out() << "Drucker-Prager report:" << endl;
  System::out() << "Total number of points           = " << newHist_.size() << endl;
  System::out() << "Points under apex return         = " << NumApex << endl;
  System::out() << "Points under plastic deformation = " << NumPlast << endl;
  System::out() << "----------------------------------------" << endl;

  //System::out() << "Erik: Commit\n";
  newHist_.swap ( preHist_ );
  latestHist_ = &preHist_;
}


//-----------------------------------------------------------------------
//   Yield
//-----------------------------------------------------------------------

void Drucker::Yield

    ( const double&           epspeq,
            double&           sigY,
            double&           dsigY )
{  
  sigY  = tau_ + h_ * epspeq;
  dsigY = h_;
  
  // when full softening is reached, limit the yield strength
  if ( sigY  < limsig_ )
  {
    //System::warn() << "FULL SOFTENING...\n";
    sigY  = limsig_;
    dsigY = 0.0;
  }

}

//-----------------------------------------------------------------------
//   YieldSurf
//-----------------------------------------------------------------------

void Drucker::YieldSurf

    ( double&              sigeq,
      Vec6&                n,
      Vec6&                m,
      Tuple<double,6,6>&   dnormdsig,
      const Vec6&          sig )
      
{
  // Standard Von Mises part
  sigeq = sqrt( std::abs(3./2.*dot(sig,tmatmul(P_,sig))) );
  
  // Von Mises derivative to stress
  n = 3./2.*1./sigeq * tmatmul(P_,sig);

  // Von Mises second derivative to stress
  for ( int i=0; i<6; i++ )
  {
    for ( int j=0; j<6; j++ )
    {
      dnormdsig(i,j) = (3./(2.*sigeq)) * P_(i,j) - 1./sigeq * n[i] * n[j];
    }
  }
  
  // Adjust for zero equivalent stress situations
  if ( sigeq < 1e-14 )
  {
    //System::debug() << "NORM TO ZERO...\n";
    n = 0.0;
    dnormdsig = 0.0;
  }

  // Add Drucker-Prager parts
  double I1 = (sig[0] + sig[1] + sig[2]);
  sigeq += alpha0_ * 1./3. * I1;
  
  Vec6 dI1ds;
  dI1ds = 0.0;
  dI1ds[0] = dI1ds[1] = dI1ds[2] = 1./3.;

  // compute df/dsig and dg/dsig
  m = n + alpha1_ * dI1ds;
  n = n + alpha0_ * dI1ds;
}

//-----------------------------------------------------------------------
//   giveHistory
//-----------------------------------------------------------------------

double Drucker::giveHistory ( const idx_t ip ) const
{
  return (*latestHist_)[ip].epspeq;
}

//-----------------------------------------------------------------------
//   allocPoints
//-----------------------------------------------------------------------

void Drucker::allocPoints

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

double Drucker::giveDissipation ( const idx_t ipoint ) const
{   
  return (*latestHist_)[ipoint].dissipation;
}

//-----------------------------------------------------------------------
//   Hist_ constructor
//-----------------------------------------------------------------------

Drucker::Hist_::Hist_ () : epspeq ( 0. ), dissipation ( 0. ), sigeq (0.), p (0.)
{
  epsp = 0.0;
  eps  = 0.0;
  sig  = 0.0;
  loading = 0;
}

//-----------------------------------------------------------------------
//   Hist_ print function
//-----------------------------------------------------------------------

void Drucker::Hist_::print () const

{
  System::out() << "epsp " << epsp << ", epspeq " << epspeq << endl;
}

// -------------------------------------------------------------------
//  Hist_::toVector
// -------------------------------------------------------------------

inline void Drucker::Hist_::toVector

 ( const Vector&  vec ) const

{
  vec[0] = sig[0];
  vec[1] = sig[1];
  vec[2] = sig[2];
  vec[3] = sig[3];
  vec[4] = sig[4];
  vec[5] = sig[5];
  vec[6] = epspeq;
  vec[7] = loading;
  vec[8] = sigeq;
  vec[9] = p;
}

