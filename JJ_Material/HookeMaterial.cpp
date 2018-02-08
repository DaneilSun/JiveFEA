/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens the isotropic elastic material
 *  This represents the material at a point in space.
 *  It is implemented in such a way that can be used for any
 *  discretisation methods, say finite elements, EFG and so on.
 *  
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *
 */
#include <jem/base/array/utilities.h>
#include <jem/base/Array.h>

#include <jem/base/limits.h>
#include <jem/base/Error.h>
#include <jem/base/PrecheckException.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/LUSolver.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/utilities.h>
#include <jem/base/System.h>

#include "utilities.h"
#include "HookeMaterial.h"

using namespace jem;
using jem::numeric::dotProduct;
using jem::numeric::LUSolver;
using jem::numeric::matmul;

const double one_third = 0.3333333333333333;

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  HookeMaterial::YOUNG_PROP   = "young";
const char*  HookeMaterial::POISSON_PROP = "poisson";
const char*  HookeMaterial::RHO_PROP     = "rho";
const char*  HookeMaterial::STATE_PROP   = "state";

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


HookeMaterial::HookeMaterial 

  ( idx_t rank_, const Properties& globdat )
    : Material ( rank_, globdat )
{
  JEM_PRECHECK ( rank_ >= 1 && rank_ <= 3 );

  young_   = 1.0;
  poisson_ = 1.0;
  rho_     = 1.0;

  stiffMat_ .resize ( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  stiffMat_ = 0.0;
}


HookeMaterial::~HookeMaterial ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void HookeMaterial::configure ( const Properties& props )
{
  props.get ( young_,   YOUNG_PROP,   0.0, maxOf( young_ ) );
  props.get ( poisson_, POISSON_PROP, 0.0, 0.5 );

  if ( rank_ == 1 )
  {
    props.get ( area_, "area" );
  }

  props.find ( rho_, RHO_PROP, 0.0, maxOf( rho_ ) );

  // read problem type, plane stress ...

  if ( rank_ == 2  )
  {
    props.get( stateString_, STATE_PROP );
    //System::out() << "Hooke " << stateString_ << "\n";
    if      ( stateString_ == "PLANE_STRAIN" )
    {
      state_ = PlaneStrain;
    }
    else if ( stateString_ == "PLANE_STRESS" )
    {
      state_ = PlaneStress;
    }
    else if ( stateString_ == "AXISYMMETRIC" )
    {
      state_ = AxiSymmetric;
    }
  }
  // compute the elastic moduli, only once time

  computeStiffMat_ ();
  
  // Build P and Q matrix, used in Mises and MisesE
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
//   getConfig
//-----------------------------------------------------------------------


void HookeMaterial::getConfig ( const Properties& conf ) const
{
  conf.set ( YOUNG_PROP,   young_   );
  conf.set ( POISSON_PROP, poisson_ );

  if ( rank_ == 2 )
  {
    conf.set ( STATE_PROP, stateString_ );
  }
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void HookeMaterial::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      idx_t                 ipoint )

{
  
  // compute the elastic moduli

  stiff = stiffMat_;

  // compute the stress vector

  matmul ( stress, stiffMat_, strain );
  
  // Store history variables
  // This step was added to the Hooke material model. The stress level can be
  // computed simply from strain and stiffness, for output reasons the values 
  // are however stored, so they do not have to be recomputed. This particular
  // material model does not gain in speed by this, but more elaborate models
  // do. To keep models general, also the current model stores history variables.
  Tuple<double,6> sig, eps;
  sig = fill3DStress(stress);
  eps = fill3DStrain(strain);

  newHist_[ipoint].sig   = sig;
  newHist_[ipoint].eps   = eps;
  newHist_[ipoint].sigeq = Mises(sig);
  newHist_[ipoint].p     = pressure(sig);
  
  latestHist_ = &newHist_;
}

//-----------------------------------------------------------------------
//   getStiffMat
//-----------------------------------------------------------------------


Matrix HookeMaterial::getStiffMat() const
{
  return stiffMat_;
}
  
//-----------------------------------------------------------------------
//   computeStiffMat_
//-----------------------------------------------------------------------


void   HookeMaterial::computeStiffMat_ () 
{
  const idx_t   n  = STRAIN_COUNTS[rank_];

  const double  e  = young_;
  const double  nu = poisson_;


  if      ( rank_ == 1 )
  {
    stiffMat_(0,0) = e * area_;
  }
  else if ( rank_ == 3 )
  {
    const double  a = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double  b = 0.5 * (1.0 - 2.0 * nu);
    const double  c = 1.0 - nu;

    stiffMat_(0,0) = a * c;
    stiffMat_(0,1) = stiffMat_(0,2) = a * nu;
    stiffMat_(1,1) = stiffMat_(0,0);
    stiffMat_(1,2) = stiffMat_(0,1);
    stiffMat_(2,2) = stiffMat_(0,0);
    stiffMat_(3,3) = a * b;
    stiffMat_(4,4) = stiffMat_(3,3);
    stiffMat_(5,5) = stiffMat_(3,3);

    // Copy lower triangle of the stress-strain matrix.

    for ( idx_t i = 0; i < n; i++ )
    {
      for ( idx_t j = 0; j < i; j++ )
      {
        stiffMat_(i,j) = stiffMat_(j,i);
      }
    }
  }
  else if ( state_ == PlaneStrain )
  {
    const double  a = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double  b = 0.5 * (1.0 - 2.0 * nu);
    const double  c = 1.0 - nu;

    stiffMat_(0,0) = a * c;
    stiffMat_(0,1) = stiffMat_(1,0) = a * nu;
    stiffMat_(1,1) = a * c;
    stiffMat_(2,2) = a * b;
  }
  else if ( state_ == PlaneStress )
  {
    const double  a = e / (1.0 - nu * nu);
    stiffMat_(0,0) = a;
    stiffMat_(0,1) = stiffMat_(1,0) = a * nu;
    stiffMat_(1,1) = a;
    stiffMat_(2,2) = a * 0.5 * (1.0 - nu);
  }
  else if ( state_ == AxiSymmetric )
  {
    stiffMat_ .resize ( STRAIN_COUNTS[rank_]+1, STRAIN_COUNTS[rank_]+1 );
    stiffMat_ = 0.0;
  
    const double  a = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double  b = 0.5 * (1.0 - 2.0 * nu);
    const double  c = 1.0 - nu;

    stiffMat_(0,0) = a * c;
    stiffMat_(0,1) = stiffMat_(1,0) = a * nu;
    stiffMat_(1,1) = a * c;
    stiffMat_(2,2) = a * b;
    stiffMat_(3,3) = a * c;
    stiffMat_(3,0) = stiffMat_(0,3) = a * nu;
    stiffMat_(3,1) = stiffMat_(1,3) = a * nu;
  }
  else
  {
    throw Error ( JEM_FUNC, "unexpected rank: " + String ( rank_ ) );
  }

}

//-----------------------------------------------------------------------
//   fill3DStress
//-----------------------------------------------------------------------

Tuple<double,6> HookeMaterial::fill3DStress

  ( const Vector&    v3 ) const

{
  if ( v3.size() == 3 )
  {
    double sig_zz = state_ == PlaneStress
                  ? 0.
                  : poisson_ * ( v3[0] + v3[1] );
    
    return fillFrom2D_ ( v3, sig_zz );
  }
  else if ( v3.size() == 4 )
  {
    // axisymmetric case
    double sig_zz = v3[3];
    return fillFrom2D_ ( v3, sig_zz );
  }
  else
  {
    return fillFrom3D_ ( v3 );
  }
}

//-----------------------------------------------------------------------
//   fill3DStrain
//-----------------------------------------------------------------------

Tuple<double,6> HookeMaterial::fill3DStrain

  ( const Vector&    v3 ) const

{
  

  if ( v3.size() == 3 )
  {
    double eps_zz = state_ == PlaneStress
                  ? -poisson_ / (1.-poisson_) * (v3[0]+v3[1])
                  : 0.;

    return fillFrom2D_ ( v3, eps_zz );
  }
  else if ( v3.size() == 4 )
  {
    // axisymmetric case
    double eps_zz = v3[3];
    return fillFrom2D_ ( v3, eps_zz );
  }
  else
  {
    return fillFrom3D_ ( v3 );
  }
}

//-----------------------------------------------------------------------
//   reduce3DStrainGrad
//-----------------------------------------------------------------------

void HookeMaterial::reduce3DStrainGrad

  ( const Vector&          v3,
    const Tuple<double,6>& v6 ) const 

{
  // reduce a full 3D tuple with gradients wrt strain
  // to a vector that is either 2D or 3D

  throw Error ( JEM_FUNC, "Untested!" );

  if ( v3.size() == 3 )
  {
    double zzContrib = 0.;

    if ( state_ == PlaneStress )
    {
      zzContrib = -poisson_ / (1.-poisson_) * v6[2];
    }
    reduceTo2D_ ( v3, zzContrib, v6 );
  }
  else
  {
    reduceTo3D_ ( v3, v6 );
  }
}

//-----------------------------------------------------------------------
//   reduce3DStressGrad
//-----------------------------------------------------------------------

void HookeMaterial::reduce3DStressGrad

  ( const Vector&          v3,
    const Tuple<double,6>& v6 ) const 

{
  // reduce a full 3D tuple with gradients wrt stress
  // to a vector that is either 2D or 3D

  throw Error ( JEM_FUNC, "Untested!" );

  if ( v3.size() == 3 )
  {
    double zzContrib = 0.;

    if ( state_ == PlaneStrain )
    {
      zzContrib = -poisson_ / (1.-poisson_) * v6[2];
    }
    reduceTo2D_ ( v3, zzContrib, v6 );
  }
  else
  {
    reduceTo3D_ ( v3, v6 );
  }
}

//-----------------------------------------------------------------------
//   reduce3DVector
//-----------------------------------------------------------------------

void HookeMaterial::reduce3DVector

  ( const Vector&          v,
    const Tuple<double,6>& t ) const 

{
  // reduce a full 3D tuple to a 2D or 3D vector
  // (works the same for both stress and strain)

  if ( v.size() == 3 )
  {
    v[0] = t[0];
    v[1] = t[1];
    v[2] = t[3];
  }
  else if ( v.size() == 4 )
  {
    v[0] = t[0];
    v[1] = t[1];
    v[2] = t[3];
    v[3] = t[2];
  }
  else
  {
    v[0] = t[0];
    v[1] = t[1];
    v[2] = t[2];
    v[3] = t[3];
    v[4] = t[4];
    v[5] = t[5];
  }
}

//-----------------------------------------------------------------------
//   reduce3DMatrix
//-----------------------------------------------------------------------

void HookeMaterial::reduce3DMatrix

  ( const Matrix&            m,
    const Tuple<double,6,6>& t ) const 

{
  if ( m.size(0) == 3 )
  {
    if ( state_ == PlaneStrain )
    {
      select2DMatrix ( m, t );
    }
    else
    {
      double d;
      Tuple<double,6,6> tmp66 = t;
      LUSolver::invert ( tmp66, d );
      select2DMatrix ( m, tmp66 );
      LUSolver::invert ( m, d );
    }
  }
  else if ( m.size(0) == 4 )
  {
    //throw Error ( JEM_FUNC, "Untested!" );
    select2DMatrix ( m, t );
    
    for ( int i=0; i<2; i++ )
    {
      m(3,i) = t(i,2);
      m(i,3) = t(2,i);
    }
    m(3,3) = t(2,2);
  }
  else
  {
    for ( idx_t i = 0; i < 6; ++i )
    {
      for ( idx_t j = 0; j < 6; ++j )
      {
        m(i,j) = t(i,j);
      }
    }
  }
}

//-----------------------------------------------------------------------
//   select2DMatrix
//-----------------------------------------------------------------------

void HookeMaterial::select2DMatrix

  ( const Matrix&            m,
    const Tuple<double,6,6>& t ) const 

{
  // selecting [xx, yy, xy] components from [xx, yy, zz, xy, yz, xz]

  for ( idx_t i = 0; i < 2; ++i )
  {
    m(i,2) = t(i,3);
    m(2,i) = t(3,i); 

    for ( idx_t j = 0; j < 2; ++j )
    {
      m(i,j) = t(i,j);
    }
  }
  m(2,2) = t(3,3);
}

//-----------------------------------------------------------------------
//   compute Von Mises equivalent stress
//-----------------------------------------------------------------------

double HookeMaterial::Mises

    ( const Tuple<double,6>& t6 )
{
  double sigeq = 0.0;
  Vector v6(6);
  reduce3DVector(v6,t6);
  
  sigeq = sqrt( std::abs(3./2.*dot(v6,matmul(P_,v6))) );

  return sigeq;
}

//-----------------------------------------------------------------------
//   compute Von Mises equivalent strain
//-----------------------------------------------------------------------

double HookeMaterial::MisesE

    ( const Tuple<double,6>& t6 )
{
  double epseq = 0.0;
  Vector v6(6);
  reduce3DVector(v6,t6);
  
  epseq = sqrt( std::abs(2./3.*dot(v6,matmul(Q_,v6))) );

  return epseq;
}

//-----------------------------------------------------------------------
//   split into pressure
//-----------------------------------------------------------------------

double HookeMaterial::pressure

  ( const Tuple<double,6>& t6 ) const

{
  return -(t6[0] + t6[1] + t6[2])*1./3.;
}

//-----------------------------------------------------------------------
//   commit
//-----------------------------------------------------------------------

void HookeMaterial::commit () 

{
  System::out() << "Hooke: Commit\n";
  newHist_.swap ( preHist_ );
  latestHist_ = &preHist_;
}

//-----------------------------------------------------------------------
//   Hist_ constructor
//-----------------------------------------------------------------------

HookeMaterial::Hist_::Hist_ () 
{

  sig   = 0.;
  eps   = 0.;
  sigeq = 0.;
  p     = 0.;

}

//-----------------------------------------------------------------------
//   getHistory
//-----------------------------------------------------------------------

void HookeMaterial::getHistory

  ( const Vector&  hvals,
    const idx_t    mpoint ) const

{
  (*latestHist_)[mpoint].toVector ( hvals );
}

//-----------------------------------------------------------------------
//   getStress
//-----------------------------------------------------------------------

void HookeMaterial::getStress

  ( const Vector&  stress,
    const idx_t    mpoint ) const

{
  reduce3DVector(stress, (*latestHist_)[mpoint].sig );
}

//-----------------------------------------------------------------------
//   getStrain
//-----------------------------------------------------------------------

void HookeMaterial::getStrain

  ( const Vector&  strain,
    const idx_t    mpoint ) const

{
  reduce3DVector(strain, (*latestHist_)[mpoint].eps );
}

// -------------------------------------------------------------------
//  Hist_::toVector
// -------------------------------------------------------------------

inline void HookeMaterial::Hist_::toVector

 ( const Vector&  vec ) const

{
  vec[0] = sig[0];
  vec[1] = sig[1];
  vec[2] = sig[2];
  vec[3] = sig[3];
  vec[4] = sig[4];
  vec[5] = sig[5];
  vec[6] = 0.0;//epspeq;
  vec[7] = 0.0;//dissipation;
  vec[8] = sigeq;
  vec[9] = p;

}

//-----------------------------------------------------------------------
//   allocPoints
//-----------------------------------------------------------------------

void HookeMaterial::allocPoints

  ( idx_t count )

{
  //System::out() << "Erik: allocPoints\n";

  for ( idx_t i = 0; i < count; ++i )
  {
    preHist_.pushBack ( Hist_() );
    newHist_.pushBack ( Hist_() );
  }
}
