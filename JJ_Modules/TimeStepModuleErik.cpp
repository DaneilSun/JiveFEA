/*
 * 
 *  Copyright (C) 2016 TU Delft. All rights reserved.
 *  
 *  This module implements an explicit solver for JemJive. A central difference
 *  scheme is used. The current Module allows the user to prescribe displacement
 *  velocity or acceleration as a boundary condition. It is also possible to 
 *  release
 *
 *  Author : E.C.Simons
 *  Date   : December 2016
 *
 */

#include <jem/base/System.h>
#include <jem/base/limits.h>
#include <jem/base/Array.h>
#include <jem/base/ArithmeticException.h>
#include <jem/util/Event.h>
#include <jem/util/Properties.h>
#include <jive/util/Globdat.h>
#include <jive/util/ItemSet.h>
#include <jive/util/DofSpace.h>
#include <jive/util/Constraints.h>
#include <jive/algebra/DiagMatrixObject.h>
#include <jive/algebra/LumpedMatrixBuilder.h>
#include <jive/algebra/SparseMatrixBuilder.h>
#include <jive/model/Model.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/app/ModuleFactory.h>

#include "Names.h"
#include "TimeStepModuleErik.h"


using namespace jem;

using jive::util::Globdat;
using jive::model::Actions;
using jive::model::ActionParams;
using jive::model::StateVector;



//=======================================================================
//   class TimeStepModuleErik
//=======================================================================

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


TimeStepModuleErik::TimeStepModuleErik ( const String& name ) :

  Super ( name )

{
  dtime_  =  1.0;
  istep_  = -1;
  BCtype_ = 0; // default = displacement boundary conditions
}


TimeStepModuleErik::~TimeStepModuleErik ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void TimeStepModuleErik::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps = props.findProps ( myName_ );

  myProps.find ( dtime_,  PropNames::DTIME, 1.0e-20, 1.0e20 );
  myProps.find ( BCtype_, "BCtype", 0, 2 );

  myProps.find ( TermGroups_ , "TermGroups" );
  myProps.find ( TermDofs_   , "TermDofs"   );
  myProps.find ( TermSteps_  , "TermSteps"  );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void TimeStepModuleErik::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );

  myConf.set ( PropNames::DTIME, dtime_ );
  myConf.set ( "BCtype",  BCtype_ );

  myConf.set ( "TermGroups", TermGroups_ );
  myConf.set ( "TermDofs"  , TermDofs_   );
  myConf.set ( "TermSteps" , TermSteps_  );
}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status TimeStepModuleErik::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  using jem::util::connect;
  
  Properties  myProps = props.findProps ( myName_ );
  const String  context = getContext ();

  double        t;
  int           i;

  // returning an address of an object pointing to class ElasticModel
  
  model_ = Model   ::get ( globdat, context );
  

  // returning an address of an object pointing to class DofSpace
  
  dofs_  = DofSpace::get ( globdat, context );
  nodes_ = NodeSet::find ( globdat );

  istep_ = -1;
  

  // Invalidate the state of this module when the DofSpace changes.

  connect ( dofs_->newSizeEvent,  this, &Self::invalidate_ );
  connect ( dofs_->newOrderEvent, this, &Self::invalidate_ );


  // Initialize the global simulation time and the time step number.

  i = 0;
  t = 0.0;

  if ( ! globdat.find( t, Globdat::TIME ) )
  {
    globdat.set ( Globdat::TIME, t );
  }

  if ( ! globdat.find( i, Globdat::TIME_STEP ) )
  {
    globdat.set ( Globdat::TIME_STEP, i );
  }

  globdat.set ( Globdat::OLD_TIME,      t );
  globdat.set ( Globdat::OLD_TIME_STEP, i );

  System::warn() << "TimeStepModuleErik: Not sure if the boundary types are implemented correctly... Please check first before using intensively.\n";

  return OK;
}


//-----------------------------------------------------------------------
//   run
//-----------------------------------------------------------------------


Module::Status TimeStepModuleErik::run ( const Properties& globdat )
{
  System::out() << "TimeStepModuleErik: run called!\n";
  using jive::util::Constraints;
  using jive::util::setSlaveDofs;
  using jive::model::STATE0;
  using jive::model::STATE1;
  using jive::model::STATE2;
  using jive::algebra::SparseMatrixBuilder;

  if ( model_ == NIL )
  {
    return DONE;
  }

  if ( istep_ < 0 )
  {
    restart_ ( globdat );
  }

  //System::out() << "Time run start\n";

  Ref<Constraints>  cons;
  Ref<SparseMatrixBuilder>  mbuilder;

  const int     dofCount = rmass_.size ();
  const double  dt2      = dtime_ * dtime_;

  Vector        fint ( dofCount );
  Vector        fext ( dofCount );

  Vector        u    ( dofCount );
  Vector        u2   ( dofCount );
  Vector        v    ( dofCount );
  Vector        a    ( dofCount );
  Vector        u0, u1, v0, v1, a0, a1;

  Properties    params;

  double        t;
  int           i;

  // Get the internal and external force vectors.

  fint = 0.0;
  fext = 0.0;

  params.set ( ActionParams::INT_VECTOR, fint );
  params.set ( PropertyNames::DTIME, dtime_ );
  
  model_->takeAction ( Actions::ADVANCE, params, globdat );

  model_->takeAction ( Actions::GET_INT_VECTOR, params, globdat );
 

  //System::out()<<fint<<"\n";

  params.set ( ActionParams::EXT_VECTOR, fext );

  model_->takeAction ( Actions::GET_EXT_VECTOR, params, globdat );

  // Get the state vectors and compute the state vector at the next
  // time step.

  StateVector::get    ( u0, STATE0, dofs_, globdat );
  StateVector::getOld ( u1, STATE0, dofs_, globdat );
  StateVector::get    ( v0, STATE1, dofs_, globdat );
  StateVector::get    ( a0, STATE2, dofs_, globdat );

  u  = (fext - fint) * dt2 * rmass_ + 2.0 * u0 - u1;
  //System::out() << u1 << "\n" << u0 << "\n" << u << "\n";
  
  u2 = u1;
  u1 = u0;
  u0 = u;

  // Update the simulation time and time step.

  globdat.get ( t, Globdat::TIME );
  globdat.get ( i, Globdat::TIME_STEP );

  globdat.set ( Globdat::OLD_TIME,      t );
  globdat.set ( Globdat::OLD_TIME_STEP, i );
  globdat.set ( Globdat::TIME,          t + dtime_ );
  globdat.set ( Globdat::TIME_STEP,     i + 1 );

  istep_++;

  model_->takeAction ( Actions::COMMIT, params, globdat );

  System::out() << "on time: " << t << " in step: " << i << "\n\n";

  

  // Set the prescribed values if necessary.
  cons = Constraints::find ( dofs_, globdat );
  params.set ( ActionParams::CONSTRAINTS, cons );
  
  // set the constraints given by the models. This includes LoadScaleModel etc.
  model_->takeAction ( Actions::GET_CONSTRAINTS, params, globdat );

  
  // Check if any constraints need to be removed for the current time step.
  Assignable<NodeGroup> group;
  IdxVector             inodes;
  
  for ( idx_t ig = 0; ig < TermSteps_.size(); ig++ )
  {
  
      if ( i >= TermSteps_[ig] )
      {
      // find nodeGroups of which the BC has to be deleted at some point
      group  = NodeGroup::get ( TermGroups_[ig], nodes_, globdat, getContext() );

      // find node indices
      idx_t nn = group.size();
      inodes . resize ( nn );
      inodes = group.getIndices ();
    
      // find dof types
      idx_t itype  = dofs_->findType ( TermDofs_[ig] );
    
      for ( idx_t in = 0; in < nn; in++ )
      {
        // find dof index and remove from the constraints
        idx_t idof = dofs_->getDofIndex ( inodes[in], itype );
        cons->eraseConstraint(idof);
      }
    }
  }

  
  // Prescribe values on the remaining constraints
  if ( cons != NIL )
  {
  
    if ( BCtype_ == 0 )
    {
      // applied displacement
      setSlaveDofs ( u0, *cons );
    }
    else if ( BCtype_ == 1 )
    {
      // applied velocity
      
      // compute current guess for velocity
      v = 1./dtime_ * ( u0 - u1 );
      setSlaveDofs ( v, *cons );
    
      // compute constrained displacement for new time step
      u0 = v * dtime_ + u1;
    }
    else if ( BCtype_ == 2 )
    {
      // applied acceleration
      
      // compute current guess for acceleration
      a = 1./dt2 * ( u0 - 2*u1 + u2 );
      setSlaveDofs ( a, *cons );
      
      // compute constrained displacement for new time step
      u0 = dt2 * a + 2 * u1 - u2;
    }
    else
    {
      throw IllegalInputException (
        getContext (),
        "Unknown prescribed boundary type! Please chose displacement (0), velocity (1) or acceleration (2)."
        );
    }
  }

  
  // compute actual acceleration and velocity
  v = 1./dtime_ * ( u0 - u1 );
  a = 1./dt2 * ( u0 - 2*u1 + u2 );
  
  v0 = v;
  a0 = a;

  return OK;
}


//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void TimeStepModuleErik::shutdown ( const Properties& globdat )
{
  model_ = NIL;
}


//-----------------------------------------------------------------------
//   restart_
//-----------------------------------------------------------------------


void TimeStepModuleErik::restart_ ( const Properties& globdat )
{
  System::out() << "TimeStepModuleErik: restart called!\n";
  
  
  
  
  
  using jive::algebra::LumpedMatrixBuilder;
  using jive::model::STATE0;
  using jive::model::STATE1;
  using jive::model::STATE2;

  Ref<LumpedMatrixBuilder>  mbuilder;

  const String  context  = getContext     ();
  const int     dofCount = dofs_->dofCount ();

  Properties    params;

  Vector        u0, u1, a0, v0;

  Vector        fext(dofCount);

  fext = 0.0;

  // Assemble the lumped mass matrix.

  mbuilder = newInstance<LumpedMatrixBuilder> ();

  mbuilder->setOptions ( 0 );
  mbuilder->setSize    ( dofCount );
  mbuilder->clear      ();

 // params.set , the mbuilder is stored or set under the name MATRIX2.

  params.set ( ActionParams::MATRIX2, mbuilder );

  model_->takeAction ( Actions::GET_MATRIX2, params, globdat );

  mbuilder->updateMatrix ();

  // Store and invert the lumped mass matrix.

  rmass_.resize ( dofCount );

  rmass_ = mbuilder->getDiagMatrix()->getValues ();

  for ( int i = 0; i < dofCount; i++ )
  {
    if ( isTiny( rmass_[i] ) )
    {
      throw ArithmeticException (
        context,
        "singular mass matrix"
      );
    }

    rmass_[i] = 1.0 / rmass_[i];
  }


  params.set ( ActionParams::EXT_VECTOR, fext );

  model_->takeAction ( Actions::GET_EXT_VECTOR, params, globdat );

  // Compute the old state vector.

  StateVector::get    ( u0, STATE0, dofs_, globdat );
  StateVector::getOld ( u1, STATE0, dofs_, globdat );
  StateVector::get    ( v0, STATE1, dofs_, globdat );
  StateVector::get    ( a0, STATE2, dofs_, globdat );


  a0 = fext*rmass_;

  u1 = u0 - dtime_ * v0 + (0.5 * dtime_ * dtime_) * a0;
}


//-----------------------------------------------------------------------
//   invalidate_
//-----------------------------------------------------------------------


void TimeStepModuleErik::invalidate_ ()
{
  istep_ = -1;
}

//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------

Ref<Module>  TimeStepModuleErik::makeNew

  ( const String&           name,
    const Properties&       conf,
    const Properties&       props,
    const Properties&       globdat)

{
  return newInstance<Self> ( name );
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   declareFlexArclenModule
//-----------------------------------------------------------------------

void declareTimeStepModuleErik ()
{
  using jive::app::ModuleFactory;

  ModuleFactory::declare ( "stepErik",
                         & TimeStepModuleErik::makeNew );
}
