/*
 * 
 *  Copyright (C) 2016 TU Delft. All rights reserved.
 *  
 *  This module implements an explicit solver for JemJive. A central difference
 *  scheme is used.
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
#include "TimeStepModule.h"


using namespace jem;

using jive::util::Globdat;
using jive::model::Actions;
using jive::model::ActionParams;
using jive::model::StateVector;



//=======================================================================
//   class TimeStepModule
//=======================================================================

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


TimeStepModule::TimeStepModule ( const String& name ) :

  Super ( name )

{
  dtime_  =  1.0;
  istep_  = -1;
  nodeID_ =  0;
  inode_  = -1;
}


TimeStepModule::~TimeStepModule ()
{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status TimeStepModule::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  using jem::util::connect;

  const String  context = getContext ();

  double        t;
  int           i;

  // returning an address of an object pointing to class ElasticModel
  model_ = Model   ::get ( globdat, context );

  // returning an address of an object pointing to class DofSpace
  dofs_  = DofSpace::get ( globdat, context );

  istep_ = -1;

  // Invalidate the state of this module when the DofSpace changes.

  connect ( dofs_->newSizeEvent,  this, &Self::invalidate_ );
  connect ( dofs_->newOrderEvent, this, &Self::invalidate_ );

  // Get the index of the node to be monitored.

  inode_ = dofs_->getItems()->findItem ( nodeID_ );

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

  return OK;
}


//-----------------------------------------------------------------------
//   run
//-----------------------------------------------------------------------


Module::Status TimeStepModule::run ( const Properties& globdat )
{
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
  Vector        u0, u1, v0, v1, a0, a1;

  Properties    params;

  double        t;
  int           i;

  // Get the internal and external force vectors.

  fint = 0.0;
  fext = 0.0;

  params.set ( ActionParams::INT_VECTOR, fint );

  model_->takeAction ( Actions::GET_INT_VECTOR, params, globdat );

  //System::out()<<fint<<"\n";

  params.set ( ActionParams::EXT_VECTOR, fext );

  model_->takeAction ( Actions::GET_EXT_VECTOR, params, globdat );

  // Get the state vectors and compute the state vector at the next
  // time step.

  StateVector::get    ( u0, STATE0, dofs_, globdat );
  StateVector::getOld ( u1, STATE0, dofs_, globdat );

  u  = (fext - fint) * dt2 * rmass_ + 2.0 * u0 - u1;
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

  if ( cons != NIL )
  {
    params.set ( ActionParams::CONSTRAINTS, cons );

    model_->takeAction ( Actions::GET_CONSTRAINTS, params, globdat );

    setSlaveDofs ( u0, *cons );
  }

  // Save the displacements in the selected node if it exists.

  if ( inode_ >= 0 )
  {
    const int   typeCount = dofs_->typeCount      ();
    Properties  myVars    = Globdat::getVariables ( myName_, globdat );

    for ( i = 0; i < typeCount; i++ )
    {
      String  typeName = dofs_->getTypeName ( i );
      int     idof     = dofs_->getDofIndex ( inode_, i );

      myVars.set ( typeName, u[idof] );
    }
  }
  
  //System::out() << "Time run end\n";
  return OK;
}


//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void TimeStepModule::shutdown ( const Properties& globdat )
{
  model_ = NIL;
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void TimeStepModule::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps = props.findProps ( myName_ );

  myProps.find ( dtime_,  PropNames::DTIME, 1.0e-20, 1.0e20 );
  myProps.find ( nodeID_, PropNames::NODE );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void TimeStepModule::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );

  myConf.set ( PropNames::DTIME, dtime_ );
  myConf.set ( PropNames::NODE,  nodeID_ );
}


//-----------------------------------------------------------------------
//   restart_
//-----------------------------------------------------------------------


void TimeStepModule::restart_ ( const Properties& globdat )
{
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


void TimeStepModule::invalidate_ ()
{
  istep_ = -1;
}

//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------

Ref<Module>  TimeStepModule::makeNew

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

void declareTimeStepModule ()
{
  using jive::app::ModuleFactory;

  ModuleFactory::declare ( "step",
                         & TimeStepModule::makeNew );
}
