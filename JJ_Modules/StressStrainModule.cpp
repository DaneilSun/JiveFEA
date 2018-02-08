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
#include "StressStrainModule.h"


using namespace jem;

using jive::util::Globdat;
using jive::model::Actions;
using jive::model::ActionParams;
using jive::model::StateVector;



//=======================================================================
//   class StressStrainModule
//=======================================================================

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


StressStrainModule::StressStrainModule ( const String& name ) :

  Super ( name )

{}


StressStrainModule::~StressStrainModule ()
{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status StressStrainModule::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  using jem::util::connect;

  const String  context = getContext ();

  // returning an address of an object pointing to class ElasticModel
  model_ = Model   ::get ( globdat, context );
  
  return OK;
}


//-----------------------------------------------------------------------
//   run
//-----------------------------------------------------------------------


Module::Status StressStrainModule::run ( const Properties& globdat )
{
  if ( model_ == NIL )
  {
    return DONE;
  }
  
  Properties  params;
  params.clear();
  model_->takeAction ( "GET_STRESSSTRAIN", params, globdat );
  
  return OK;
  
}


//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void StressStrainModule::shutdown ( const Properties& globdat )
{
  model_ = NIL;
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void StressStrainModule::configure

  ( const Properties&  props,
    const Properties&  globdat )

{}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void StressStrainModule::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{}

//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------

Ref<Module>  StressStrainModule::makeNew

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

void declareStressStrainModule ()
{
  using jive::app::ModuleFactory;

  ModuleFactory::declare ( "ss",
                         & StressStrainModule::makeNew );
}

