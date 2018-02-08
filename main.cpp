
#include <jive/app/declare.h>
#include <jive/app/ChainModule.h>
#include <jive/app/ControlModule.h>
#include <jive/app/OutputModule.h>
#include <jive/app/ReportModule.h>
#include <jive/app/Application.h>
#include <jive/app/InfoModule.h>
#include <jive/app/UserconfModule.h>
#include <jive/geom/declare.h>
#include <jive/model/declare.h>
#include <jive/femodel/declare.h>
#include <jive/fem/declare.h>
#include <jive/fem/InitModule.h>
#include <jive/fem/InputModule.h>
#include <jive/fem/ShapeModule.h>
#include <jive/implict/declare.h>
#include <jive/gl/declare.h>
#include <jive/model/LoadScaleModel.h>
#include <jive/implict/NonlinModule.h>
#include <jive/implict/NewmarkModule.h>
#include <jive/implict/declare.h>

#include "models.h"
#include "modules.h"
#include "TimeStepModule.h"

using namespace jem;

using jive::app::Application;
using jive::app::Module;
using jive::app::ChainModule;
using jive::app::ControlModule;
using jive::app::ReportModule;
using jive::app::UserconfModule;
using jive::app::InfoModule;
using jive::app::OutputModule;
using jive::fem::InitModule;
using jive::fem::InputModule;
using jive::fem::ShapeModule;
using jive::model::LoadScaleModel;
using jive::implict::NonlinModule;
using jive::implict::NewmarkModule;

//-----------------------------------------------------------------------
//   mainModule
//-----------------------------------------------------------------------


Ref<Module> mainModule ()
{
  Ref<ChainModule>   chain = newInstance<ChainModule> ();
  Ref<ControlModule> ctrl;

  // Register the required models and other classes.

  jive::app::declareModules   ();
  jive::fem::declareMBuilders ();
  jive::model::declareModels  ();
  jive::gl::declareModules    ();
  jive::implict::declareModels    ();
  jive::implict::declareModules   ();
  
  declareModels  ();
  declareModules ();

  // Define the main module chain.

  chain->pushBack ( newInstance<InputModule>    ( "input"    ) );
  chain->pushBack ( newInstance<ShapeModule>    ( "shape"    ) );
  chain->pushBack ( newInstance<InitModule>     ( "init"     ) );
  chain->pushBack ( newInstance<InfoModule>     ( "info"     ) );

  chain->pushBack ( newInstance<UserconfModule> ("usermodules") );

  // Wrap the entire chain in a ControlModule.
  ctrl = newInstance<ControlModule> ( "control" );

  ctrl ->runWhile ("i < 1000");
  chain->pushBack ( ctrl );  

  return newInstance<ReportModule> ( "report", chain );
}


//-----------------------------------------------------------------------
//   main
//-----------------------------------------------------------------------


int main ( int argc, char** argv )
{
  return Application::exec ( argc, argv, & mainModule );
}
