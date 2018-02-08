/*
 * 
 *  Copyright (C) 2010 TU Delft. All rights reserved.
 *  
 *  This class implements a model for dirichlet boundary conditions.
 *
 *  Author:  F.P. van der Meer, F.P.vanderMeer@tudelft.nl
 *  Date:    May 2010
 *
 */

#include <jem/base/System.h>
#include <jem/base/Float.h>
#include <jem/numeric/utilities.h>
#include <jem/io/PrintWriter.h>
#include <jem/io/FileWriter.h>
#include <jive/fem/NodeGroup.h>
#include <jive/model/Actions.h>
#include <jive/model/ModelFactory.h>
#include <jive/util/XDofSpace.h>

#include "DirichletModel.h"
#include "SolverNames.h"

using jem::io::endl;


using jem::io::PrintWriter;
using jem::io::FileWriter;
using jive::fem::NodeGroup;
using jive::model::Actions;
using jive::util::XDofSpace;


//=======================================================================
//   class DirichletModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  DirichletModel::TYPE_NAME       = "Dirichlet";

const char*  DirichletModel::DISP_INCR_PROP  = "dispIncr";
const char*  DirichletModel::DISP_RATE_PROP  = "dispRate";
const char*  DirichletModel::INIT_DISP_PROP  = "initDisp";
const char*  DirichletModel::MAX_DISP_PROP   = "maxDisp";
const char*  DirichletModel::NODES_PROP      = "nodeGroups";
const char*  DirichletModel::DOF_PROP        = "dofs";
const char*  DirichletModel::FACTORS_PROP    = "factors";
const char*  DirichletModel::LOADED_PROP     = "loaded";

const idx_t  DirichletModel::U_LOAD_         = 1 << 0;

// Added by Erik:
const char*  DirichletModel::TURN_DISP_PROP   = "turnDisp";
const char*  DirichletModel::HOLD_DISP_PROP   = "holdDisp";

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


DirichletModel::DirichletModel

  ( const String&      name,
    const Ref<Model>&  child ) :

    Super  ( name  )

{
  dispScale0_  = 0.;
  dispIncr0_   = 0.0;
  initDisp_    = 0.0;
  maxDispVal_  = Float::MAX_VALUE;
  turnDispVal_ = Float::MAX_VALUE;
  holdDispVal_ = Float::MAX_VALUE;
  turnBool_    = false;
  unloadBool_  = false;
  method_      = INCREMENT;

}


DirichletModel::~DirichletModel ()
{}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool DirichletModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  // initialization

  if ( action == Actions::INIT )
  {
    init_ ( globdat );

    return true;
  }

  // apply displacement increment

  if ( action == Actions::GET_CONSTRAINTS )
  {
    applyConstraints_ ( params, globdat );

    return true;
  }

  // check state

  if ( action == SolverNames::CHECK_COMMIT || action == Actions::CHECK_COMMIT )
  {
    System::debug() << "DIRI: CHECK_COMMIT\n";
    checkCommit_ ( params, globdat );

    return true;
  }

  // proceed to next time step

  if ( action == Actions::COMMIT )
  {
    commit_ ( params, globdat );

    return true;
  }

  // advance to next time step

  if ( action == Actions::ADVANCE )
  {
    globdat.set ( "var.accepted", true );

    advance_ ( globdat );

    return true;
  }

  // adapt step size

  else if ( action == SolverNames::SET_STEP_SIZE )
  {
    System::debug() << "DIRI: SET_STEP_SIZE\n";
    setDT_ ( params );
  }
  
  if ( action == "WRITE_RESTART" )
  {
    // then write new files
    Ref<PrintWriter> oFile = newInstance<PrintWriter> (newInstance<FileWriter> ( "Diri.restart" ));
    oFile->nformat.setScientific(true);
    oFile->nformat.setFractionDigits (15);
    oFile->setPageWidth(10000);

    *oFile << dispScale_ << endl;
    
    oFile->close();
    
    system("mkdir NewRestart");
    system("cp *.restart NewRestart");
    
    system("mkdir Restart");
    system("mv *.restart Restart/");
    
  }
  
  
  return false;
}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void DirichletModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps = props.findProps ( myName_ );

  double maxD = Float::MAX_VALUE;

  if ( myProps.find ( dispRate_, DISP_RATE_PROP ) )
  {
    myProps.find ( initDisp_,  INIT_DISP_PROP );

    method_     = RATE;
    dispScale0_ = dispScale_ = 0.;
  }
  else
  {
    myProps.find ( initDisp_,  INIT_DISP_PROP );
    
    //
    int out = system("ls Diri.restart");
    if ( out == 0 )
    {
      Ref<FileReader> iFile = newInstance<FileReader> ( "Diri.restart" );
      initDisp_ = iFile->parseFloat();
    }
    //
    
    
    myProps.get ( dispIncr0_, DISP_INCR_PROP );

    method_     = INCREMENT;
    dispIncr_   = dispIncr0_;
    dispScale0_ = dispScale_ = initDisp_ - dispIncr0_;  // cancel first increment
  }

  myProps.find ( maxDispVal_ , MAX_DISP_PROP , 0.0, maxD );
  myProps.find ( turnDispVal_, TURN_DISP_PROP, 0.0, maxD );
  myProps.find ( holdDispVal_, HOLD_DISP_PROP, 0.0, maxD );
  
  System::debug() << "ERIK " << maxDispVal_ << "\n";
  System::debug() << "ERIK " << turnDispVal_ << "\n";

  myProps.get( nodeGroups_, NODES_PROP );
  ngroups_ = nodeGroups_.size ( );

  myProps.get( dofTypes_, DOF_PROP );

  if ( dofTypes_.size() != ngroups_ )
  {
    throw IllegalInputException ( JEM_FUNC,
          "dofTypes must have the same length as nodeGroups" );
  }

  if ( myProps.find ( factors_, FACTORS_PROP ) )
  { 
    if ( factors_.size() != ngroups_ )
    {
      throw IllegalInputException ( JEM_FUNC,
            "dofTypes must have the same length as nodeGroups" );
    }
  }
  else
  {

    idx_t loaded;

    factors_.resize ( ngroups_ );

    factors_ = 0.;

    if ( myProps.find( loaded, LOADED_PROP, -1, ngroups_-1 ) )
    {
      factors_[loaded] = 1.;
    }

  }
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void DirichletModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );

  if ( method_ == INCREMENT )
  {
    myConf.set ( DISP_INCR_PROP, dispIncr0_  );
    myConf.set ( INIT_DISP_PROP, initDisp_   );
  }
  else
  {
    myConf.set ( DISP_RATE_PROP, dispRate_   );
  }

  myConf.set ( MAX_DISP_PROP ,  maxDispVal_  );
  myConf.set ( TURN_DISP_PROP,  turnDispVal_ );
  myConf.set ( HOLD_DISP_PROP,  holdDispVal_ );

  myConf.set ( NODES_PROP,     nodeGroups_   );
  myConf.set ( DOF_PROP,       dofTypes_     );
  myConf.set ( FACTORS_PROP,   factors_      );
}



//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------


Ref<Model> DirichletModel::makeNew

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  return newInstance<Self> ( name );
}

//-----------------------------------------------------------------------
//   init_
//-----------------------------------------------------------------------


void DirichletModel::init_ ( const Properties& globdat )
{
  // Get nodes, then dofs of nodes, and constraints of dofs

  nodes_ = NodeSet::find    ( globdat );
  dofs_  = XDofSpace::get   ( nodes_.getData(), globdat );
  cons_  = Constraints::get ( dofs_, globdat );
  
    
}

//-----------------------------------------------------------------------
//   advance_
//-----------------------------------------------------------------------

void DirichletModel::advance_

  ( const Properties&  globdat )

{
  if ( method_ == RATE && std::abs(dispIncr_) < Float::EPSILON )
  {
    System::warn() << myName_ << " zero increment in RATE mode."
      << " It seems the time increment has not been set." << endl;
  }

  dispScale_   = dispScale0_ + dispIncr_;
  
  // Added by Erik
  if ( ( std::abs(dispScale_) > turnDispVal_ ) && ( turnBool_ == false ) ) 
  {
    System::warn() << "Turning Point reached, flipping dispIncr_!!!" << endl;
    dispIncr_ *= -1.;
    turnBool_ = true;
  }
  else if ( std::abs(dispScale_) > holdDispVal_ )
  {
    System::warn() << "Holding Point reached, stopped increasing dispIncr_!!!" << endl;
    dispIncr_ *= 0.;
    turnBool_ = true;
  }
  else if (( std::abs(dispScale_*dispScale0_) != dispScale_*dispScale0_ ) &&  turnBool_ )
  {
    System::warn() << "Full unloading after flip, Flag for termination!!!" << dispIncr_ << endl;
    unloadBool_ = true;
  }
  // End added by Erik

  System::debug() << "New displacement factor " << dispScale_ << endl;
}

//-----------------------------------------------------------------------
//   applyConstraints_
//-----------------------------------------------------------------------

void DirichletModel::applyConstraints_

  ( const Properties&  params,
    const Properties&  globdat )

{
  idx_t                 nn;
  Assignable<NodeGroup> group;
  IdxVector             inodes;
  String                context;

  // loop over node groups

  for ( idx_t ig = 0; ig < ngroups_; ig++ )
  {
    group  = NodeGroup::get ( nodeGroups_[ig], nodes_, globdat, context );

    nn     = group.size();
    inodes . resize ( nn );
    inodes = group.getIndices ();

    idx_t itype  = dofs_->findType ( dofTypes_[ig] );

    double val = dispScale_ * factors_[ig];

    // apply constraint

    for ( idx_t in = 0; in < nn; in++ )
    {
      idx_t idof = dofs_->getDofIndex ( inodes[in], itype );
      cons_->addConstraint ( idof, val );
    }
  }

  // compress for more efficient storage

  cons_->compress();
}

//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------

void DirichletModel::checkCommit_

  ( const Properties&  params,
    const Properties&  globdat )

{
  // terminate the computation if displacement exceeds maximum.
  // be careful with this!
  
  System::debug() << "DIRI, maxDispVal_ = " << maxDispVal_ << ", dispScale_ = " << dispScale_ << endl;

  if ( std::abs ( dispScale_ ) > maxDispVal_ ) 
  {
    System::warn() << myName_ << " says: TERMINATE because "
      << " disp > maxDispVal." << endl;

    params.set ( SolverNames::TERMINATE, "sure" );
  }
  else if ( unloadBool_ ) 
  {
    System::warn() << myName_ << " says: TERMINATE because "
      << " loading/unloading cycle is completed." << endl;

    params.set ( SolverNames::TERMINATE, "sure" );
  }
}

//-----------------------------------------------------------------------
//   commit_
//-----------------------------------------------------------------------


void DirichletModel::commit_

  ( const Properties&  params,
    const Properties&  globdat )

{
  // store converged boundary quantities

  dispScale0_  = dispScale_;
}

//-----------------------------------------------------------------------
//   setDT_
//-----------------------------------------------------------------------

void DirichletModel::setDT_

  ( const Properties&  params )

{
  double       dt;
  double       dt0;

  params.get ( dt,  SolverNames::STEP_SIZE   );
  params.get ( dt0, SolverNames::STEP_SIZE_0 );

  if ( method_ == RATE )
  {
    dispIncr_ = dispRate_ * dt;
  }
  else
  {
    // Added By Erik
    double chk0, chk1, chk2;
    chk0 = chk1 = chk2 = 1.0;
    chk0 = jem::numeric::sign(chk2,dispIncr0_);
    chk1 = jem::numeric::sign(chk2,dispIncr_);
    (chk0 != chk1 ) ? dispIncr0_ *= -1.0 : 0.0;
    // End Added By Erik
  
    // rate dependent with an initial increment is also supported
    // the displacement rate is then implicit input as dispIncr0/dt0
    dispIncr_  = dispIncr0_ * dt / dt0;
    
    System::out() << "Setting step size to " << dispIncr_ << " ("
    << dt / dt0 * 100 << "\% of initial step size)" << endl;
  }
  
  // Erik 2016-10-14: There is a bug in the program when used in combination
  // with the AdaptiveStepModule. Both in the current Model and the other 
  // Module can a step size be defined. The actual initial step size used by 
  // the program is the one defined in the DirichletModel. The step size from 
  // the AdaptiveStep is merely used to scale the step size by the Model...
  
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   declareDirichletModel
//-----------------------------------------------------------------------


void declareDirichletModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( DirichletModel::TYPE_NAME,
                          & DirichletModel::makeNew );
}
