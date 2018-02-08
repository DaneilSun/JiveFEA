//=======================================================================
//
// Model that prints load-displacement data for nodegroup to file:
// The sum of the nodal loads and the average of the nodal displacements
//
// Adapted from: MonitorModel, 
// FPM, 29-1-2008
// 
// Data is computed in GET_MATRIX0 and GET_INT_VECTOR
//         and optionally written to file in COMMIT
//         (rather using the SampleModule to write to file recommended)
//
//=======================================================================

#include <jem/io/PrintWriter.h>
#include <jem/io/FileWriter.h>
#include <jem/base/array/select.h>
#include <jem/base/array/utilities.h>
#include <jem/base/System.h>
#include <jem/util/Properties.h>
#include <jive/util/Globdat.h>
#include <jive/util/ItemSet.h>
#include <jive/util/DofSpace.h>
#include <jive/util/Assignable.h>
#include <jive/model/Actions.h>
#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/StateVector.h>
#include <jive/fem/ElementSet.h>
#include <jive/fem/NodeGroup.h>

#include "models.h"


using namespace jem;
using jem::io::endl;

using jem::util::Properties;
using jem::io::PrintWriter;
using jem::io::FileWriter;
using jive::Vector;
using jive::IntVector;
using jive::IntMatrix;
using jive::StringVector;
using jive::util::DofSpace;
using jive::util::Assignable;
using jive::model::Model;
using jive::fem::NodeSet;
using jive::fem::NodeGroup;
using jive::fem::ElementSet;


//=======================================================================
//   class LoadDispModel
//=======================================================================

class LoadDispModel : public Model
{
 public:

  static const char*        NODES_PROP;
  static const char*        GROUP_PROP;
  static const char*        TYPES_PROP;

                            LoadDispModel

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       conf,
      const Properties&       globdat )      const;

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat );

  void                      updateIDofs_

    ();


 protected:

  virtual                  ~LoadDispModel  ();


 private:

  Ref<DofSpace>             dofs_;
  Assignable<ElementSet>    elems_;
  Assignable<NodeSet>       nodes_;

  IntVector                 inodes_;
  IntMatrix                 idofs_;
  int                       nn_;
  int                       rank_;
  String                    groupName_;
  Assignable<NodeGroup>     ngroup_;
  int                       typeCount_;
  IntVector                 types_;
  StringVector              typeNames_;
};


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  LoadDispModel::NODES_PROP     = "nodes";
const char*  LoadDispModel::GROUP_PROP     = "group";
const char*  LoadDispModel::TYPES_PROP     = "types";


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


LoadDispModel::LoadDispModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Model ( name )

{
  const String       context = getContext ();

  elems_ = ElementSet::get ( globdat, context );
  nodes_ = elems_.getNodes ();
  dofs_  = DofSpace::get ( nodes_.getData(), globdat, context );
  nn_ = -1;
  rank_ = nodes_.rank();
  typeCount_ = dofs_->typeCount();
  groupName_ = "";
}


LoadDispModel::~LoadDispModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void LoadDispModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  using jive::util::Globdat;

  Properties    myProps = props.findProps ( myName_ );

  const String  context = getContext ();

  IntVector     nodeIDs;

  Properties    myVars = Globdat::getVariables ( myName_, globdat );

  if ( myProps.find( nodeIDs, NODES_PROP ) )
  {
    // Get node and dof numbers specified with 'nodes = IntVector'
    
    nn_ = nodeIDs.size();

    inodes_.resize ( nn_ );
    idofs_.resize ( nn_ , types_.size() );

    for ( int i = 0; i < nn_; i++ )
    {
      inodes_[i] = nodes_.findNode ( nodeIDs[i] );

      for ( int j = 0; j < types_.size(); ++j )
      {
        idofs_(i,j) = dofs_->getDofIndex ( inodes_[i], types_[j] );
      }
    }
  }

  if ( myProps.find( groupName_, GROUP_PROP ) )
  {
    // Get node and dof numbers specified with 'group = String'
    
    ngroup_ = NodeGroup::get ( groupName_, nodes_, globdat, context );
  }
  
  if ( myProps.find( typeNames_, TYPES_PROP ) )
  {
    types_.resize ( typeNames_.size() );
  }
  else
  {
    types_.resize ( rank_ );
    typeNames_.resize ( rank_ );

    for ( int i = 0; i < rank_; ++i )
    {
      typeNames_[i] = dofs_->getTypeName ( i );
    }
  }

  for ( int i = 0; i < typeNames_.size(); ++i )
  {
    types_[i] = dofs_->getTypeIndex ( typeNames_[i] );
  }
}

//-----------------------------------------------------------------------
//   updateIDofs_
//-----------------------------------------------------------------------


void LoadDispModel::updateIDofs_

  ()

{
  nn_ = ngroup_.size();

  inodes_.resize ( nn_ );
  idofs_ .resize ( nn_ , types_.size() );

  inodes_ = ngroup_.getIndices();

  for ( int i = 0; i < nn_; i++ )
  {
    for ( int j = 0; j < types_.size(); ++j )
    {
      idofs_(i,j) = dofs_->findDofIndex ( inodes_[i], types_[j] );
    }
  }
}
  

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void LoadDispModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );

  myConf.set ( NODES_PROP, inodes_ );

  myConf.set ( GROUP_PROP, groupName_ );

  myConf.set ( TYPES_PROP, typeNames_ );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool LoadDispModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::util::Globdat;
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  if ( nn_ < 0 && nodes_ == NIL )
  {
    return false;
  }

  Vector      load ( types_.size() ) ;
  Vector      disp ( types_.size() ) ;

  if ( action == Actions::GET_MATRIX0 ||
       action == Actions::GET_INT_VECTOR )
  {
    System::out() << "LoadDisp: get_matrix0 or get_int_vector\n";
    // get global data
    
    Properties  myVars = Globdat::getVariables ( myName_, globdat );

    // get internal force and solution vector
    
    Vector      fint;
    Vector      state;

    params.get ( fint, ActionParams::INT_VECTOR );
    StateVector::get ( state, dofs_, globdat );
    updateIDofs_();


    // compute cumulative load and average displacement
    for ( int i = 0; i < types_.size(); i++ )
    {
      load[i] = sum ( select ( fint , idofs_(ALL,i) ) );
      disp[i] = sum ( select ( state , idofs_(ALL,i) ) ) / nn_;
    }

    // store in globdat

    myVars.set ( "load" , load );
    myVars.set ( "disp" , disp );

    return true;
  }
  return false;
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newLoadDispModel
//-----------------------------------------------------------------------


static Ref<Model>     newLoadDispModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<LoadDispModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareLoadDispModel
//-----------------------------------------------------------------------


void declareLoadDispModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "LoadDisp", & newLoadDispModel );
}
