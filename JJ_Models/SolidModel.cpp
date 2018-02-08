
#include <jem/base/Array.h>
#include <jem/base/IllegalInputException.h>
#include <jem/base/System.h>
#include <jem/numeric/algebra/LUSolver.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/io/PrintWriter.h>
#include <jem/io/FileWriter.h>

#include <jive/geom/ShapeTable.h>
#include <jive/util/DummyItemSet.h>
#include <jive/util/Globdat.h>
#include <jive/util/Printer.h>
#include <jive/geom/Quad.h>
#include <jive/geom/Triangle.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/model/ModelFactory.h>

#include "models.h"
#include "SolidModel.h"

using jem::io::PrintWriter;
using jem::io::FileWriter;
using jem::numeric::inverse;
using jem::numeric::det;
using jem::System;
using jive::util::Printer;
using jive::model::StateVector;
using jive::model::ActionParams;
using jive::util::DummyItemSet;

//=======================================================================
//   class SolidModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char*  SolidModel::MATERIAL_PROP   = "material";
const char*  SolidModel::RHO_PROP        = "rho";
const char*  SolidModel::INTEL_PROP      = "ischeme_element";
const char*  SolidModel::INTBC_PROP      = "ischeme_boundary";
const char*  SolidModel::STRAIN_PROP     = "Large_strain";
const char*  SolidModel::OPPL_PROP       = "op_plot";

const char*  SolidModel::DOF_TYPE_NAME_1 = "u";
const char*  SolidModel::DOF_TYPE_NAME_2 = "v";



//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


SolidModel::SolidModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Model ( name )

{  
 
  // Get the properties for this object. Note that myName_ is
  // protected member of the Model class.

  Properties  myConf  = conf .makeProps ( myName_ );
  Properties  myProps = props.getProps  ( myName_ );
  
  const String  context = getContext ();
  
  // Initialize the model parameters.
  
  myProps.get ( rho_, RHO_PROP );
  myConf .set ( RHO_PROP, rho_ );
  myProps.get ( intel_, INTEL_PROP );
  myConf .set ( INTEL_PROP, intel_ );
  myProps.get ( intbc_, INTBC_PROP );
  myConf .set ( INTBC_PROP, intbc_ );
  
  oppl_ = 0;
  myProps.find ( oppl_, OPPL_PROP );
  myConf .set  ( OPPL_PROP, oppl_     );
  
  strainForm_ = false;
  myProps.find ( strainForm_, STRAIN_PROP );
  myConf .set ( STRAIN_PROP, strainForm_ );
  
  // Get the elements and the nodes from the global database.
  egroup_ = ElementGroup::get ( myConf, myProps, globdat, context );
  elems_  = egroup_.getElements ();
  nodes_  = elems_.getNodes ();
  
  // Initialize the shape object and the DOFs.

  initShape_ ();
  initDofs_  ( globdat );

  getShapeGrads_    = getShapeGradsFunc        ( nodes_.rank() );
  getShapeFuncs_    = getShapeFuncsFunc        ( nodes_.rank() );
  getNLShapeGrads_  = getNLShapeGradsFunc      ( nodes_.rank() );

  // Make sure that each element has the same number of nodes as the
  // shape object.

  elems_.checkSomeElements (
    context,
    egroup_.getIndices (),
    shape_->nodeCount  ()
  );

  // Initialize the material model and allocate history variables.
  
  material_ = newMaterial ( MATERIAL_PROP, myConf, myProps, globdat );
  material_-> allocPoints ( elems_.size() * (shape_->ipointCount()) );

  material_->configure ( myProps );
  material_->getConfig ( myConf );

  // Initialize values
  
  rank_     = shape_->globalRank  ();
  ndCount_  = shape_->nodeCount   ();
  ipCount_  = shape_->ipointCount ();
  dofCount_ = dofTypes_.size()*ndCount_;
  strCount_ = STRAIN_COUNTS[rank_];
  
  // Check for proper material state
  state_    = material_->findState ();
  if ( state_ != "PLANE_STRAIN"  )
  {
    throw IllegalInputException (
      getContext (),
      "Unknown stress definition, only Plane Strain is currently allowed!"
    );
  };
  
  // Creates a dummy ItemSet to store data from integration points
  Ref<DummyItemSet>  items =
    newInstance<DummyItemSet> ( "ipoints", "ipoint" );

  items->addItems ( elems_.size() * shape_->ipointCount() );
  items->store ( globdat );
  
}


SolidModel::~SolidModel ()
{}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool SolidModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;

  if ( action == Actions::GET_INT_VECTOR )
  {
    // Assemble the internal force vector.
    System::out() << "Solid: get_int_vector\n";

    Vector  fint;

    params.get  ( fint, ActionParams::INT_VECTOR );
    getMatrix0_ ( NIL,  fint, globdat );

    return true;
  }

  if ( action == Actions::GET_MATRIX0 )
  {
    // Assemble the global stiffness matrix.
    System::out() << "Solid: get_matrix0\n";
    
    Ref<MBuilder>  mbld;
    Vector         fint;

    params.get  ( mbld, ActionParams::MATRIX0 );
    params.get  ( fint, ActionParams::INT_VECTOR );
    getMatrix0_ ( mbld, fint, globdat );

    return true;
  }

  if ( action == Actions::GET_MATRIX2 )
  {
    // Assemble the global mass matrix.
    System::out() << "Solid: get_matrix2\n";
    

    Ref<MBuilder>  mbld;

    params.get  ( mbld, ActionParams::MATRIX2 );
    getMatrix2_ ( mbld, globdat );

    return true;
  }
  
  if ( action == Actions::COMMIT )
  {
    // Commit and save history variables for next step
    System::out() << "Solid: commit\n";
    material_->commit ();
    
    return true;
  }

  if ( action == Actions::GET_TABLE )
  { 
    // Write to table for output
    System::out() << "Solid: get_table\n";
    
    return getTable_ ( params, globdat );
  }
  
  if ( action == "GET_STRESSSTRAIN" )
  { 
    // Write to table for output
    System::out() << "Solid: get_StressStrain\n";
    
    getStressStrain_( globdat );
    
    return true;//getTable_ ( params, globdat );
  } 
  
  // Unsupported action.
  return false;
}


//-----------------------------------------------------------------------
//   initShape_
//-----------------------------------------------------------------------


void SolidModel::initShape_ ()
{
  using jive::geom::Quad4;
  using jive::geom::Triangle3;
  using jive::MatrixVector;
  const int  rank     = nodes_.rank ();
  const int  maxNodes = elems_.maxElemNodeCount ();

  // Make sure that all elements have the same number of nodes.

  elems_.checkElements ( getContext(), maxNodes );

  if      ( rank == 2 && maxNodes == 3 )
  {
    //shape_ = Triangle3::getShape ( myName_ + ".shape" );
    shape_ = Triangle3::getShape ( myName_ + ".shape", intel_, intbc_ );
  }
  else if ( rank == 2 && maxNodes == 4 )
  {
    //shape_ = Quad4::getShape     ( myName_ + ".shape" );
    shape_ = Quad4::getShape     ( myName_ + ".shape", intel_ , intbc_ );
  }
  else
  {
    throw IllegalInputException (
      getContext (),
      "invalid finite element mesh"
    );
  }
}


//-----------------------------------------------------------------------
//   initDofs_
//-----------------------------------------------------------------------


void SolidModel::initDofs_ ( const Properties&  globdat )
{
  IntVector  inodes ( nodes_.size() );

  // Get the DofSpace from the global database.

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );

  // Assign the number of dof, add new DOF types and store the index.
  // The current code implements problem with two unknowns per node. 
  
  dofTypes_.resize ( 2 );

  dofTypes_[0] = dofs_->addType ( DOF_TYPE_NAME_1 );
  dofTypes_[1] = dofs_->addType ( DOF_TYPE_NAME_2 );

  // Add the node DOFs.

  inodes = iarray ( inodes.size() );

  dofs_->addDofs ( inodes, dofTypes_ );

}


//-----------------------------------------------------------------------
//   getMatrix0_
//-----------------------------------------------------------------------


void SolidModel::getMatrix0_

  ( Ref<MBuilder>      mbld,
    const Vector&      fint,
    const Properties&  globdat )

{
  using jem::numeric::matmul;
  using jive::model::StateVector;

  IntVector   ielems     = egroup_.getIndices  ();
  const int   ielemCount = ielems.size         ();
  
  Matrix      B       (strCount_, dofTypes_.size() * ndCount_ );
  Matrix      C       (strCount_,strCount_);
  Vector      stress  (strCount_);
  Vector      strain  (strCount_);
  B = 0.0;
  C = 0.0;
  stress = 0.0;
  strain = 0.0;
  
  Cubix       grads   ( rank_,    ndCount_, ipCount_ );
  Matrix      coords  ( rank_,    ndCount_ );
  Matrix      elmat   ( dofCount_, dofCount_ );
  Vector      elvec   ( dofCount_ );
  Vector      elforce ( dofCount_ );
  Vector      weights ( ipCount_ );
  IntVector   inodes  ( ndCount_ );
  IntVector   idofs   ( dofCount_ );
  
  Vector      state;
  
  idx_t ipoint = 0;
  

  // Get the state vector containing the current solution.

  StateVector::get ( state, dofs_, globdat );

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.
    int  ielem = ielems[ie];
    
    // Get the nodes and coordinates of this element.

    elems_.getElemNodes  ( inodes, ielem );
    nodes_.getSomeCoords ( coords, inodes );
    
    
    // Get the DOFs attached to this element.

    dofs_->getDofIndices ( idofs, inodes, dofTypes_ );

    // Get the shape function gradients and the integration point
    // weights.
    shape_->getShapeGradients ( grads, weights, coords );

    
    // Initialize element stiffness matrix and displacement vector.
    elmat   = 0.0;
    elforce = 0.0;
    elvec   = select ( state, idofs );
    
    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      // Construct the B-matrix for this integration point.
      getShapeGrads_ ( B, grads(ALL,ALL,ip) );
      strain = matmul( B, elvec );
       
      // find the element stiffness matrix
      material_->update ( stress, C, strain, ipoint++ );
      
      // Construct the element stiffness matrix
      elmat   += matmul(B.transpose(), matmul (C,B) ) * weights[ip];
      elforce += matmul(B.transpose(), stress) * weights[ip];
    }
    
    // build global internal force vector
    select ( fint, idofs ) += elforce;

    // Add the element matrix to the global matrix.
    if ( mbld != NIL )
    {
      mbld->addBlock ( idofs, idofs, elmat );
    }

  }
  
  
  

}

//-----------------------------------------------------------------------
//   getMatrix2_
//-----------------------------------------------------------------------


void SolidModel::getMatrix2_

  ( Ref<MBuilder>      mbld,
    const Properties&  globdat )

{
  //System::out() << "called: Mass matrix" << "\n";
  using jem::numeric::matmul;
  using jive::model::StateVector;
  
  IntVector   ielems     = egroup_.getIndices  ();
  const int   ielemCount = ielems.size         ();
  
  Matrix      N_all   ( ndCount_   , ipCount_ );
  Matrix      N       ( rank_, dofCount_   );
  N = 0.0;

  Cubix       grads   ( rank_,    ndCount_, ipCount_ );
  Matrix      coords  ( rank_,    ndCount_ );
  Matrix      elmat   ( dofCount_, dofCount_ );
  Vector      elvec   ( dofCount_ );
  Vector      weights ( ipCount_ );
  IntVector   inodes  ( ndCount_ );
  IntVector   idofs   ( dofCount_ );
  /*
  Vector      state;
  // Get the state vector containing the current solution.
  StateVector::get ( state, dofs_, globdat );
  */
  
  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.
    int  ielem = ielems[ie];
    
    // Get the nodes and coordinates of this element.

    elems_.getElemNodes  ( inodes, ielem );
    nodes_.getSomeCoords ( coords, inodes );

    // Get the DOFs attached to this element.

    dofs_->getDofIndices ( idofs, inodes, dofTypes_ );

    // Get the shape function values and the integration point
    // weights.
    N_all = shape_->getShapeFunctions ();
    shape_->getIntegrationWeights (weights, coords);

    // Assemble the element matrix.

    elmat = 0.0;

    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      // Construct the N-matrix for this integration point.
      getShapeFuncs_ ( N, N_all(ALL,ip) );
      
      // Construct the element weight matrix for this integration point.
      elmat += rho_ * matmul(N.transpose(), N ) * weights[ip];
    }

    if ( mbld != NIL )
    {
      // Add the element matrix to the global matrix.

      mbld->addBlock ( idofs, idofs, elmat );
    }
  }
}


//-----------------------------------------------------------------------
//   getTable_
//-----------------------------------------------------------------------


bool SolidModel::getTable_

  ( const Properties&  params,
    const Properties&  globdat )

{
  //System::out() << "called: getTable_\n";

  Ref<XTable>  table;
  String       name;

  // Get the name of the requested table.

  params.get ( name, ActionParams::TABLE_NAME );
  
  // Check whether the requested table is supported by this model.

  if ( name == "stress" )
  {
    params.get ( table, ActionParams::TABLE );
    
    // Check whether the table is associated with the nodes.

    if ( table->getRowItems() == elems_.getNodes().getData() )
    {
      Vector  weights;

      // Get the table weights. These will be used later to scale the
      // rows of the table.

      params.get ( weights, ActionParams::TABLE_WEIGHTS );
      getStress_ ( *table, weights, globdat );

      return true;
    }
    else if ( table->getRowItems() == elems_.getData() )
    {
      throw IllegalInputException (
        getContext (),
        "Element stress output not implemented yet, please choose 'nodes'."
        );
    }

  }
  else if ( name == "stressGP" ) 
  {
    params.get ( table, ActionParams::TABLE );

    if ( table->rowCount() == elems_.size() * shape_->ipointCount() )
    {
      // Get the table weights. These will be used later to scale the
      // rows of the table.
      Vector  weights;
      params.get ( weights, ActionParams::TABLE_WEIGHTS );
      getStressGP_ ( *table, weights, globdat );
    }
    else
    {
      throw IllegalInputException (
        getContext (),
        "Table size incorrect! Check SolidModel::getTable_"
        );
    }
    
    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   getStress_
//-----------------------------------------------------------------------


void SolidModel::getStress_

  ( XTable&            table,
    const Vector&      weights,
    const Properties&  globdat )

{
  // Initiate some variables
  IntVector   ielems     = egroup_.getIndices  ();
  const int   ielemCount = ielems.size         ();
  IntVector   inodes   ( ndCount_ );
  Matrix      elstrss  ( ndCount_, strCount_+3 );
  Vector Hist(10);
  int ipoint = 0;
  
  // Add columns to the table.
  IntVector   jcols;
  jcols.resize ( strCount_+3 );
  jcols[0] = table.addColumn ( "sigma_xx" );
  jcols[1] = table.addColumn ( "sigma_yy" );
  jcols[2] = table.addColumn ( "sigma_xy" );
  jcols[3] = table.addColumn ( "sigma_eq" );
  jcols[4] = table.addColumn ( "epspeq" );
  jcols[5] = table.addColumn ( "pressure" );
  
  // Loop over elements and integration point to obtain data. 
  for ( int ie = 0; ie < ielemCount; ie++ )
  { 
    // Get the global element index.
    int  ielem = ielems[ie];
    
    // Get the nodes of this element.
    elems_.getElemNodes  ( inodes, ielem );
    
    // zero fill element nodal stress array
    elstrss = 0.0;
    
    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      ipoint = ie * ipCount_ + ip;
      Hist = 0.0;
      material_->getHistory( Hist, ipoint );
  
      for ( int inode = 0; inode < ndCount_; inode++ )
      {
        elstrss(inode,0) += Hist[0]/ipCount_;
        elstrss(inode,1) += Hist[1]/ipCount_;
        elstrss(inode,2) += Hist[3]/ipCount_;
        elstrss(inode,3) += Hist[8]/ipCount_;
        elstrss(inode,4) += Hist[6]/ipCount_;
        elstrss(inode,5) += Hist[9]/ipCount_;
      }
  
    }
    
    // Add internal variables to table
    table.addBlock ( inodes, jcols, elstrss );
    
    // Increment the table weights associated with the element nodes.
    select ( weights, inodes ) += 1.0;
    
  }

}


//-----------------------------------------------------------------------
//   getStressGP_
//-----------------------------------------------------------------------

void SolidModel::getStressGP_

  ( XTable&            table,
    const Vector&      weights,
    const Properties&  globdat )

{
  // Initiate some variables
  IntVector   ielems     = egroup_.getIndices  ();
  const int   ielemCount = ielems.size         ();
  IntVector   inodes   ( ndCount_ );
  Matrix      coords   ( rank_, ndCount_ );
  Matrix      elstrss  ( shape_->ipointCount(),2 + 6 + 3 );
  Vector Hist(10);
  IntVector tmpipoint(shape_->ipointCount());
  tmpipoint = 0;
  idx_t ipoint = 0;
  
  Matrix ipcoords (rank_, shape_->ipointCount ());
  ipcoords = 0.0;
  
  
  // Add columns to the table.
  IntVector   jcols;
  jcols.resize ( 2 + 6 + 3 );
  
  jcols[0]  = table.addColumn ( "x-coord" );
  jcols[1]  = table.addColumn ( "y-coord" );
  jcols[2]  = table.addColumn ( "sigma_xx" );
  jcols[3]  = table.addColumn ( "sigma_yy" );
  jcols[4]  = table.addColumn ( "sigma_zz" );
  jcols[5]  = table.addColumn ( "sigma_xy" );
  jcols[6]  = table.addColumn ( "sigma_yz" );
  jcols[7]  = table.addColumn ( "sigma_xz" );
  jcols[8]  = table.addColumn ( "sigma_eq" );
  jcols[9]  = table.addColumn ( "epspeq" );
  jcols[10] = table.addColumn ( "pressure" );
  
  // Loop over elements and integration point to obtain data. 
  for ( int ie = 0; ie < ielemCount; ie++ )
  { 
    // Get the global element index.
    int  ielem = ielems[ie];
    
    // Get the nodes of this element.
    elems_.getElemNodes  ( inodes, ielem );
    nodes_.getSomeCoords ( coords, inodes );
    
    // zero fill element nodal stress array
    elstrss = 0.0;
    
    // Obtain global coordinates for integration points
    shape_->getGlobalIntegrationPoints( ipcoords, coords );
    
    for ( int ip = 0; ip < ipCount_; ip++, ipoint++ )
    {
      
      tmpipoint[ip] = ipoint;
      
      // obtain history
      Hist = 0.0;
      material_->getHistory( Hist, ipoint );
      
      // Store coordinates
      elstrss(ip,0) = ipcoords(0,ip);
      elstrss(ip,1) = ipcoords(1,ip);
      
      // Store stress variables
      elstrss(ip,2)  = Hist[0];
      elstrss(ip,3)  = Hist[1];
      elstrss(ip,4)  = Hist[2];
      elstrss(ip,5)  = Hist[3];
      elstrss(ip,6)  = Hist[4];
      elstrss(ip,7)  = Hist[5];
      
      // store internal variables
      elstrss(ip,8)  = Hist[8];
      elstrss(ip,9)  = Hist[6];
      elstrss(ip,10) = Hist[9];
    
    }
    
    // Add internal variables to table
    table.addBlock ( tmpipoint, jcols, elstrss );
    
    // Increment the table weights associated with the element nodes.
    select ( weights, tmpipoint ) += 1.0;

  }

}

//-----------------------------------------------------------------------
//   getStressStrain_
//-----------------------------------------------------------------------

void SolidModel::getStressStrain_

  ( const Properties&  globdat )

{
  System::out() << "Solid: getStressStrain_\n";
  
  // Get some stuff from Globdat
  using jive::util::Globdat;
  Properties  myVars = Globdat::getVariables ( myName_, globdat );
  
  // Initiate some variables
  IntVector   ielems     = egroup_.getIndices  ();
  const int   ielemCount = ielems.size         ();
  IntVector   inodes   ( ndCount_ );
  
  Vector      elstrss    ( 6 );
  Vector      elstrn     ( 6 );
  Vector      internal   ( 3 );
  elstrss     = 0.0;
  elstrn      = 0.0;
  internal    = 0.0;
  
  Vector Hist(10);
  Vector Hist_strain(6);
  
  int ipoint = 0;
    
  // Get step number
  int i = 0.0;
  globdat.get ( i, Globdat::TIME_STEP );
  
  // Loop over elements and integration point to obtain data. 
  for ( int ie = 0; ie < ielemCount; ie++ )
  { 
    // Get the global element index.
    int  ielem = ielems[ie];
    
    // Check if current element 
    if ( ielem != oppl_ )  
    {
      continue;
    }
    
    // Get the nodes of this element.
    elems_.getElemNodes  ( inodes, ielem );
    
    // zero fill element nodal stress array
    elstrss   = 0.0;
    elstrn    = 0.0;
    internal  = 0.0;
    
    // Loop over integration points
    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      ipoint = ie * ipCount_ + ip;
      Hist = 0.0;
      Hist_strain = 0.0;
      material_->getHistory( Hist, ipoint );
      material_->getStrain(Hist_strain, ipoint );
      
      // stresses
      elstrss[0] += Hist[0]/ipCount_;
      elstrss[1] += Hist[1]/ipCount_;
      elstrss[2] += Hist[2]/ipCount_;
      elstrss[3] += Hist[3]/ipCount_;
      elstrss[4] += Hist[4]/ipCount_;
      elstrss[5] += Hist[5]/ipCount_;
      
      // strain
      elstrn[0]  += Hist_strain[0]/ipCount_;
      elstrn[1]  += Hist_strain[1]/ipCount_;
      elstrn[2]  += Hist_strain[2]/ipCount_;//0.0;
      elstrn[3]  += Hist_strain[3]/ipCount_;
      elstrn[4]  += Hist_strain[4]/ipCount_;//0.0;
      elstrn[5]  += Hist_strain[5]/ipCount_;//0.0;
      
      // other values
      internal[0] += Hist[8]/ipCount_;
      internal[1] += Hist[6]/ipCount_;
      internal[2] += Hist[9]/ipCount_;
    }
    
    // Add internal variables to globdat. These can be stored through 
    // a SampleModule
    myVars.set ( "time" , i );
    myVars.set ( "stress"   , elstrss );
    myVars.set ( "strain"   , elstrn  );
    myVars.set ( "internal" , internal );
    
  }
  
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newSolidModel
//-----------------------------------------------------------------------


Ref<Model>            newSolidModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  // Return an instance of the SolidModel class.

  return newInstance<SolidModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareSolidModel
//-----------------------------------------------------------------------


void declareSolidModel ()
{
  using jive::model::ModelFactory;

  // Register the SolidModel with the ModelFactory.

  ModelFactory::declare ( "Solid", newSolidModel );
}
