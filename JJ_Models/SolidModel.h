#ifndef SOLID_MODEL_H
#define SOLID_MODEL_H

#include <jive/algebra/MatrixBuilder.h>
#include <jive/fem/ElementGroup.h>
#include <jive/fem/ElementSet.h>
#include <jive/geom/InternalShape.h>
#include <jive/model/Model.h>
#include <jive/util/Assignable.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/XTable.h>

#include "Array.h"
#include "Material.h"
#include "utilities.h"
#include "utilitiesLarge.h"

using namespace jem;

using jem::util::Properties;
using jive::util::XTable;
using jive::util::XDofSpace;
using jive::util::Assignable;
using jive::model::Model;
using jive::geom::InternalShape;
using jive::fem::NodeSet;
using jive::fem::ElementSet;
using jive::fem::ElementGroup;
using jive::algebra::MBuilder;


//-----------------------------------------------------------------------
//   class SolidModel
//-----------------------------------------------------------------------


class SolidModel : public Model
{
 public:

  static const char*      MATERIAL_PROP;
  static const char*      RHO_PROP;
  static const char*      DOF_TYPE_NAME_1;
  static const char*      DOF_TYPE_NAME_2;
  static const char*      INTEL_PROP;
  static const char*      INTBC_PROP;
  static const char*      ELSET_PROP;
  static const char*      STRAIN_PROP;
  static const char*      OPPL_PROP;


                          SolidModel

    ( const String&         name,
      const Properties&     conf,
      const Properties&     props,
      const Properties&     globdat );

  virtual bool            takeAction

    ( const String&         action,
      const Properties&     params,
      const Properties&     globdat );


 protected:

  virtual                ~SolidModel   ();


 private:

  void                    initShape_     ();

  void                    initDofs_

    ( const Properties&     globdat );

  void                    getMatrix0_

    ( Ref<MBuilder>         mbld,
      const Vector&         fint,
      const Properties&     globdat );

  void                    getMatrix2_

    ( Ref<MBuilder>         mbld,
      const Properties&     globdat );

  bool                    getTable_

  ( const Properties&  params,
    const Properties&  globdat );

  void                    getStress_

  ( XTable&            table,
    const Vector&      weights,
    const Properties&  globdat );
    
  void                    getStressGP_

  ( XTable&            table,
    const Vector&      weights,
    const Properties&  globdat );
    
  void                    getStressStrain_

  ( const Properties&  globdat );


 private:

  Assignable<ElementGroup> egroup_;
  Assignable<NodeSet>     nodes_;
  Assignable<ElementSet>  elems_;
  Ref<InternalShape>      shape_;

  Ref<XDofSpace>          dofs_;
  IntVector               dofTypes_;

  Ref<Material>           material_;
  ShapeGradsFunc          getShapeGrads_;
  ShapeFuncsFunc          getShapeFuncs_;
  NLShapeGradsFunc        getNLShapeGrads_;
  double                  rho_;
  String                  elset_;
  
  String                  intel_;
  String                  intbc_;
  String                  mat_;
  bool                    strainForm_;
  
  int                     rank_;
  int                     ndCount_;
  int                     ipCount_;
  int                     dofCount_;
  int                     strCount_;
  String                  state_;
  
  int                     oppl_;
};


#endif
