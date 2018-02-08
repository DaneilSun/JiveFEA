#ifndef MATERIAL_H
#define MATERIAL_H

#include <jem/base/System.h>
#include <jem/util/Properties.h>
#include <jem/base/array/utilities.h>
#include <jem/base/array/intrinsics.h>
#include <jive/Array.h>
#include <jem/base/Object.h>

using jem::idx_t;
using jem::System;
using jem::Tuple;
using jem::io::endl;
using jem::Ref;
using jem::String;
using jem::Object;
using jem::util::Properties;
using jive::Vector;
using jive::Matrix;
using jive::IntVector;
using jive::BoolVector;
using jive::StringVector;

//-----------------------------------------------------------------------
//   class Material
//-----------------------------------------------------------------------

// The Material class represents a material model. Its main task is to
// calculate the strain-stress stiffness matrix.
// At least 'update' should be implemented by any derived class

class Material : public Object
{
 public:


  // General material model functions
  
  explicit                Material

    ( const idx_t           rank,
      const Properties&     globdat );

  virtual void            configure

    ( const Properties&     props );

  virtual void            getConfig

    ( const Properties&     conf )            const;

  virtual void            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      const idx_t           ip ) = 0;

  virtual void            commit ();

  virtual void            cancel ();

  virtual bool            despair ();

  virtual bool            endDespair ();
  
  
  // History related parameters

  virtual void            allocPoints

    ( const idx_t           count );

  virtual void            allocPoints

    ( const idx_t           count,
      const Matrix&         transfer,
      const IntVector&      oldPoints );

  virtual void            deallocPoints 
    
    ( const idx_t           count );

  idx_t                   getHistoryCount   () const;

  StringVector            getHistoryNames   () const;

  void                    getHistoryNames   
    
    ( const StringVector&   hnames ) const;

  virtual void            getHistory

    ( const Vector&         hvals,
      const idx_t           mpoint ) const;
      
  virtual void            getStress

    ( const Vector&         stress,
      const idx_t           mpoint ) const;
      
  virtual void            getStrain

    ( const Vector&         strain,
      const idx_t           mpoint ) const;

  virtual double          giveHistory       ( const idx_t point  ) const;
  
  
  // Internal state variable related functions
  
  inline virtual idx_t    isLoading         ( idx_t point  ) const; 
  inline virtual idx_t    wasLoading        ( idx_t point  ) const;
  inline virtual idx_t    pointCount  ()    const;
  inline virtual double   giveDissipation   ( const idx_t point  ) const;
  inline virtual double   giveTheta   ()    const;
  inline virtual double   givePoisson ()    const;
  inline virtual String   findState   ()    const;
  
  
  bool                    isViscous ()      const;

  void                    setDT

    ( const double          dt );

  inline virtual Vector   getThermStrain  () const;

  virtual void            updateThermStrain

    ( const Properties&     globdat );
    
  
  // Stress related functions
  
  /*  
  virtual Matrix          getStiffMat

    ( )                                      const = 0;
    */
  virtual Tuple<double,6> fill3DStrain

    ( const Vector&         v3 )             const = 0;

  virtual Tuple<double,6> fill3DStress

    ( const Vector&         v3 )             const = 0;

 protected:

  virtual                ~Material      ();

  Tuple<double,6>         fillFrom2D_

    ( const Vector&         v3,
      const double          zzval ) const;

  Tuple<double,6>         fillFrom3D_

    ( const Vector&         v6 ) const;

  void                    reduceTo2D_

    ( const Vector&          v3,
      const double           zz,
      const Tuple<double,6>& v6 )            const;

  void                    reduceTo3D_

    ( const Vector&          v3,
      const Tuple<double,6>& v6 )            const;


  // Variables //////////////////////////////////////////////////

 protected:

  idx_t                  rank_;

  Vector                 thermStrain_;
  bool                   desperateMode_;
  BoolVector             hasSwitched_;
  BoolVector             useSecant_;

  bool                   viscous_;
  double                 dt_;

  idx_t                  id_;
  static idx_t           counter_;

  StringVector           historyNames_;
};

//-----------------------------------------------------------------------
//   implementation
//-----------------------------------------------------------------------
// The functions below are the default implementations. The actual material
// models should redefine the functions if they are used. Otherwise zero, false
// or NAN may be the default output for the function.

inline double Material::giveTheta () const
{
  // default implementation: no angle

  return NAN;
}

inline double Material::giveHistory ( const idx_t ipoint ) const
{   
  // default implementation: no history
  return 0.;
}

inline double Material::giveDissipation ( const idx_t ipoint ) const
{   
  // default implementation: no dissipation
  return 0.;
}

inline idx_t Material::pointCount () const
{
  return 0;
}

inline idx_t Material::isLoading ( idx_t ipoint ) const
{
  return false;
}

inline idx_t Material::wasLoading ( idx_t ipoint ) const
{
  return false;
}

inline Vector Material::getThermStrain () const
{
  return thermStrain_;
}

inline bool Material::isViscous () const
{
  return viscous_;
}

inline void Material::setDT ( const double dt ) 
{
  dt_ = dt;
}

inline double Material::givePoisson () const
{
  return 0.;
}

inline String Material::findState () const
{
  return "None!\n";
}

//-----------------------------------------------------------------------
//   related functions
//-----------------------------------------------------------------------


Ref<Material>             newMaterial

  ( const String&           name,
    const Properties&       conf,
    const Properties&       props,
    const Properties&       globdat );

#endif 
