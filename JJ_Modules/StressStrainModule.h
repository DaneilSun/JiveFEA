#ifndef STRESSSTRAIN_MODULE_H
#define STRESSSTRAIN_MODULE_H

#include <jive/app/Module.h>
#include <jive/util/DofSpace.h>

#include "import.h"
#include "Array.h"


//-----------------------------------------------------------------------
//   class StressStrainModule
//-----------------------------------------------------------------------


class StressStrainModule : public Module
{
 public:

  typedef Module            Super;
  typedef StressStrainModule    Self;


  explicit                  StressStrainModule

    ( const String&           name = "ss" );

  virtual Status            init

    ( const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual Status            run

    ( const Properties&       globdat );

  virtual void              shutdown

    ( const Properties&       globdat );

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       props,
      const Properties&       globdat )        const;
      
  static Ref<Module>       makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat);


 protected:

  virtual                  ~StressStrainModule  ();

 private:

  Ref<Model>                model_;

};


#endif
