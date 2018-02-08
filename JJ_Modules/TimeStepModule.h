#ifndef TIME_STEP_MODULE_H
#define TIME_STEP_MODULE_H

#include <jive/app/Module.h>
#include <jive/util/DofSpace.h>

#include "import.h"
#include "Array.h"


//-----------------------------------------------------------------------
//   class TimeStepModule
//-----------------------------------------------------------------------


class TimeStepModule : public Module
{
 public:

  typedef Module            Super;
  typedef TimeStepModule    Self;


  explicit                  TimeStepModule

    ( const String&           name = "step" );

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

  virtual                  ~TimeStepModule  ();


 private:

  void                      restart_

    ( const Properties&       globdat );

  void                      invalidate_     ();


 private:

  Ref<Model>                model_;
  Ref<DofSpace>             dofs_;
  Vector                    rmass_;
  double                    dtime_;
  int                       istep_;
  int                       nodeID_;
  int                       inode_;

};


#endif
