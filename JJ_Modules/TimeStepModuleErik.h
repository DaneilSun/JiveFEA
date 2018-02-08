#ifndef TIME_STEP_MODULE_ERIK_H
#define TIME_STEP_MODULE_ERIK_H

#include <jive/app/Module.h>
#include <jive/fem/ElementGroup.h>
#include <jive/fem/ElementSet.h>
#include <jive/fem/NodeGroup.h>
#include <jive/util/Assignable.h>
#include <jive/util/DofSpace.h>

#include "import.h"
#include "Array.h"

using jive::fem::ElementGroup;
using jive::fem::ElementSet;
using jive::fem::NodeGroup;
using jive::fem::NodeSet;
using jive::util::Assignable;
using jive::StringVector;

//-----------------------------------------------------------------------
//   class TimeStepModuleErik
//-----------------------------------------------------------------------


class TimeStepModuleErik : public Module
{
 public:

  typedef Module            Super;
  typedef TimeStepModuleErik    Self;


  explicit                  TimeStepModuleErik

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

  virtual                  ~TimeStepModuleErik  ();


 private:

  void                      restart_

    ( const Properties&       globdat );

  void                      invalidate_     ();


 private:

  
  Ref<Model>                model_;
  Ref<DofSpace>             dofs_;
  Assignable<NodeSet>       nodes_;
  
  Vector                    rmass_;
  double                    dtime_;
  int                       istep_;
  
  int                       BCtype_; // 0 = displacement, 1 = velocity, 2 = acceleration
  
  StringVector              TermGroups_;
  StringVector              TermDofs_;
  IntVector                 TermSteps_;
  
};


#endif
