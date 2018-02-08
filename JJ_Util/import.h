#ifndef IMPORT_H
#define IMPORT_H

//-----------------------------------------------------------------------
//   forward declarations
//-----------------------------------------------------------------------


namespace jem
{
  namespace util
  {
    class Properties;
  }
}

namespace jive
{
  namespace util
  {
    class DofSpace;
  }

  namespace app
  {
    class Module;
  }

  namespace model
  {
    class Model;
  }
}


//-----------------------------------------------------------------------
//   using statements/declarations
//-----------------------------------------------------------------------


using namespace jem;
using jem::util::Properties;
using jive::app::Module;
using jive::model::Model;
using jive::util::DofSpace;


#endif

