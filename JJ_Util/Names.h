#ifndef NAMES_H
#define NAMES_H


//-----------------------------------------------------------------------
//   class Names
//-----------------------------------------------------------------------


class Names
{
 public:

  static const char*    DOFS[3];

};


//-----------------------------------------------------------------------
//   class PropertyNames
//-----------------------------------------------------------------------


class PropertyNames
{
 public:

  static const char*    DTIME;
  static const char*    EXPR;
  static const char*    NODE;
  static const char*    YOUNG_PROP;
  static const char*    POISSON_PROP;
  static const char*    RHO;
  static const char*    SAMPLE_ELEM;
  static const char*    SHAPE;
  static const char*    THICKNESS;
  static const char*    DT_ASSEMBLY;
};


typedef PropertyNames   PropNames;



#endif
