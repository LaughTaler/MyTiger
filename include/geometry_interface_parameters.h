#pragma once
#ifndef _READ_MODEL_H_
#define _READ_MODEL_H_
#include "parameters.h"
#include <sstream>

namespace TiGER {




struct readModelParameters :public APIPara
    {
    public:
      readModelParameters()
          : heal(39), getname_flag(false), bool_flag(true), divideclosed_flag(true), bool_para(1)
      {};
        DEFINE_SETTER(bool, getname_flag,readModelParameters, false);
        DEFINE_SETTER(bool, bool_flag,readModelParameters, true);
        DEFINE_SETTER(bool, divideclosed_flag, readModelParameters,true);
        DEFINE_SETTER(int, heal, readModelParameters,39);       
        DEFINE_SETTER(int, bool_para, readModelParameters, 1);       
    };


    struct generateStlParameters : public APIPara {
       public:
        generateStlParameters()
            : Deflection(0.001),
              Angle(5),
              Relative(true),
              InParallel(true),
              MinSize(1e-7),
              InternalVerticesMode(true),
              ControlSurfaceDeflection(true)    
        {};
        DEFINE_SETTER(double, Deflection, generateStlParameters, 0.001);
        DEFINE_SETTER(double, Angle, generateStlParameters, 5);
        DEFINE_SETTER(double, MinSize, generateStlParameters, 1e-7);
        DEFINE_SETTER(bool, Relative, generateStlParameters, true);
        DEFINE_SETTER(bool, InParallel, generateStlParameters, true);
        DEFINE_SETTER(bool, InternalVerticesMode, generateStlParameters, true);
        DEFINE_SETTER(bool, ControlSurfaceDeflection, generateStlParameters,
                      true);
    };

  struct generateStlDTParameters : public APIPara {
       public:
        generateStlDTParameters()
            : Deflection(0.001),
              Angle(5),
              Layer(5), sizealpha(0.7)
               , demension(10)
               , delta(0.1)
        {};
        DEFINE_SETTER(double, Deflection, generateStlDTParameters, 0.001);
        DEFINE_SETTER(double, Angle, generateStlDTParameters, 5);
        DEFINE_SETTER(int, Layer, generateStlDTParameters, 5);
        DEFINE_SETTER(double, sizealpha, generateStlDTParameters, 0.7);
        DEFINE_SETTER(int, demension, generateStlDTParameters, 10);
        DEFINE_SETTER(int, delta, generateStlDTParameters, 0.1);
    };

     struct CurveBParameters : public APIPara {
       public:
        CurveBParameters()
            : Deflection(0.001),
              Angle(5),
              Layer(5)
        {};
        DEFINE_SETTER(double, Deflection, CurveBParameters, 0.001);
        DEFINE_SETTER(double, Angle, CurveBParameters, 5);
        DEFINE_SETTER(int, Layer, CurveBParameters, 5);
    };
}

#endif