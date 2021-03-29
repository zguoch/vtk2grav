#ifndef STDAFX_H
#define STDAFX_H
namespace FORWARD
{
    //========================== constants =====================================
    // If other parameter use SI unit, the gravity result will in mGal unit
    #define G 6.6743E-6
    #define PI 3.1415926535897932384626433832795
    #define EPS 1E-10
    //==========================           =====================================

    struct STRUCT_SITE
    {
        double x,y,z;
    };
    
}


#endif