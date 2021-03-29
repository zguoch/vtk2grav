/**
 * @file Cube.h 
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Gravity forward of a Cube body with uniform density.
 * @version 0.1
 * @date 2021-03-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef Cube_H
#define Cube_H

#include<cmath>
#include<vector>
#include "stdfunc.h"

namespace FORWARD
{
    struct STRUCT_CUBE
    {
        std::vector<double> bound;
        double density;
    };
    
    class cCube
    {
    private:
        /**
         * @brief Bounding box of Cube, [xmin, xmax, ymin, ymax, zmin, zmax]
         * 
         */
        // const double* m_bound; //read only of what it pointed
        // const double m_density;
        // /**
        //  * @brief Position of a calculation point, [x0, y0, z0]
        //  * 
        //  */
        // const double* m_calSite; //read only of what it pointed
    public:
        cCube(/* args */);
        ~cCube();
        /**
         * @brief Calculate gravity of a cube at site with (x,y,z) coordinate. 
         * Upward Z axis is positive.
         * 
         * @param grav Gravity [mGal]
         * @param bound Bounding box of a cube body, [xmin, xmax, ymin, ymax, zmin, zmax] in unit [m]
         * @param density Density of the cube [kg/m3]
         * @param site Calculation point position [x0, y0, z0] in unit [m]
         * @return int 
         */
        double calSinglePoint(const STRUCT_CUBE cube, STRUCT_SITE site);
        double calSinglePoint(const std::vector<STRUCT_CUBE > cubes, STRUCT_SITE site);
    };

}

#endif