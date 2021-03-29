#include "Cube.h"
namespace FORWARD
{

    cCube::cCube(/* args */)
    // :m_bound(nullptr),
    // m_density(0),
    // m_calSite(nullptr)
    {
        
    }

    cCube::~cCube()
    {
    }
    
    double cCube::calSinglePoint(const STRUCT_CUBE cube, STRUCT_SITE site)
    {
        double x1=cube.bound[0],x2=cube.bound[1];
        double y1=cube.bound[2],y2=cube.bound[3];
        double zmin=cube.bound[4], zmax=cube.bound[5];
        // tmp var
        double l_top2site,l_bot2site;
        double tempx1,tempx2,tempy1,tempy2,G122,G121,G112,G111,G222,G221,G212,G211;
        double p222,p221,p212,p211,p122,p121,p112,p111;	
        double p2,p1;
        double x,y;
        l_top2site=fabs(site.z - cube.bound[5]); //distance from top of cube to calculate point
        l_bot2site=fabs(site.z - cube.bound[4]);//distance from bottom of cube to calculate point
        p2=l_bot2site*l_bot2site;
        p1=l_top2site*l_top2site;
        y=site.y;
        x=site.x;
        tempx1=x1-x;
        tempx2=x2-x;

        tempy1=y1-y;
        tempy2=y2-y;
    
        p222=sqrt(pow(tempx2,2)+pow(tempy2,2)+p2);
        p221=sqrt(pow(tempx2,2)+pow(tempy2,2)+p1);
        p212=sqrt(pow(tempx2,2)+pow(tempy1,2)+p2);
        p211=sqrt(pow(tempx2,2)+pow(tempy1,2)+p1);

        p122=sqrt(pow(tempx1,2)+pow(tempy2,2)+p2);
        p121=sqrt(pow(tempx1,2)+pow(tempy2,2)+p1);
        p112=sqrt(pow(tempx1,2)+pow(tempy1,2)+p2);
        p111=sqrt(pow(tempx1,2)+pow(tempy1,2)+p1);

        G222=tempx2*log(tempy2+p222)+tempy2*log(tempx2+p222)+l_bot2site*atan2((l_bot2site*p222),(tempx2*tempy2));
        G221=tempx2*log(tempy2+p221)+tempy2*log(tempx2+p221)+l_top2site*atan2((l_top2site*p221),(tempx2*tempy2));
        G212=tempx2*log(tempy1+p212)+tempy1*log(tempx2+p212)+l_bot2site*atan2((l_bot2site*p212),(tempx2*tempy1));
        G211=tempx2*log(tempy1+p211)+tempy1*log(tempx2+p211)+l_top2site*atan2((l_top2site*p211),(tempx2*tempy1));

        G122=tempx1*log(tempy2+p122)+tempy2*log(tempx1+p122)+l_bot2site*atan2((l_bot2site*p122),(tempx1*tempy2));
        G121=tempx1*log(tempy2+p121)+tempy2*log(tempx1+p121)+l_top2site*atan2((l_top2site*p121),(tempx1*tempy2));
        G112=tempx1*log(tempy1+p112)+tempy1*log(tempx1+p112)+l_bot2site*atan2((l_bot2site*p112),(tempx1*tempy1));
        G111=tempx1*log(tempy1+p111)+tempy1*log(tempx1+p111)+l_top2site*atan2((l_top2site*p111),(tempx1*tempy1));
        
        return -G*cube.density*(G222+G211+G121+G112-G221-G212-G122-G111); //mGal
    }

    double cCube::calSinglePoint(const std::vector<STRUCT_CUBE > cubes, STRUCT_SITE site)
    {
        double grav = 0;
        for (size_t i = 0; i < cubes.size(); i++)
        {
            grav += calSinglePoint(cubes[i], site);
        }
        return grav;
    }
}