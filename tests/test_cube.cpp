#include "Cube.h"
#include <iostream>
#include <fstream>

void test_Cube_Gravity_SinglePoint();
void test_Cube_Gravity_MultiPoints(std::string fname_xyz);

int main()
{
    std::cout<<"=========== Test cube gravity ==============="<<std::endl;
    std::cout<<" -- test single point"<<std::endl;
    test_Cube_Gravity_SinglePoint();
    std::cout<<" -- test multi-points"<<std::endl;
    test_Cube_Gravity_MultiPoints("sites.xyz");
    std::cout<<"=========== Test cube gravity end ==========="<<std::endl;
}
void test_Cube_Gravity_MultiPoints(std::string fname_xyz)
{
    std::ifstream fin(fname_xyz);
    FORWARD::STRUCT_SITE site;
    std::vector<FORWARD::STRUCT_SITE> sites;
    if(!fin)
    {
        std::cout<<"Open file failed: "<<fname_xyz<<std::endl;
        return;
    }
    while(!fin.eof())
    {
        fin>>site.x>>site.y>>site.z;
        sites.push_back(site);
    }
    // calculate
    std::vector<double> grav;
    FORWARD::cCube cube_forward;
    FORWARD::STRUCT_CUBE cube;
    std::vector<double> bound{-2.5, 2.5, -2.5, 2.5, -6, -3};
    double density = 1000;
    cube.bound = bound;
    cube.density=density;
    for (size_t i = 0; i < sites.size(); i++)
    {
        grav.push_back(cube_forward.calSinglePoint(cube,sites[i]));
    }
    for (size_t i = 0; i < grav.size(); i++)
    {
        std::cout<<sites[i].x<<" "<<sites[i].y<<" "<<sites[i].z<<" "<<grav[i]<<std::endl;
    }
    
}
void test_Cube_Gravity_SinglePoint()
{
    FORWARD::cCube cube_foward;
    FORWARD::STRUCT_CUBE cube;
    std::vector<double> bound{-2.5, 2.5, -2.5, 2.5, -6, -3};
    double density = 1000;
    FORWARD::STRUCT_SITE site={1,-2,3};
    cube.bound = bound;
    cube.density=density;
    double grav = cube_foward.calSinglePoint(cube, site);

    // info
    std::cout<<"Cube bound (m): ";
    for (int i = 0; i < 6; i++)std::cout<<bound[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Cube density (kg/m3): "<<density<<std::endl;
    std::cout<<"Calculation point (m): ";
    std::cout<<site.x<<" "<<site.y<<" "<<site.z<<std::endl;
    std::cout<<"Cube gravity (mGal): "<<grav<<std::endl;
}