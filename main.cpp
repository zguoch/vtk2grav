#include <vtkSmartPointer.h>

#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkAppendFilter.h>
#include <vtkPointDataToCellData.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

#include <iostream>
#include <fstream>
#include <string>
#include "Cube.h"
#include "MultiProgressBar.h"
#include "omp.h"
using namespace std;
vtkSmartPointer<vtkUnstructuredGrid> ReadUnstructuredGrid(std::string const& fileName, std::vector<std::string> fieldNames);
void vtkUnstructuredGrid2Cubes(std::vector<FORWARD::STRUCT_CUBE>& cubes, vtkUnstructuredGrid* usg);
void readSites_xyz(std::string fname_xyz, std::vector<FORWARD::STRUCT_SITE>& sites);

int main(int argc, char *argv[])
{
    std::cout<<" hello vtk2gm"<<std::endl;
    if(argc!=3)
    {
        cout<<"You have to give the vtu file and xyz file"<<endl;
        exit(0);
    }
    // 1. read vtu file
    std::vector<string> filedNames{"T","density"};
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = ReadUnstructuredGrid(std::string(argv[1]), filedNames);
    // 2. convert unstructuredGrid FORWARD::STRUCT_CUBE 
    std::vector<FORWARD::STRUCT_CUBE> cubes;
    vtkUnstructuredGrid2Cubes(cubes,unstructuredGrid);
    // 3. read calculation point sites
    std::vector<FORWARD::STRUCT_SITE> sites;
    readSites_xyz(argv[2], sites);
    // 4. calculate
    std::vector<double> grav;
    // initialize as 0
    for (size_t i = 0; i < sites.size(); i++)grav.push_back(0);
    FORWARD::cCube cube_forward;
    MultiProgressBar multibar(cubes.size(),COLOR_BAR_BLUE);
    int nThreads = omp_get_max_threads();
    omp_set_num_threads(nThreads);
    cout<<"Use "<<nThreads<<" threads "<<endl;
    #pragma omp parallel for shared(sites, cubes, cube_forward, grav)
    for (size_t k = 0; k < cubes.size(); k++)
    {
        for (size_t i = 0; i < sites.size(); i++)
        {
            grav[i]+=cube_forward.calSinglePoint(cubes[k],sites[i]);
        }
        #pragma omp critical
        multibar.Update();
    }
    
    // 5. write to file
    ofstream fout("gravity.txt");
    for (size_t i = 0; i < grav.size(); i++)
    {
        fout<<grav[i]<<endl;
    }
    fout.close();
    return 0;
}
void readSites_xyz(std::string fname_xyz, std::vector<FORWARD::STRUCT_SITE>& sites)
{
    sites.clear();

    std::ifstream fin(fname_xyz);
    FORWARD::STRUCT_SITE site;
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
}
void vtkUnstructuredGrid2Cubes(std::vector<FORWARD::STRUCT_CUBE>& cubes, vtkUnstructuredGrid* usg)
{
    cubes.clear();
    FORWARD::STRUCT_CUBE cube;
    vtkSmartPointer<vtkPointDataToCellData> p2c = vtkSmartPointer<vtkPointDataToCellData>::New();
    p2c->SetInputData(usg);
    p2c->PassPointDataOn();
    p2c->Update();
    vtkSmartPointer<vtkDataArray> density = p2c->GetOutput()->GetCellData()->GetArray("density");
    if(!density)
    {
        cout<<"Density doesn't exist in the vtu file"<<endl;
        exit(0);
    }
    int nCells = p2c->GetOutput()->GetNumberOfCells();
    vtkSmartPointer<vtkCell> cell;
    vtkSmartPointer<vtkPoints> points;
    int cellType;
    for (int i = 0; i < nCells; i++)
    {
        cell = usg->GetCell(i);
        cellType = cell->GetCellType();
        if(cellType==VTK_VOXEL)
        {
            cell->GetPoints();
            cube.density = density->GetTuple1(i);
            // for VTK_VOXEL(11) and VTK_HEXAHEDRON(12), regard cell as cube even though VTK_HEXAHEDRON is not a cube, but actually in ASPECT modeling results, VTK_HEXAHEDRON is cube
            // so only need the first point and last point 
            points = cell->GetPoints();
            // coordinates of the cube vertex
            cube.bound={points->GetPoint(0)[0], points->GetPoint(7)[0],
                        points->GetPoint(0)[1], points->GetPoint(7)[1],
                        points->GetPoint(0)[2], points->GetPoint(7)[2]}; //[xmin,xmax,ymin,ymax,zmin,zmax]
            cubes.push_back(cube);
        }else if(cellType==VTK_HEXAHEDRON)
        {
            cell->GetPoints();
            cube.density = density->GetTuple1(i);
            // for VTK_VOXEL(11) and VTK_HEXAHEDRON(12), regard cell as cube even though VTK_HEXAHEDRON is not a cube, but actually in ASPECT modeling results, VTK_HEXAHEDRON is cube
            // so only need the first point and last point 
            points = cell->GetPoints();
            // coordinates of the cube vertex
            cube.bound={points->GetPoint(0)[0], points->GetPoint(6)[0],
                        points->GetPoint(0)[1], points->GetPoint(6)[1],
                        points->GetPoint(0)[2], points->GetPoint(6)[2]}; //[xmin,xmax,ymin,ymax,zmin,zmax]
            cubes.push_back(cube);
        }else
        {
            cout<<"Unsupported cell type: "<<cellType<<endl;
            cout<<"The "<<i<<"th cell will be skipped!"<<endl;
        }
    }
}
vtkSmartPointer<vtkUnstructuredGrid> ReadUnstructuredGrid(std::string const& fileName, std::vector<std::string> fieldNames)
{
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
    std::string extension = "";
    if (fileName.find_last_of(".") != std::string::npos)
    {
    extension = fileName.substr(fileName.find_last_of("."));
    }
    // std::transform(extension.begin(), extension.end(),
    //                 extension.begin(), ::tolower);

    if (extension == ".vtu")
    {
    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName (fileName.c_str());
    reader->Update();
    unstructuredGrid = reader->GetOutput();
    }
    else if (extension == ".vtk")
    {
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName (fileName.c_str());
    reader->Update();
    unstructuredGrid = reader->GetOutput();
    }
    else
    {
    std::cout<<"The file extension is not recognized: "<<extension<<std::endl;
    }
    // remove not used field
    int fieldNum = unstructuredGrid->GetPointData()->GetNumberOfArrays();
    vector<string> rmFieldNames;
    for (int i = 0; i < fieldNum; i++){
        bool needRemove = true;
        string arrayName(unstructuredGrid->GetPointData()->GetArrayName(i));
        for (int j = 0; j < fieldNames.size(); j++)
        {
            if(arrayName == fieldNames[j])needRemove=false;
        }
        if(needRemove)rmFieldNames.push_back(arrayName);
    }
    for (int i = 0; i < rmFieldNames.size(); i++)
    {
        unstructuredGrid->GetPointData()->RemoveArray(rmFieldNames[i].c_str());
    }
    return unstructuredGrid;
}

