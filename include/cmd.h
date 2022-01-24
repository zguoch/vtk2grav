#ifndef CMD_H
#define CMD_H
#include "version.h"
#include "getopt.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;

// ================ VTK =====================
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkAppendFilter.h>
#include <vtkPointDataToCellData.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
// ==========================================
#include "omp.h"
#include <iostream>
#include <fstream>
#include <string>
#include "Cube.h"
#include "MultiProgressBar.h"

#define VERSION_MAJOR 1
#define VERSION_MINOR 0

#ifdef _WIN32
    #include "windows.h"
    #define BLACK			0
    #define BLUE			1
    #define GREEN			2
    #define CYAN			3
    #define RED				4
    #define MAGENTA			5
    #define BROWN			6
    #define LIGHTGRAY		7
    #define DARKGRAY		8
    #define LIGHTBLUE		9
    #define LIGHTGREEN		10
    #define LIGHTCYAN		11
    #define LIGHTRED		12
    #define LIGHTMAGENTA	13
    #define YELLOW			14
    #define WHITE			15
    static HANDLE   m_hConsole=GetStdHandle(STD_OUTPUT_HANDLE);
    static WORD     m_currentConsoleAttr;
    static CONSOLE_SCREEN_BUFFER_INFO csbi;
    #define COLOR_PURPLE ""
    #define COLOR_RED "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (RED & 0x0F) );cout<<""
    #define COLOR_GREEN "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (GREEN & 0x0F) );cout<<""
    #define COLOR_YELLOW "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (YELLOW & 0x0F) );cout<<""
    #define COLOR_BLUE "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (BLUE & 0x0F) );cout<<""
    #define COLOR_DEFAULT "";SetConsoleTextAttribute(m_hConsole, m_currentConsoleAttr );cout<<""
    #define ERROR_COUT "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (RED & 0x0F) );cout<<"Error: "<<COLOR_DEFAULT
    #define WARN_COUT "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (YELLOW & 0x0F) );cout<<"Warning: "<<COLOR_DEFAULT
#else
    // define color, this seems only work on MacOS and linux, doesn't work on windows
    #define ERROR_COUT "["<<"\033[31mError: "<<"\033[0m] "
    #define WARN_COUT "["<<"\033[33mWarning: "<<"\033[0m] "
    #define COLOR_PURPLE "\033[35m"
    #define COLOR_RED "\033[31m"
    #define COLOR_GREEN "\033[32m"
    #define COLOR_YELLOW "\033[33m"
    #define COLOR_BLUE "\033[34m"
    #define COLOR_DEFAULT "\033[0m"
#endif


namespace CMD
{                                                                                            
    bool bash_run(int argc, char** argv);
    void readSites_xyz(std::string fname_xyz, std::vector<FORWARD::STRUCT_SITE>& sites);
    /**
     * @brief Make sure cube bound as [xmin,xmax,ymin,ymax,zmin,zmax]
     * 
     * @param cube 
     */
    void makeSureCube(FORWARD::STRUCT_CUBE& cube);
    // void vtkUnstructuredGrid2Cubes(std::vector<FORWARD::STRUCT_CUBE>& cubes, vtkUnstructuredGrid* usg,double refT, double refRho, double alpha);
    void vtkUnstructuredGrid2Cubes(std::vector<FORWARD::STRUCT_CUBE>& cubes, vtkUnstructuredGrid* usg,std::string fieldName_rho,double rho0);
    vtkSmartPointer<vtkUnstructuredGrid> ReadUnstructuredGrid(std::string const& fileName, std::vector<std::string> fieldNames);
    // calculation mode: 0d, 1d, 2d, 3d
    #define CALCULATION_MODE_SINGLEPOINT 0
    #define CALCULATION_MODE_ONEDIMENSION 1
    #define CALCULATION_MODE_TWODIMENSION 2
    #define CALCULATION_MODE_THREEDIMENSION 3
    
    // variable selection
    #define VARIABLE_SELECTION_PTX 0
    #define VARIABLE_SELECTION_PHX 1
    #define VARIABLE_SELECTION_T 2
    #define VARIABLE_SELECTION_P 3
    #define VARIABLE_SELECTION_X 4
    #define VARIABLE_SELECTION_H 5
    #define VARIABLE_SELECTION_PT 6
    #define VARIABLE_SELECTION_PX 7
    #define VARIABLE_SELECTION_TX 8
    #define VARIABLE_SELECTION_PH 9
    #define VARIABLE_SELECTION_HX 10
    
    template<typename T>
    struct STRUCT_ARG
    {
        T value;
        bool have;
    };
    
    class cCMDarg
    {
    private:
        // bool m_haveZ0;
        STRUCT_ARG<int> m_threadNumOMP;
        STRUCT_ARG<int> m_InfoLevel;//level of output information
        STRUCT_ARG<double> m_rho0; //reference density(density - reference density = density contrast). SI unit
        STRUCT_ARG<string> m_FieldName_Density; //field name of density in vtu file
        STRUCT_ARG<string> m_vtu_T, m_xyz_sites, m_outFile;
        vector<string> m_valueR_str;
        STRUCT_ARG<bool> m_extractOnly; //only extract T and/or calculate rho, don't do gravity calculation. This is used to extract T and calculate gravity on server, because the ASPECT output with a lot of fields, other fields is not usefull for gravity calculation.
    public:
        cCMDarg(/* args */);
        ~cCMDarg();
    public:
        int m_CalculationMode;
        int m_VariableSelection;
    public:
        bool Parse(int argc, char** argv); //Parse arguments
        bool Validate(); // validate arguments and print corresponding error information
        /**
         * @brief Calculate gravity of thermal structure. Using \f$ \rho = \rho_0[1-\alpha(T - T_0)] \f$
         * 
         * @param vtuFile_T VTU file containing point data field "T"
         * @param xyzFile_sites xyz coordinate file of gravity calculation sites.
         * @param outputFile_grav Output file of gravity (gravity in unit mGal)
         * @param refT reference temperature (\f$ T_0 \f$, the unit is the same as "T" in vtu file)
         * @param refRho reference density (\f$ \rho_0 \f$, in unit of \f$ kg/m^3 \f$)
         * @param alpha thermal expansion (\f$ \alpha \f$, in unit of \f$ 1/K \f$)
         * @return true 
         * @return false 
         */
        // bool Temperature2Gravity(string vtuFile_T, string xyzFile_sites, string outputFile_grav, 
        //     double refT, double refRho, double alpha);
        bool getCubes_VTU(string vtuFile, string fieldName_rho,std::vector<FORWARD::STRUCT_CUBE>& cubes,double rho0);
        bool Gravity_Cubes(std::vector<FORWARD::STRUCT_CUBE> cubes,string xyzFile_sites, string outputFile_grav);
        bool calGravity(string vtuFile, string fieldName_rho, string xyzFile_sites, string outputFile_grav, double rho0);
    private:
        bool GetOptionValue(int opt, char* optarg, STRUCT_ARG<double>& arg);
        template<typename T>
        vector<T> linspace(T xmin, T xmax, T dx);
        inline void Error(string errInfo){
            cout<<ERROR_COUT<<errInfo<<COLOR_DEFAULT<<endl;
            exit(0);
        };
        inline void Warn(string warnInfo, int level=1){
            if(m_InfoLevel.value<=0)return; //Silent Mode
            if(level<m_InfoLevel.value)
            cout<<WARN_COUT<<warnInfo<<COLOR_DEFAULT<<endl;
        };
        inline void Info(string info, string tag="Info", int level=1){
            if(m_InfoLevel.value<=0)return; //Silent Mode
            if(level<m_InfoLevel.value)
            cout<<"["<<tag<<": ]"<<info<<COLOR_DEFAULT<<endl;
        };
    };
    bool isNum(string str);
    vector<string> string_split(string s, string delimiter);
    
    static void StartText()
    {
        //30: black  31:red  32:green  33:yellow  34:blue  35:purple  36:darkgreen
        cout<<COLOR_YELLOW;       //print text in yellow color
        cout << "***************************************************\n";
        cout << "*                 program vtk2grav                *\n";
        cout << "*                 ~~~~~~~ ~~~~~~~                 *\n";
        cout << "*  Version: "<<vtk2grav_VERSION<<"                     *\n";
        cout << "*                                                 *\n";
        cout << "*  Calculate gravity of 3D thermal structure      *\n";
        cout << "*  - Upward Z-axis indicates positive direction   *\n";
        cout << "*  - VTK cell type: VTK_VOXELand  VTK_HEXAHEDRO   *\n";
        cout << "*  unit:                                          *\n";
        cout << "*      temperature - deg.C,                       *\n";
        cout << "*      density - kg/m3,                           *\n";
        cout << "*      gravity - mGal                             *\n";
        cout << "*                                                 *\n";
        cout << "* (c) Zhikui Guo, GEOMAR, "<<vtk2grav_DATE<<", Kiel        *\n";
        cout << "*                                                 *\n";
        cout << "***************************************************\n";
        cout << "\n";
        cout<<COLOR_DEFAULT;
                                                                                                                                                                                                                        
    }
    static void helpINFO()
    {
        string version=vtk2grav_VERSION;
        string author="Zhikui Guo";
        string locus="GEOMAR, Germany";
        string email="zguo@geomar.de";
        unsigned int wordWidth=20;
        // time_t now=time(0);
        // char* now_str=ctime(&now);
        string now_str=vtk2grav_DATE;

        //30: black  31:red  32:green  33:yellow  34:blue  35:purple  36:darkgreen
        cout<<"========================== vtk2grav ==========================="<<std::endl;
        cout<<"vtk2grav, calculate gravity of 3D thermal structure."<<std::endl;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Author "<<COLOR_GREEN<<author<<COLOR_DEFAULT<<std::endl;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Locus "<<COLOR_GREEN<<locus<<COLOR_DEFAULT<<std::endl;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Date "<<COLOR_GREEN<<now_str<<COLOR_DEFAULT<<std::endl;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Version "<<COLOR_GREEN<<version<<COLOR_DEFAULT<<std::endl;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Email "<<COLOR_GREEN<<email<<COLOR_DEFAULT<<std::endl;
        cout<<"============================================================"<<std::endl;
        cout<<COLOR_BLUE<<"Usage: vtk2grav [options]"<<COLOR_DEFAULT<<std::endl;
        cout<<"options:"<<std::endl;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -h "<<COLOR_BLUE<<"List descriptions of usage and available arguments"<<COLOR_DEFAULT<<std::endl;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -v "<<COLOR_BLUE<<"Print swEOS version number"<<COLOR_DEFAULT<<std::endl;
        cout<<COLOR_RED<<setw(wordWidth)<<setiosflags(ios::left)<<"  -i "<<COLOR_BLUE<<"Input vtk file of temperature field."<<COLOR_DEFAULT<<std::endl;
        // cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -E "<<COLOR_BLUE<<"Only extract T and save to a new vtu file"<<COLOR_DEFAULT<<std::endl;
        cout<<COLOR_RED<<setw(wordWidth)<<setiosflags(ios::left)<<"  -p "<<COLOR_BLUE<<"Position file of calculation points (xyz)."<<COLOR_DEFAULT<<std::endl;
        // cout<<COLOR_GREEN<<setw(wordWidth)<<setiosflags(ios::left)<<"  -T "<<COLOR_BLUE<<"Set reference temperature(deg.C), e.g. -T 1300."<<COLOR_DEFAULT<<std::endl;
        cout<<COLOR_GREEN<<setw(wordWidth)<<setiosflags(ios::left)<<"  -D "<<COLOR_BLUE<<"Set reference density(kg/m3) which will be subtracted from the input density field, e.g. -D 3300. Default is 0"<<COLOR_DEFAULT<<std::endl;
        cout<<COLOR_GREEN<<setw(wordWidth)<<setiosflags(ios::left)<<"  -F "<<COLOR_BLUE<<"Specify field name of density in the vtu file, default is density."<<COLOR_DEFAULT<<std::endl;
        // cout<<COLOR_GREEN<<setw(wordWidth)<<setiosflags(ios::left)<<"  -A "<<COLOR_BLUE<<"Set thermal expansion coeff., e.g. -A 1E5."<<COLOR_DEFAULT<<std::endl;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -o "<<COLOR_BLUE<<"Set output file, e.g. -o density.txt "<<COLOR_DEFAULT<<std::endl;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -t "<<COLOR_BLUE<<"Set number of thread for parallel computing."<<COLOR_DEFAULT<<std::endl;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -V "<<COLOR_BLUE<<"Set output information level, e.g. -V 0 means silent mode."<<COLOR_DEFAULT<<std::endl;
        cout<<COLOR_GREEN<<"Example: "<<COLOR_BLUE<<"vtk2grav -i T_density.vtu -F density -D 3300  -p sites.xyz -o grav.txt"<<COLOR_DEFAULT<<std::endl;
    }
    static void StartText_artASCII()
    {
        // ANSI shadow
        cout<<COLOR_GREEN<<"██╗   ██╗████████╗██╗  ██╗██████╗  ██████╗ ██████╗  █████╗ ██╗   ██╗\n"
        <<"██║   ██║╚══██╔══╝██║ ██╔╝╚════██╗██╔════╝ ██╔══██╗██╔══██╗██║   ██║\n"
        <<"██║   ██║   ██║   █████╔╝  █████╔╝██║  ███╗██████╔╝███████║██║   ██║\n"
        <<"╚██╗ ██╔╝   ██║   ██╔═██╗ ██╔═══╝ ██║   ██║██╔══██╗██╔══██║╚██╗ ██╔╝\n"
        <<" ╚████╔╝    ██║   ██║  ██╗███████╗╚██████╔╝██║  ██║██║  ██║ ╚████╔╝ \n"
        <<"  ╚═══╝     ╚═╝   ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝  ╚═══╝  \n"
        <<COLOR_DEFAULT<<std::endl;                                                                                                                                                                             
    }
}
#endif
