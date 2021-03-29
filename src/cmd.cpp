#include "cmd.h"
// using namespace CMD;
namespace CMD
{
  bool bash_run(int argc, char** argv)
  {
    #ifdef _WIN32
      // set terminal as black(bg)+white(fg) model
      system("color 07"); //see https://www.geeksforgeeks.org/how-to-print-colored-text-in-c/
      GetConsoleScreenBufferInfo(m_hConsole, &csbi); 
      m_currentConsoleAttr = csbi.wAttributes;
      int width = (int)(csbi.srWindow.Right-csbi.srWindow.Left+1);
      // int height = (int)(csbi.srWindow.Bottom-csbi.srWindow.Top+1);
      if(width>119)
      {
          StartText_artASCII();
      }else
      {
          StartText();
      }
    #else
      struct winsize w;
      ioctl(0, TIOCGWINSZ, &w);
      if(w.ws_col>119)
      {
          StartText_artASCII();
      }else
      {
          StartText();
      }
    #endif
    
    // helpINFO();
    //parse arguments and check 
    cCMDarg arg;
    if(!arg.Parse(argc, argv)) return false;
    if(!arg.Validate()) return false;
    return true;
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
  void vtkUnstructuredGrid2Cubes(std::vector<FORWARD::STRUCT_CUBE>& cubes, vtkUnstructuredGrid* usg, double refT, double refRho, double alpha)
  {
      cubes.clear();
      FORWARD::STRUCT_CUBE cube;
      vtkSmartPointer<vtkPointDataToCellData> p2c = vtkSmartPointer<vtkPointDataToCellData>::New();
      p2c->SetInputData(usg);
      p2c->PassPointDataOn();
      p2c->Update();
      

      vtkSmartPointer<vtkDataArray> Array_T = p2c->GetOutput()->GetCellData()->GetArray("T");
      if(!Array_T)
      {
          cout<<"ERROR: T doesn't exist in the vtu file"<<endl;
          exit(0);
      }
      int nCells = p2c->GetOutput()->GetNumberOfCells();
      vtkSmartPointer<vtkCell> cell;
      vtkSmartPointer<vtkPoints> points;
      int cellType;
      int numUnsupportedCell = 0;
      double T;
      
      vtkSmartPointer<vtkDoubleArray> Array_density = vtkSmartPointer<vtkDoubleArray>::New();
      Array_density->SetNumberOfValues(nCells);
      Array_density->SetName("density");
      for (int i = 0; i < nCells; i++)
      {
          cell = usg->GetCell(i);
          cellType = cell->GetCellType();
          cube.density = 0;
          if(cellType==VTK_VOXEL)
          {
              cell->GetPoints();
              // calculate density from T
              // cube.density = density->GetTuple1(i);
              T = Array_T->GetTuple1(i);
              cube.density = refRho*alpha*(T - refT);
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
              // cube.density = density->GetTuple1(i);
              T = Array_T->GetTuple1(i);
              cube.density = refRho*alpha*(T - refT);
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
              numUnsupportedCell ++;
              if(numUnsupportedCell>100)
              {
                std::cout<<"More than 100 unsupported cells exist. Please check the vtu file of temperature, only support VTK_VOXEL and VTK_HEXAHEDRON (ASPECT)"<<endl;
                exit(0);
              }
          }
          Array_density->SetValue( i, cube.density); //used to save density to a vtu file
      }
      p2c->GetOutput()->GetCellData()->AddArray(Array_density);
      // write to file 
      // Write file
      vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
      vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
      writer->SetFileName("T_density.vtu");
      writer->SetInputData(p2c->GetOutput());
      writer->Write();
  }
  vtkSmartPointer<vtkUnstructuredGrid> ReadUnstructuredGrid(std::string const& fileName, std::vector<std::string> fieldNames)
  {
      vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
      std::string extension = "";
      if (fileName.find_last_of(".") != std::string::npos)
      {
      extension = fileName.substr(fileName.find_last_of("."));
      }
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

  cCMDarg::cCMDarg(/* args */)
  :m_threadNumOMP{omp_get_max_threads(), false},
  m_InfoLevel{5, false},
  m_refT{1300,false}, 
  m_refRho{3300,false},
  m_alpha{1E-5,false},
  m_outFile{"density.txt",false}
  {
    
  }

  cCMDarg::~cCMDarg()
  {
  }
  bool isNum(string str)
  {
      stringstream sin(str);
      double d;
      char c;
      if(!(sin >> d))
          return false;
      if (sin >> c)
          return false;
      return true;
  }
  bool cCMDarg::GetOptionValue(int opt, char* optarg, STRUCT_ARG<double>& arg)
  {
    string optarg_str=optarg;
    if(isNum(optarg_str))
    {
      arg.value=atof(optarg);
    }else
    { 
      char optCh=opt;
      cout<<ERROR_COUT<<"Option of -"<<optCh<<" argument is empty or cannot be recognized"<<endl;
      return false;
    }
    arg.have = true;
    return true;
  }
  
  bool cCMDarg::Parse(int argc, char** argv)
  {
    if(argc<2)return false; //there is no arguments
    int opt; 
    const char *optstring = "i:o:p:T:D:A:t:V:vh"; // set argument templete
    int option_index = 0;
    int valid_args=0;
    STRUCT_ARG<double> doubleOptValue;
    while ((opt = getopt(argc, argv, optstring)) != -1) 
    {
      if(opt!='?')
      {
        valid_args++;
      }else
      {
        cout<<WARN_COUT<<"Unknown option "<<argv[optind-1]<<endl;
      }
      switch (opt)
      {
      case 'h':
        helpINFO();
        exit(0);
        break;
      case 'v':
        cout<<"Version: "<<VERSION_MAJOR<<"."<<VERSION_MINOR<<endl;
        exit(0);
        break;
      case 'i':  //input vtu file m_valueO=optarg;
        m_vtu_T.value = optarg;
        m_vtu_T.have = true;
        break;
      case 'p':  //calculation point position file, xyz
        m_xyz_sites.value = optarg;
        m_xyz_sites.have = true;
        break;
      case 'o':  //density output file, single column
        m_outFile.value = optarg;
        m_outFile.have = true;
        break;
      case 'T':  //reference temperature (deg.C)
        if(!GetOptionValue(opt, optarg, m_refT))return false;
        break;
      case 'D':  //reference density (kg/m3)
        if(!GetOptionValue(opt, optarg, m_refRho))return false;
        break;
      case 'A':  //Coefficient of thermal expansion, alpha
        if(!GetOptionValue(opt, optarg, m_alpha))return false;
        break;
      case 't':  //thread number
        if(!GetOptionValue(opt, optarg, doubleOptValue))return false;
        m_threadNumOMP.value = (int)(doubleOptValue.value);
        if(m_threadNumOMP.value>omp_get_max_threads())m_threadNumOMP.value=omp_get_max_threads();
        if(m_threadNumOMP.value<1)m_threadNumOMP.value=1;
        m_threadNumOMP.have = true;
        break;
      case 'V':  //output information level
        if(!GetOptionValue(opt, optarg, doubleOptValue))return false;
        m_InfoLevel.value = (int)(doubleOptValue.value);
        m_InfoLevel.have = true;
        break;
      default:
        break;
      }
    }
    return true;
  }
  
  bool cCMDarg::Validate()
  {
    //check required arguments
    if(!m_vtu_T.have)
    {
      Error("must specify temperature vtu file by -i argument");
    }
    if(!m_xyz_sites.have)
    {
      Error("must specify calculation points position file (xyz) by -p argument");
    }
    // if(!m_refT.have)
    // {
    //   Warn("The reference temperature is not specified by -T argument, use default value: "+to_string(m_refT.value));
    // }
    // if(!m_refRho.have)
    // {
    //   Warn("The reference density is not specified by -D argument, use default value: "+to_string(m_refRho.value));
    // }
    // if(!m_alpha.have)
    // {
    //   Warn("The thermal expansion (alpha) is not specified by -A argument, use default value: "+to_string(m_alpha.value));
    // }
    // if(!m_outFile.have)
    // {
    //   Warn("The output file name is not specified by -o argument, use default filename: "+m_outFile.value);
    // }
    string modelInfo = "Thread number: "+to_string(m_threadNumOMP.value);
    modelInfo +=", Ref. T = "+to_string(m_refT.value)+" deg.C";
    modelInfo +=", Ref. rho = "+to_string(m_refRho.value)+" kg/m3";
    modelInfo +=", alpha = "+to_string(m_alpha.value);
    Info(modelInfo);
    // run 
    Temperature2Gravity(m_vtu_T.value, m_xyz_sites.value, m_outFile.value, m_refT.value, m_refRho.value,m_alpha.value);
    return true;
  }
  bool cCMDarg::Temperature2Gravity(string vtuFile_T, string xyzFile_sites, string outputFile_grav, 
            double refT, double refRho, double alpha)
  {
    // 1. read vtu file
    std::vector<string> filedNames{"T"};
    Info("Reading vtu file ...");
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = ReadUnstructuredGrid(vtuFile_T, filedNames);
    // 2. convert unstructuredGrid FORWARD::STRUCT_CUBE 
    std::vector<FORWARD::STRUCT_CUBE> cubes;
    vtkUnstructuredGrid2Cubes(cubes,unstructuredGrid, refT, refRho, alpha);
    // 3. read calculation point sites
    Info("Reading xyz file ...");
    std::vector<FORWARD::STRUCT_SITE> sites;
    readSites_xyz(xyzFile_sites, sites);
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
    ofstream fout(outputFile_grav);
    for (size_t i = 0; i < grav.size(); i++)
    {
        fout<<grav[i]<<endl;
    }
    fout.close();
    return true;
  }
  template<typename T>
  vector<T> cCMDarg::linspace(T xmin, T xmax, T dx)
  {
    vector<T> tmp;
    for (T x = xmin; x <= xmax; x=x+dx)
    {
      tmp.push_back(x);
    }
    return tmp;
  }
  vector<string> string_split (string s, string delimiter) {
      size_t pos_start = 0, pos_end, delim_len = delimiter.length();
      string token;
      vector<string> res;
      while ((pos_end = s.find (delimiter, pos_start)) != string::npos) 
      {
          token = s.substr (pos_start, pos_end - pos_start);
          pos_start = pos_end + delim_len;
          res.push_back (token);
      }
      res.push_back (s.substr (pos_start));
      return res;
  }
}