# 0. CMake Minimum version
cmake_minimum_required(VERSION 3.3...3.12 FATAL_ERROR)


set(PROGRAM_NAME vtk2grav)

project(vtk2grav LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(vtk2grav_MAJOR_VERSION 1)
set(vtk2grav_MINOR_VERSION 0)
set(vtk2grav_PATCH_VERSION 4)
set(vtk2grav_EXTRA_VERSION "")
set(vtk2grav_EXTRA_VERSION_TEXI "${vtk2grav_EXTRA_VERSION}")
if(NOT vtk2grav_RELEASE)
  find_package(Git)
  if(GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ERROR_QUIET
                    OUTPUT_VARIABLE GIT_COMMIT_HASH
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()
  if(GIT_COMMIT_HASH)
    set(vtk2grav_EXTRA_VERSION "${vtk2grav_EXTRA_VERSION}-git-${GIT_COMMIT_HASH}")
  endif()
  set(vtk2grav_EXTRA_VERSION_TEXI "${vtk2grav_EXTRA_VERSION_TEXI} (development version)")
endif()
set(vtk2grav_VERSION "${vtk2grav_MAJOR_VERSION}.${vtk2grav_MINOR_VERSION}")
set(vtk2grav_VERSION "${vtk2grav_VERSION}.${vtk2grav_PATCH_VERSION}${vtk2grav_EXTRA_VERSION}")
set(vtk2grav_SHORT_LICENSE "GNU General Public License")

if(NOT DATE)
  string(TIMESTAMP DATE "%Y-%m-%d")
endif()
set(vtk2grav_DATE "${DATE}")

if(NOT vtk2grav_HOST)
  execute_process(COMMAND hostname OUTPUT_VARIABLE HOSTNAME
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT HOSTNAME)
    set(HOSTNAME "unknown")
  endif()
  set(vtk2grav_HOST "${HOSTNAME}")
endif()

if(NOT vtk2grav_PACKAGER)
  execute_process(COMMAND whoami OUTPUT_VARIABLE PACKAGER
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT PACKAGER)
    set(PACKAGER "unknown")
  endif()
  string(REPLACE "\\" " " PACKAGER ${PACKAGER})
  set(vtk2grav_PACKAGER "${PACKAGER}")
endif()

if(APPLE)
  set(vtk2grav_OS "MacOSX")
elseif(CYGWIN OR MSYS)
  # detect if we use the MinGW compilers on Cygwin - if we do, handle the build
  # as a pure Windows build
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
     CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpmachine
                    OUTPUT_VARIABLE CXX_COMPILER_MACHINE
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(CXX_COMPILER_MACHINE MATCHES "mingw")
      set(vtk2grav_OS "Windows")
      set(WIN32 1)
      add_definitions(-DWIN32)
    endif()
  endif()
else()
  set(vtk2grav_OS "${CMAKE_SYSTEM_NAME}")
endif()

# Using OpenMP
include_directories(/usr/local/Cellar/libomp/11.0.0/include /usr/local/Cellar/libomp/11.1.0/include)
find_package(OpenMP REQUIRED)
if(APPLE)
    set(CMAKE_C_COMPILER clang)
    set(CMAKE_CXX_COMPILER clang++)

    if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY omp)
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY omp)
    endif()
    link_directories("/usr/local/lib")
endif()
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
#windows vs2017 doesn't need opm library link 
set(LIBomp "")
if (APPLE)
  set(LIBomp "omp")
endif()

# find vtk
find_package(VTK COMPONENTS 
  vtkCommonColor
  vtkCommonCore
  vtkCommonDataModel
  vtkFiltersCore
  vtkFiltersSources
  vtkIOLegacy
  vtkIOXML
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2 QUIET
)
if (NOT VTK_FOUND)
  message("Skipping ReadAllUnstructuredGridTypes: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
message(STATUS "VTK_LIBRARIES: ${VTK_LIBRARIES}")

#2. Source code dir
include_directories(include)
include_directories(Forward/include)
aux_source_directory(src SRC_MAIN)
add_subdirectory(src)
# Add forward library
add_subdirectory(Forward/src)

# add Executable target
add_executable(${PROGRAM_NAME} ${SRC_MAIN} main.cpp)
# use lib of forward functions
target_link_libraries(${PROGRAM_NAME} forward ${LIBomp} ${VTK_LIBRARIES})
# target_link_libraries(ReadAllUnstructuredGridTypes PRIVATE ${VTK_LIBRARIES})
# vtk_module_autoinit is needed
# vtk_module_autoinit(
#   TARGETS ${PROGRAM_NAME}
#   MODULES ${VTK_LIBRARIES}
#   )

# Add test 
set(EnableTest ON)
if(EnableTest)
  ENABLE_TESTING()
  ADD_SUBDIRECTORY(tests)
  # test 1
  SET(testNames "test_cube")
  foreach(testName ${testNames})
    add_executable(${testName} tests/${testName}.cpp)
    target_link_libraries(${testName} forward ${LIBomp})
    add_test(NAME ${testName} COMMAND ${testName})
  endforeach()
endif(EnableTest)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
               ${CMAKE_CURRENT_SOURCE_DIR}/include/version.h)

########################################################
###                   INSTALL                        ###
########################################################
if(WIN32)
if(ENABLE_OS_SPECIFIC_INSTALL)
  set(vtk2grav_BIN .)
  set(vtk2grav_LIB .)
  set(vtk2grav_DOC .)
  set(vtk2grav_MAN .)
  set(vtk2grav_INC .)
else()
  include(GNUInstallDirs)
  set(vtk2grav_BIN ${CMAKE_INSTALL_BINDIR})
  set(vtk2grav_LIB ${CMAKE_INSTALL_LIBDIR})
  set(vtk2grav_DOC ${CMAKE_INSTALL_DOCDIR})
  set(vtk2grav_MAN ${CMAKE_INSTALL_MANDIR}/man1)
  set(vtk2grav_INC ${CMAKE_INSTALL_INCLUDEDIR})
endif()
if(CYGWIN)
  unix2dos(vtk2grav_API)
  if(ENABLE_PRIVATE_API)
    unix2dos(vtk2grav_PRIVATE_API)
  endif()
    unix2dos(WELCOME_FILE)
    unix2dos(SDK_FILE)
    unix2dos(LICENSE_FILE)
    unix2dos(CREDITS_FILE)
    unix2dos(CHANGELOG_FILE)
    unix2dos(TUTORIAL_FILES)
    foreach(DIR ${DEMOS_DIRS})
      file(GLOB DEMO_FILES ${DIR}/?*.*)
      unix2dos(DEMO_FILES)
    endforeach()
  endif()
elseif(APPLE AND ENABLE_OS_SPECIFIC_INSTALL)
  # set these so that the files get installed nicely in the MacOSX
  # .app bundle
  message("macOS install path config")
  set(vtk2grav_BIN ../MacOS)
  set(vtk2grav_LIB ../MacOS)
  set(vtk2grav_DOC ../../..)  
  set(vtk2grav_MAN ../../..)
  set(vtk2grav_INC ../MacOS)
elseif(APPLE AND NOT ENABLE_OS_SPECIFIC_INSTALL)
  # set these so that the files get installed nicely in the MacOSX
  # .app bundle
  set(vtk2grav_BIN .)
  set(vtk2grav_LIB ./lib)
  set(vtk2grav_DOC ./doc)
  set(vtk2grav_MAN ./man)
  set(vtk2grav_INC ./include)
  else()
  include(GNUInstallDirs)
  set(vtk2grav_BIN ${CMAKE_INSTALL_BINDIR})
  set(vtk2grav_LIB ${CMAKE_INSTALL_LIBDIR})
  set(vtk2grav_DOC ${CMAKE_INSTALL_DOCDIR})
  set(vtk2grav_MAN ${CMAKE_INSTALL_MANDIR}/man1)
  set(vtk2grav_INC ${CMAKE_INSTALL_INCLUDEDIR})
endif() 

install (TARGETS ${PROGRAM_NAME} DESTINATION ${vtk2grav_BIN})
if(WIN32)
  # install vtk libs
  FILE(GLOB VTK_LIBS "${VTK_DIR}/../../../bin/*.dll")
  install (FILES ${VTK_LIBS} DESTINATION ${vtk2grav_BIN})
  install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Release/ DESTINATION ${vtk2grav_BIN})
  install (DIRECTORY ${PROJECT_SOURCE_DIR}/../commandline/build/Release/ DESTINATION "commandline")
else()
  FILE(GLOB VTK_LIBS "${VTK_DIR}/../../*vtk*.*")
  install (FILES ${VTK_LIBS} DESTINATION ${vtk2grav_LIB})
endif()

########################################################
###                   CPACK                          ###
########################################################
set(CPACK_PACKAGE_VENDOR "Zhikui Guo (zguo@geomar.de)")
set(CPACK_PACKAGE_VERSION_MAJOR ${vtk2grav_MAJOR_VERSION})
set(CPACK_PACKAGE_VERSION_MINOR ${vtk2grav_MINOR_VERSION})
set(CPACK_PACKAGE_VERSION_PATCH ${vtk2grav_PATCH_VERSION})
set(CPACK_PACKAGE_DESCRIPTION_FILE ${WELCOME_FILE})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "swEOS, a desktop application and API for salt-water (H2ONaCl) EOS and thermodynamic properties.")
if(vtk2grav_EXTRA_VERSION MATCHES "-git.*") # so that we'll overwrite the archives
  set(CPACK_PACKAGE_FILE_NAME ${PROGRAM_NAME}-${vtk2grav_OS}-Installer)
  set(CPACK_SOURCE_PACKAGE_FILE_NAME ${PROGRAM_NAME}-git-source)
else()
  set(CPACK_PACKAGE_FILE_NAME ${PROGRAM_NAME}-${vtk2grav_VERSION}-${vtk2grav_OS})
  set(CPACK_SOURCE_PACKAGE_FILE_NAME ${PROGRAM_NAME}-${vtk2grav_VERSION}-source)
endif()
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROGRAM_NAME}")
set(CPACK_RESOURCE_FILE_LICENSE ${LICENSE_FILE})
set(CPACK_RESOURCE_FILE_README ${WELCOME_FILE})
set(CPACK_RESOURCE_FILE_WELCOME ${WELCOME_FILE})
set(CPACK_PACKAGE_CONTACT "zguo@geomar.de")
set(CPACK_PACKAGE_EXECUTABLE "${PROGRAM_NAME}")
if(ENABLE_PACKAGE_STRIP)
  set(CPACK_STRIP_FILES TRUE)
else()
  set(CPACK_STRIP_FILES FALSE)
endif()
set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_CURRENT_BINARY_DIR}" )
if (WIN32)
    if (USE_WIX_TOOLSET)
        set(CPACK_GENERATOR "WIX") # this need WiX Tooset installed and a path to candle.exe
    elseif(USE_7Z_CPACK)
    set(CPACK_GENERATOR "7Z")
    else ()
        set(CPACK_GENERATOR "NSIS") # this needs NSIS installed, and available
    endif ()
    set(CPACK_SOURCE_GENERATOR "ZIP")
    # set icon
    set(CPACK_NSIS_DISPLAY_NAME  "${PROGRAM_NAME}")
    set(CPACK_NSIS_CONTACT "zguo@geomar.de")
    set(CPACK_PACKAGE_EXECUTABLES "${PROGRAM_NAME}" "swEOS")
    set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${PROGRAM_NAME}-${vtk2grav_VERSION}")
    set(CPACK_NSIS_MENU_LINKS "https://www.scibyte.cn" "Homepage for ${PROGRAM_NAME}")
    set(CPACK_NSIS_URL_INFO_ABOUT "https://www.scibyte.cn")
    set(CPACK_NSIS_HELP_LINK "https://www.scibyte.cn")
    set(CPACK_NSIS_MUI_ICON "${PROJECT_SOURCE_DIR}/icon.ico")
    # set(CPACK_NSIS_MUI_HEADERIMAGE "${PROJECT_SOURCE_DIR}/../manual/home/static/images/vtk2grav_mac.jpg")
    set(CPACK_NSIS_WELCOME_TITLE "swEOS: salt-water (H2O-NaCl) EOS")
    set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
elseif ( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
    set(CPACK_GENERATOR "DragNDrop")
    set(CPACK_BUNDLE_APPLE_CERT_APP "Developer ID Application: Zhikui GUO (AFCJBNSQCL)")
else ()
    set(CPACK_GENERATOR "ZIP")
    set(CPACK_SOURCE_GENERATOR "TGZ") 
endif ()

include(CPack)

message(STATUS "")
message(STATUS "${PROGRAM_NAME} ${${PROGRAM_NAME}_VERSION} has been configured for ${vtk2grav_OS}")
message(STATUS "")
message(STATUS " * Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS " * C compiler: " ${CMAKE_C_COMPILER})
message(STATUS " * C++ compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS " * Install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS "")
