#cmake_minimum_required(VERSION 2.8)
#add_executable(SBS SBS.cc)
#target_link_libraries(SBS ${Geant4_LIBRARIES})
#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(APEXG4MC)
find_package(Geant4 REQUIRED)
##
## Get the ROOT config
include(FindROOT.cmake)
set(INC_ROOT ${ROOT_INCLUDE_DIR})
set(LIB_ROOT ${ROOT_LIBRARY_DIR})

#include_directories(${Geant4_INCLUDE_DIRS} ${INC_ROOT})

add_definitions(${Geant4_CXX_FLAGS})
set(CMAKE_CXX_FLAGS ${Geant4_CXX_FLAGS})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(GEANT4_INSTALL_DATA "Data" ON)
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include ${INC_ROOT})
#link_directories($

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cc ${sources} ${headers})
target_link_libraries(${PROJECT_NAME} ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
#
if(EXISTS "Septa-JB_map.table")
  file(SHA256 Septa-JB_map.table hash)
  if(hash EQUAL 60fba61d77aafc44d7ed4e38b2fb6615d66bf6cc68cdbe01640f95fbe2644036)
    message(STATUS "Septum map already downloaded")
  else()
    message("Downloading Septum field map...")
    file(DOWNLOAD http://nuclear.gla.ac.uk/~david/Septa-JB_map.table ${CMAKE_CURRENT_BINARY_DIR}/Septa-JB_map.table SHOW_PROGRESS)
  endif()
else()
  message("Downloading Septum field map...")
  file(DOWNLOAD http://nuclear.gla.ac.uk/~david/Septa-JB_map.table ${CMAKE_CURRENT_BINARY_DIR}/Septa-JB_map.table SHOW_PROGRESS)
endif()
#----------------------------------------------------------------------------
ADD_CUSTOM_TARGET(link_target ALL
                  COMMAND ${CMAKE_COMMAND} -E create_symlink ../macros ./macros)
