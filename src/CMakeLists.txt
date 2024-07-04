# Set an environment variable CMAKE_PREFIX_PATH with paths for find_library etc to search
# or add: HINTS "/path/to/lib/dir" at the end of the find_library call

add_executable(ignis.x)
target_sources(ignis.x PRIVATE main.cc
                               driver_premixed.cc 
                               driver_diffusion.cc 
                               driver_diffusion_table.cc 
                               driver_flamelet.cc 
                               ignis.cc 
                               streams.cc)

set_property(TARGET ignis.x PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)  # adds LC_RPATH to executable

find_package(radlib)
find_package(sootlib)
find_package(HighFive REQUIRED)
find_library(CANTERA NAMES cantera_shared cantera)
find_library(KINSOL sundials_kinsol)
find_library(CVODES sundials_cvodes)
find_library(YAML   yaml-cpp)

target_include_directories(ignis.x PRIVATE "$ENV{CONDA_PREFIX}/include")
target_include_directories(ignis.x PRIVATE "/usr/local/include")

target_link_directories(ignis.x PRIVATE "$ENV{CONDA_PREFIX}/lib")
target_link_directories(ignis.x PRIVATE "/usr/local/lib")

target_link_libraries(ignis.x radlib::radlib sootlib::sootlib ${CANTERA} ${YAML} ${KINSOL} ${CVODES} )
target_link_libraries(ignis.x HighFive)

install(TARGETS ignis.x DESTINATION "${PROJECT_SOURCE_DIR}/run")
