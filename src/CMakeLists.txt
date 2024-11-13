# Set an environment variable CANTERA_INCLUDE_PATH, e.g., /usr/local/cantera/include
# Set an environment variable CMAKE_PREFIX_PATH with paths for find_library etc to search
# or add: HINTS "/path/to/lib/dir" at the end of the find_library call
# NOTE: below, target include directories and link directories has usr/local before CONDA_PREFIX.
#       You should have the same order in your CMAKE_PREFIX_PATH

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
find_library(CANTERA NAMES cantera_shared cantera)
find_library(YAML   yaml-cpp)
find_package(HighFive REQUIRED)
find_library(SUNKIN  sundials_kinsol)
find_library(SUNCV   sundials_cvodes)
find_library(SUNCORE sundials_core)
find_library(SUNNV   sundials_nvecserial)

target_include_directories(ignis.x PRIVATE "$ENV{CANTERA_INCLUDE_PATH}")
target_include_directories(ignis.x PRIVATE "/usr/local/include")
target_include_directories(ignis.x PRIVATE "$ENV{CONDA_PREFIX}/include")

target_link_directories(ignis.x PRIVATE "/usr/local/lib")
target_link_directories(ignis.x PRIVATE "$ENV{CONDA_PREFIX}/lib")

target_link_libraries(ignis.x radlib::radlib sootlib::sootlib)
target_link_libraries(ignis.x ${CANTERA} ${YAML} HighFive)
target_link_libraries(ignis.x ${SUNKIN} ${SUNCV} ${SUNCORE} ${SUNNV} )

install(TARGETS ignis.x DESTINATION "${PROJECT_SOURCE_DIR}/run")
