# Set an environment variable CMAKE_PREFIX_PATH with paths for find_library etc to search
# or add: HINTS "/path/to/lib/dir" at the end of the find_library call

add_executable(fuego.x)
target_sources(fuego.x PRIVATE main.cc
                               driver_premixed.cc 
                               driver_diffusion.cc 
                               driver_diffusion_table.cc 
                               fuego.cc 
                               streams.cc)

if(APPLE)
    find_library(ACCELERATE Accelerate)
    set_property(TARGET fuego.x PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)  # adds LC_RPATH to executable
endif()

find_library(CANTERA NAMES cantera_shared cantera)
find_library(RADLIB  NAMES radlib_shared  radlib)
find_library(SOOTLIB NAMES sootlib)
find_library(KINSOL sundials_kinsol)
find_library(CVODES sundials_cvodes)
find_library(YAML   yaml-cpp)

target_include_directories(fuego.x PRIVATE "$ENV{CONDA_PREFIX}/include")
target_include_directories(fuego.x PRIVATE "/usr/local/include")
target_include_directories(fuego.x PRIVATE "/opt/radlib/include")
target_include_directories(fuego.x PRIVATE "/opt/sootlib/include")

target_link_directories(fuego.x PRIVATE "$ENV{CONDA_PREFIX}/lib")
target_link_directories(fuego.x PRIVATE "/usr/local/lib")
target_link_directories(fuego.x PRIVATE "/opt/radlib/lib")
target_link_directories(fuego.x PRIVATE "/opt/sootlib/lib")

target_link_libraries(fuego.x ${CANTERA} ${YAML} ${RADLIB} ${SOOTLIB} ${KINSOL} ${CVODES} ${ACCELERATE})

install(TARGETS fuego.x DESTINATION "${PROJECT_SOURCE_DIR}/run")
