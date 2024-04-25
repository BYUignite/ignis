# Set an environment variable CMAKE_PREFIX_PATH with paths for find_library etc to search
# or add: HINTS "/path/to/lib/dir" at the end of the find_library call

add_executable(flame.x)
target_sources(flame.x PRIVATE main.cc
                               driver_premixed.cc 
                               driver_diffusion.cc 
                               driver_diffusion_table.cc 
                               flame.cc 
                               streams.cc)

if(APPLE)
    find_library(ACCELERATE Accelerate)
    set_property(TARGET flame.x PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)  # adds LC_RPATH to executable
endif()

#<<<<<<< HEAD
find_library(CANTERA NAMES cantera_shared cantera)
find_library(RADLIB  NAMES radlib_shared  radlib)
find_library(SOOTLIB NAMES sootlib)
find_library(KINSOL sundials_kinsol)
find_library(CVODES sundials_cvodes)

target_include_directories(flame.x PRIVATE "$ENV{CONDA_PREFIX}/include")
target_include_directories(flame.x PRIVATE "/usr/local/include")
target_include_directories(flame.x PRIVATE "/opt/radlib/include")
target_include_directories(flame.x PRIVATE "/opt/sootlib/include")

target_link_directories(flame.x PRIVATE "$ENV{CONDA_PREFIX}/lib")
target_link_directories(flame.x PRIVATE "/usr/local/lib")
target_link_directories(flame.x PRIVATE "/opt/radlib/lib")
target_link_directories(flame.x PRIVATE "/opt/sootlib/lib")

target_link_libraries(flame.x ${CANTERA} ${RADLIB} ${SOOTLIB} ${KINSOL} ${CVODES} ${ACCELERATE})
#=======
#find_library(YAML yaml-cpp)
#
#find_library(KINSOL sundials_kinsol)
#find_library(CVODES sundials_cvodes)
#
#find_package(LAPACK REQUIRED)
#
#target_include_directories(flame.x PRIVATE "/usr/local/include" "/opt/radlib/include" "/opt/sootlib/include/sootlib")
#target_link_directories(   flame.x PRIVATE "/usr/local/lib"     "/opt/radlib/lib"     "/opt/sootlib/lib")
#target_link_libraries(     flame.x ${CANTERA} ${YAML} ${KINSOL} ${CVODES} ${Python_LIBRARIES} ${Python} ${ACCELERATE} radlib sootModel)
#target_link_libraries(     flame.x LAPACK::LAPACK)
#>>>>>>> premixed

install(TARGETS flame.x DESTINATION "${PROJECT_SOURCE_DIR}/run")
