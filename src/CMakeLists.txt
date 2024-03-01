# Set an environment variable CMAKE_LIBRARY_PATH with paths for find_library to search
# or add: HINTS "/path/to/lib/dir" at the end of the find_library call

add_executable(flame.x)
target_sources(flame.x PRIVATE driver.cc flame.cc flame.h solver_kinsol.h integrator_cvode.h streams.cc streams.h)

if(APPLE)
    find_library(ACCELERATE Accelerate)
    set_property(TARGET flame.x PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)  # adds LC_RPATH to executable
endif()

find_library(CANTERA NAMES cantera_shared cantera)
find_library(RADLIB  NAMES radlib_shared  radlib)
find_library(KINSOL sundials_kinsol)
find_library(CVODES sundials_cvodes)

target_include_directories(flame.x PRIVATE "$ENV{CONDA_PREFIX}/include")
target_include_directories(flame.x PRIVATE "/usr/local/include")
target_include_directories(flame.x PRIVATE "/opt/radlib/include")

target_link_directories(flame.x PRIVATE "$ENV{CONDA_PREFIX}/lib")
target_link_directories(flame.x PRIVATE "/usr/local/lib")
target_link_directories(flame.x PRIVATE "/opt/radlib/lib")

target_link_libraries(flame.x ${CANTERA} ${RADLIB} ${KINSOL} ${CVODES} ${ACCELERATE})

install(TARGETS flame.x DESTINATION "${PROJECT_SOURCE_DIR}/run")
