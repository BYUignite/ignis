
add_executable(flame.x)
target_sources(flame.x PRIVATE driver.cc flame.cc flame.h solver_kinsol.h integrator_cvode.h streams.cc streams.h)

if(APPLE)
    find_library(ACCELERATE Accelerate)
    find_library(YAML yaml-cpp)
endif()
find_package(Threads REQUIRED)

find_library(KINSOL sundials_kinsol)
find_library(CVODES sundials_cvodes)
find_library(CANTERA cantera)
find_package(Python COMPONENTS Interpreter Development)

target_include_directories(flame.x PRIVATE "/usr/local/include" "/opt/radlib/include")
target_link_directories(   flame.x PRIVATE "/usr/local/lib" "/opt/radlib/lib")
target_link_libraries(     flame.x ${CANTERA} ${YAML} ${KINSOL} ${CVODES} ${Python_LIBRARIES} ${Python} ${ACCELERATE} radlib)

install(TARGETS flame.x DESTINATION "${PROJECT_SOURCE_DIR}/run")
