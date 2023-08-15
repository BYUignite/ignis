
add_executable(flame.x)
target_sources(flame.x PRIVATE driver.cc flame.cc flame.h solver_kinsol.h integrator_cvode.h)

if(APPLE)
    find_library(ACCELERATE Accelerate)
endif()
find_package(Threads REQUIRED)

target_include_directories(flame.x PRIVATE "/usr/local/include")
target_link_directories(   flame.x PRIVATE "/usr/local/lib")
target_link_libraries(     flame.x cantera yaml-cpp sundials_kinsol sundials_cvode sundials_cvodes sundials_ida sundials_nvecserial sundials_sunlinsoldense sundials_sunlinsolband ${ACCELERATE} Threads::Threads)

install(TARGETS flame.x DESTINATION "${PROJECT_SOURCE_DIR}/run")
