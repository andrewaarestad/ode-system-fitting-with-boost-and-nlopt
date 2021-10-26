macro(print_all_variables)
    message(STATUS "print_all_variables------------------------------------------{")
    get_cmake_property(_variableNames VARIABLES)
    foreach (_variableName ${_variableNames})
        message(STATUS "${_variableName}=${${_variableName}}")
    endforeach()
    message(STATUS "print_all_variables------------------------------------------}")
endmacro()

cmake_minimum_required(VERSION 3.20)
project(twoEquationModel)

set(CMAKE_CXX_STANDARD 14)
set(BOOST_ROOT boost_1_77_0)

find_package(Boost 1.77.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# ASSUMES NLOPT IS INSTALLED in /user/local/include, e.g. with `brew install nlopt`
include_directories(/usr/local/include)


add_executable(twoEquationModel main.cpp)
add_executable(harmonicExample harmonic_example.cpp)

target_link_libraries(twoEquationModel "/usr/local/lib/libnlopt.dylib")
target_link_libraries(twoEquationModel ${PYTHON_LIBRARIES})

print_all_variables()