cmake_minimum_required(VERSION 3.15)

# timestamp
# string(TIMESTAMP COMPILE_TIME %Y%m%d-%H%M%S)

# specify the c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# the project name and version
project(tJMF VERSION 1.0.0)
# configure_file(HubbardNU_MFConfig.h.in HubbardNu_MFConfig.h)

# Debug configuraton
# set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

# source code files
set(SRC_LIST main.cpp tJSBMF.cpp)

# add the self-constructed library
# add_subdirectory(MathFunctions)

# add executable
set(JSONCPP_PREFIX "/opt/homebrew/Cellar/jsoncpp/1.9.5")
# find_package(Eigen3 3.4 REQUIRED)
find_package(GSL REQUIRED)
# find_package(jsoncpp CONFIG REQUIRED)

add_executable(${PROJECT_NAME} ${SRC_LIST})
link_directories(${JSONCPP_PREFIX}/lib)
target_link_libraries(${PROJECT_NAME} PUBLIC GSL::gsl ${JSONCPP_PREFIX}/lib/libjsoncpp.dylib)

# add the binary tree to the search path for include files
# so that we will find HelloCMakeConfig.h
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR} ${GSL_INCLUDE_DIRS} ${JSONCPP_PREFIX}/include)