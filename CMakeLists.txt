cmake_minimum_required(VERSION 3.23)
project(Trigonometric_Solution_LaPaz_1949)

set(CMAKE_CXX_STANDARD 23)


include_directories(${PROJECT_SOURCE_DIR}/excerpt)
link_directories(${PROJECT_SOURCE_DIR}/excerpt)
add_executable(${CMAKE_PROJECT_NAME} main.cpp ${PROJECT_SOURCE_DIR}/excerpt/excerpt.cpp)
