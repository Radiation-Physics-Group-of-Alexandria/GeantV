cmake_minimum_required(VERSION 2.8)
project(NudyTest)

#add_executable(fission fission.cpp)
#target_link_libraries(fission Nudy)

add_executable(inelastic inelastic.cpp)
target_link_libraries(inelastic Nudy)

#add_executable(thermal thermal.cpp)
#target_link_libraries(thermal Nudy)

#add_executable(thermal thermal.cpp)
#target_link_libraries(thermal Nudy)
