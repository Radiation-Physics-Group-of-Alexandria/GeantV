include_directories(inc)
add_executable(proc_mgr src/HPCProcessManager.cxx inc/HPCProcessManager.h)
target_link_libraries(proc_mgr ${cppzmq_LIBRARY})