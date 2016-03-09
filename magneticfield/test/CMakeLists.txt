#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources )  # Common source  files
file(GLOB headers )  # Common include files
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS                "${GEANT5_CXX_FLAGS}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_DEBUG          "${GEANT5_CXX_FLAGS_DEBUG}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_MINISIZEREL    "${GEANT5_CXX_FLAGS_MINISIZEREL}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_RELEASE        "${GEANT5_CXX_FLAGS_RELEASE}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GEANT5_CXX_FLAGS_RELWITHDEBINFO}")
set_source_files_properties(${sources} PROPERTY CMAKE_EXEC_LINKER_FLAGS        "${GEANT5_EXE_LINKER_FLAGS}")
set_source_files_properties(${sources} PROPERTY COMPILE_FLAGS                  "${GEANT5_COMPILE_FLAGS}")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#

# include_directories(${PROJECT_SOURCE_DIR}/magneticfield/inc)
include_directories(${PROJECT_SOURCE_DIR}/magneticfield/inc ${PROJECT_SOURCE_DIR}/VecMagFieldRoutine/include)


## Need -lGeant_v in order to obtain the function/template Geant::cxx::ErrorHandlerImpl

add_executable(testMagFieldEquation testMagFieldEquation.cxx) ## ${sources} ${headers})
target_link_libraries(testMagFieldEquation Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Created target for testMagFieldEquation")

add_executable(testVectorMagFieldEquation testVectorMagFieldEquation.cxx) ## ${sources} ${headers})
target_link_libraries(testVectorMagFieldEquation Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#message(STATUS "Created target for testVectorMagFieldEquation")

add_executable(testTemplateMagFieldEquation testTemplateMagFieldEquation.cxx) ## ${sources} ${headers})
target_link_libraries(testTemplateMagFieldEquation Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

add_executable(testStepperFixedCashKarp testStepperFixedCashKarp.cxx) ## ${sources} ${headers})
target_link_libraries(testStepperFixedCashKarp Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Created target for testStepperFixedCashKarp")

add_executable(testTemplateStepperFixedCashKarp testTemplateStepperFixedCashKarp.cxx) ## ${sources} ${headers})
target_link_libraries(testTemplateStepperFixedCashKarp Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#message(STATUS "Created target for testTemplateStepperFixedCashKarp")

add_executable(testTemplateStepperFixedCashKarp2 testTemplateStepperFixedCashKarp2.cxx) ## ${sources} ${headers})
target_link_libraries(testTemplateStepperFixedCashKarp2 Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#message(STATUS "Created target for testTemplateStepperFixedCashKarp2")

add_executable(testIntegrationDriverCashKarp testIntegrationDriverCashKarp.cxx) ## ${sources} ${headers})
target_link_libraries(testIntegrationDriverCashKarp Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#message(STATUS "Created target for testIntegrationDriverCashKarp")


add_executable(testToyClass testToyClass.cxx) ## ${sources} ${headers})
target_link_libraries(testToyClass Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Created target for testToyClass")

add_executable(testVectorIntegrationDriver testVectorIntegrationDriver.cxx) ## ${sources} ${headers})
target_link_libraries(testVectorIntegrationDriver Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Created target for testVectorIntegrationDriver")

add_executable(BenchmarkIntegrationDriver BenchmarkIntegrationDriver.cxx) ## ${sources} ${headers})
target_link_libraries(BenchmarkIntegrationDriver Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Created target for BenchmarkIntegrationDriver")

#add_executable(TestConstBzStepper TestConstBzStepper.cxx) ## ${sources} ${headers})
#target_link_libraries(TestConstBzStepper Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})


add_executable(ReproduceError ReproduceError.cxx) ## ${sources} ${headers})
target_link_libraries(ReproduceError Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#message(STATUS "Created target for testVectorMagFieldEquation")

# add_executable(testTClassicalRK4 testTClassicalRK4.cxx) ## ${sources} ${headers})
# target_link_libraries(testTClassicalRK4 Vmagfield ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
# message(STATUS "Created target for testTClassicalRK4")

add_executable(testStepperFixed testStepperFixed.cxx) ## ${sources} ${headers})
target_link_libraries(testStepperFixed Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Created target for testStepperFixed")

add_executable(testIntegrationDriver testIntegrationDriver.cxx) ## ${sources} ${headers})
target_link_libraries(testIntegrationDriver Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Created target for testIntegrationDriver")

add_executable(testGUPropagateField testGUPropagateField.cxx) ## ${sources} ${headers})
target_link_libraries(testGUPropagateField Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Created target for testGUPropagateField")

add_executable(testConstVecFieldStepper testConstVecFieldStepper.cxx) ## ${sources} ${headers})
target_link_libraries(testConstVecFieldStepper Vmagfield -lGeant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Created target for testConstVecFieldStepper")

## add_executable(testFieldPropagator TestGUFieldPropagator.cxx) ## ${sources} ${headers})
## target_link_libraries(testFieldPropagator Vmagfield  ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
## # ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} -lGeom -lEG -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lXsec -lGeant_v)
## # add_dependencies(micromc xseclib Geantlib_v)
## message(STATUS "Created target for testFieldPropagator")

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
#add_custom_target(testFP DEPENDS testFieldPropagator)
#install(TARGETS testFieldPropagator DESTINATION bin)

add_custom_target(testMagEq DEPENDS testMagFieldEquation)
install(TARGETS testMagFieldEquation DESTINATION bin)
