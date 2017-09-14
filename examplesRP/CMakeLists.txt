#----------------------------------------------------------------------------------------------
# Add include directories, sources and headers
#
include_directories(
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/inc
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ClicRP/inc
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ExN03RP/inc
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3RP/inc
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestFastSim/inc
  ${CMAKE_SOURCE_DIR}/realphysics/dummies/inc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/management/inc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/models/EM/inc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/processes/EM/inc
  ${CMAKE_SOURCE_DIR}/material/inc
  ${CMAKE_SOURCE_DIR}/base/inc
  ${CMAKE_SOURCE_DIR}/vecprot_v2/inc
  ${CMAKE_SOURCE_DIR}/numa/inc
)

if (USE_TBB)
  include_directories(${TBB_INCLUDE_DIR})
endif()

set(sources
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/src/CaloApp.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/src/CaloDetectorConstruction.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/src/CaloAppData.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/src/CaloPrimaryGenerator.cxx
#  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/src/UserPhysicsList.cxx
#  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/src/StepMaxProcess.cxx
#  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/src/Hist.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ExN03RP/src/ExN03ApplicationRP.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3RP/src/TestEM3ApplicationRP.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/src/TestEm5.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/src/TestEm5Data.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/src/UserDetectorConstruction.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/src/UserPrimaryGenerator.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/src/UserPhysicsList.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/src/StepMaxProcess.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/src/Hist.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestFastSim/src/TestFastSim.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestFastSim/src/TestFSData.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestFastSim/src/UserDetectorConstructionFS.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestFastSim/src/UserPrimaryGeneratorFS.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestFastSim/src/UserPhysicsListFS.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestFastSim/src/StepMaxProcessFS.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestFastSim/src/HistFS.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ClicRP/src/ClicApp.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ClicRP/src/ClicDetectorConstruction.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ClicRP/src/ClicAppData.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ClicRP/src/ClicPrimaryGenerator.cxx
#  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ClicRP/src/UserPhysicsList.cxx
#  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ClicRP/src/StepMaxProcess.cxx
#  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ClicRP/src/Hist.cxx
  )

set(headers
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/inc/CaloApp.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/inc/CaloDetectorConstruction.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/inc/CaloAppData.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/inc/CaloPrimaryGenerator.h
#  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/inc/UserPhysicsList.h
#  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/inc/StepMaxProcess.h
#  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/inc/Hist.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ExN03RP/inc/ExN03ApplicationRP.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3RP/inc/TestEM3ApplicationRP.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/TestEm5.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/TestEm5Data.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/UserDetectorConstruction.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/UserPrimaryGenerator.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/UserPhysicsList.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/StepMaxProcess.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/Hist.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestFastSim/inc/TestFastSim.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestFastSim/inc/TestFSData.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestFastSim/inc/UserDetectorConstructionFS.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestFastSim/inc/UserPrimaryGeneratorFS.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestFastSim/inc/UserPhysicsListFS.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestFastSim/inc/StepMaxProcessFS.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestFastSim/inc/HistFS.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ClicRP/inc/ClicApp.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ClicRP/inc/ClicDetectorConstruction.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ClicRP/inc/ClicAppData.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ClicRP/inc/ClicPrimaryGenerator.h
#  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ClicRP/inc/UserPhysicsList.h
#  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/inc/StepMaxProcess.h
#  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/inc/Hist.h
  )

if (USE_ROOT)
  ROOT_GENERATE_DICTIONARY(G__ExamplesRP ${headers} MODULE GeantExamplesRP
    LINKDEF
    ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ExN03RP/inc/LinkDef.h
##    ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3RP/inc/LinkDef.h
##    ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/LinkDef.h
    OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
    set(sources ${sources} G__ExamplesRP.cxx)
endif()
#----------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------
# Library and dependencies
#
add_library(GeantExamplesRP SHARED ${sources})

target_link_libraries(GeantExamplesRP -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v RealPhysics Material
  ${TBB_LIBRARIES})

if (USE_ROOT)
  target_link_libraries(GeantExamplesRP ${ROOT_LIBRARIES})
endif()

if (KNC)
  target_link_libraries(GeantExamplesRP ${Vc_MIC_LIBRARIES})
endif()

add_dependencies(GeantExamplesRP Geant_v RealPhysics Material)
add_custom_target(ExamplesRPlib DEPENDS GeantExamplesRP)

#----------------------------------------------------------------------------------------------
# Install the files
#
if(APPLE)
  # for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeantExamplesRP.so DESTINATION lib)
else()
  install(TARGETS GeantExamplesRP DESTINATION lib)
endif()

if(USE_ROOT)
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeantExamplesRP.rootmap
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeantExamplesRP_rdict.pcm DESTINATION lib)
endif()

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

#----------------------------------------------------------------------------------------------
# The Geant4 version of the example applications (only if build with Geant4 was requested)
#
if (WITH_GEANT4)
  add_subdirectory(Geant4/SimplyCalo)
  add_subdirectory(Geant4/TestEm5)
endif()
