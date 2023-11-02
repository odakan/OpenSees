#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
    PRIVATE
      CTensor.cpp
	  MaterialStateVariables.cpp
	  NestedSurface.cpp
	  SRSMYSand.cpp
	  OPS_SRSMYSand.cpp
    PUBLIC
      CTensor.h
	  MaterialStateVariables.h
	  NestedSurface.h
	  SRSMYSand.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
