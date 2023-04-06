#==============================================================================#
#                                                                              #
#        OpenSees -- Open System For Earthquake Engineering Simulation         #
#                Pacific Earthquake Engineering Research Center                #
#                                                                              #
#==============================================================================#
target_sources(OPS_Material
    PRIVATE
      MaterialTensorOperations.cpp
      MaterialStateVariables.cpp
      DataDrivenNestedSurfaces.cpp
	  YieldSurfacePackage.cpp
      MultiYieldSurfaceHardeningSoftening.cpp
      VonMisesDDM.cpp
      DruckerPragerDMM.cpp
      MatsuokaNakaiDMM.cpp
      OPS_VonMisesDMM.cpp
      OPS_DruckerPragerDMM.cpp
      OPS_MatsuokaNakaiDMM.cpp
    PUBLIC
      MaterialTensorOperations.h
      MaterialStateVariables.h
      DataDrivenNestedSurfaces.h
	  YieldSurfacePackage.h
      MultiYieldSurfaceHardeningSoftening.h
      VonMisesDMM.h
      DruckerPragerDMM.h
      MatsuokaNakaiDMM.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})