# HEMCO/src/Interfaces/NUOPC/CMakeLists.txt

add_library(HCOI_NUOPC STATIC
        hcoi_nuopc_mod.F90
)

target_link_libraries(HCOI_NUOPC
    PUBLIC
        HCOI_Shared
        HEMCOBuildProperties
)

target_include_directories(HCOI_NUOPC
	INTERFACE ${HEMCO_BINARY_DIR}/mod
)

target_compile_options(HCOI_NUOPC
	PRIVATE
		""
		$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-ffree-line-length-none>
)