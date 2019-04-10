
find_package(ITK)
if(ITK_FOUND)
	include(${ITK_USE_FILE})
else(ITK_FOUND)
	message(FATAL_ERROR "Cannot build without ITK. Please set ITK_DIR.")
endif(ITK_FOUND)

add_executable(ImageCompareCommand ImageCompareCommand.cxx) 
target_link_libraries(ImageCompareCommand ${ITK_LIBRARIES})

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/testing/BTMCS_TEST)
file(COPY ${PROJECT_SOURCE_DIR}/testing/input/ForwardSolver.in DESTINATION ${PROJECT_BINARY_DIR}/testing/BTMCS_TEST)

add_test(NAME BTMCS_TEST COMMAND
    ForwardSolverDiffusion -gfileInput ${PROJECT_SOURCE_DIR}/testing/input/label_646464.img -T 89.8020514942181 -gdiffgm 1e-010 -gdiffwm 1e-009 -grho 0.03 -gp2 0 -gp1 1 -gzc 0.112713783647275 -gyc 0.130969875204521 -gxc 0.117607637039548 -log_summary -info -ksp_view -ksp_converged_reason WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/testing/BTMCS_TEST)

add_test(NAME BTMCS_CHECK1 COMMAND
	ImageCompareCommand -TestImage ${PROJECT_BINARY_DIR}/testing/BTMCS_TEST/TumorDensity.001.mhd -BaselineImage ${PROJECT_SOURCE_DIR}/testing/baseline/TumorDensity.001.mhd -ToleranceIntensity 1 -ToleranceNumberOfPixels 10)

add_test(NAME BTMCS_CHECK2 COMMAND
	ImageCompareCommand -TestImage ${PROJECT_BINARY_DIR}/testing/BTMCS_TEST/DeformationField.001.mhd -BaselineImage ${PROJECT_SOURCE_DIR}/testing/baseline/DeformationField.001.mhd -ToleranceIntensity 1 -ToleranceNumberOfPixels 10)
