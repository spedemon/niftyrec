
MACRO(BUILD_INSTALL_LIBRARY LIBRARY_NAME LIBRARY_TYPE)
        CUDA_ADD_LIBRARY(${LIBRARY_NAME} ${LIB_TYPE} ${LIBRARY_NAME}.h ${LIBRARY_NAME}.cu ${LIBRARY_NAME}_kernels.cu)
        #TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${CUDA_SDK_LIBRARIES} ${CUDA_CUT_LIBRARY} )
        INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION ${LIB_INSTALL_DIR})
        INSTALL(FILES ${LIBRARY_NAME}.h DESTINATION ${INC_INSTALL_DIR})
ENDMACRO(BUILD_INSTALL_LIBRARY)

#-----------------------------------------------------------------------------
IF(USE_CUDA)
        BUILD_INSTALL_LIBRARY(tt_project_ray_gpu ${LIB_TYPE})
        BUILD_INSTALL_LIBRARY(tt_backproject_ray_gpu ${LIB_TYPE})
ENDIF(USE_CUDA)

#-----------------------------------------------------------------------------

