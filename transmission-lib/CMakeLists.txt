
#-----------------------------------------------------------------------------
MACRO(BUILD_INSTALL_LIBRARY LIBRARY_NAME LIBRARY_TYPE LIBRARY_USE_GPU)
    ADD_LIBRARY(${LIBRARY_NAME} ${LIBRARY_TYPE} ${LIBRARY_NAME}.h ${LIBRARY_NAME}.cpp)
    IF(${LIBRARY_USE_GPU} MATCHES ON)
        TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${CUDA_LIBRARIES} ${CUDA_cufft_LIBRARY})
        CUDA_ADD_CUFFT_TO_TARGET(${LIBRARY_NAME})
        CUDA_ADD_CUBLAS_TO_TARGET(${LIBRARY_NAME})
    ELSE()
        #TARGET_LINK_LIBRARIES(${LIBRARY_NAME})
    ENDIF()
    INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION ${LIB_INSTALL_DIR})
    INSTALL(FILES ${LIBRARY_NAME}.h DESTINATION ${INC_INSTALL_DIR})
ENDMACRO(BUILD_INSTALL_LIBRARY)

#-----------------------------------------------------------------------------

BUILD_INSTALL_LIBRARY(_tt_common ${LIB_TYPE} ${USE_CUDA})
#-----------------------------------------------------------------------------
