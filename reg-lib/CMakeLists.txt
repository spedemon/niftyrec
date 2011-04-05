
#-----------------------------------------------------------------------------
SET(NAME _reg_affineTransformation)
ADD_LIBRARY(${NAME} ${LIB_TYPE} ${NAME}.h ${NAME}.cpp)
TARGET_LINK_LIBRARIES(${NAME} reg_nifti)
INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES ${NAME}.h DESTINATION ${INC_INSTALL_DIR})
#-----------------------------------------------------------------------------
SET(NAME _reg_resampling)
ADD_LIBRARY(${NAME} ${LIB_TYPE} ${NAME}.h ${NAME}.cpp)
TARGET_LINK_LIBRARIES(${NAME} reg_nifti _reg_affineTransformation)
INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES ${NAME}.h DESTINATION ${INC_INSTALL_DIR})
#-----------------------------------------------------------------------------
SET(NAME _reg_tools)
ADD_LIBRARY(${NAME} ${LIB_TYPE} ${NAME}.h ${NAME}.cpp)
TARGET_LINK_LIBRARIES(${NAME} reg_nifti)
INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES ${NAME}.h DESTINATION ${INC_INSTALL_DIR})
#-----------------------------------------------------------------------------
SET(NAME _reg_bspline)
ADD_LIBRARY(${NAME} ${LIB_TYPE} ${NAME}.h ${NAME}.cpp)
TARGET_LINK_LIBRARIES(${NAME} reg_nifti _reg_affineTransformation)
INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES ${NAME}.h DESTINATION ${INC_INSTALL_DIR})
#-----------------------------------------------------------------------------
SET(NAME _reg_mutualinformation)
ADD_LIBRARY(${NAME} ${LIB_TYPE} ${NAME}.h ${NAME}.cpp)
TARGET_LINK_LIBRARIES(${NAME} reg_nifti)
INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES ${NAME}.h DESTINATION ${INC_INSTALL_DIR})
#-----------------------------------------------------------------------------
SET(NAME _reg_ssd)
ADD_LIBRARY(${NAME} ${LIB_TYPE} ${NAME}.h ${NAME}.cpp)
TARGET_LINK_LIBRARIES(${NAME} reg_nifti)
INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES ${NAME}.h DESTINATION ${INC_INSTALL_DIR})
#-----------------------------------------------------------------------------
SET(NAME _reg_blockMatching)
ADD_LIBRARY(${NAME} ${LIB_TYPE} ${NAME}.h ${NAME}.cpp)
TARGET_LINK_LIBRARIES(${NAME} reg_nifti)
INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES ${NAME}.h DESTINATION ${INC_INSTALL_DIR})
#-----------------------------------------------------------------------------