
FILE(GLOB PY_FILES "*.py")
foreach(PY_FILE ${PY_FILES})
        configure_file(${PY_FILE} "${CMAKE_CURRENT_BINARY_DIR}"/${PY_FILE} @ONLY)
endforeach(PY_FILE)
EXECUTE_PROCESS(COMMAND python ${CMAKE_BINARY_DIR}/Python/NiftyRec/setup.py build)
         
#Pack with CMake
FILE(GLOB PY_FILES "${CMAKE_BINARY_DIR}/Python/NiftyRec/*.py")
INSTALL(FILES ${PY_FILES} DESTINATION ${PYTHON_INSTALL_DIR})

#Install when 'make install'
INSTALL(CODE "EXECUTE_PROCESS(COMMAND python ${CMAKE_BINARY_DIR}/Python/NiftyRec/setup.py install )")

