
FILE(GLOB PY_FILES "*.py")
foreach(PY_FILE ${PY_FILES})
        configure_file(${PY_FILE} "${CMAKE_CURRENT_BINARY_DIR}"/${PY_FILE} @COPYONLY)
endforeach(PY_FILE)
INSTALL(FILES ${PY_FILES} phantom.tif DESTINATION ${BIN_INSTALL_DIR}/simple_phantom_example)

