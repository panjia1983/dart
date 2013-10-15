###############################################
# apps/impulseTest
file(GLOB impulseTest_srcs "*.cpp")
file(GLOB impulseTest_hdrs "*.h")
add_executable(impulseTest ${impulseTest_srcs} ${impulseTest_hdrs})
target_link_libraries(impulseTest dart)
set_target_properties(impulseTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
