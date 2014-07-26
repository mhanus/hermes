project(vrabec_E5-6T-diffusion)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp ../../../../utils.cpp)
add_executable(${PROJECT_NAME}-debug main.cpp definitions.cpp ../../../../utils.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${HERMES_RELEASE_FLAGS})
set_property(TARGET ${PROJECT_NAME}-debug PROPERTY COMPILE_FLAGS ${HERMES_DEBUG_FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D_RELEASE})
target_link_libraries(${PROJECT_NAME}-debug ${HERMES2D_DEBUG})

