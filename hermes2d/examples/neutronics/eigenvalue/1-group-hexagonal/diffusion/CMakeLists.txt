project(1-group-hexagonal-diffusion)
add_executable(${PROJECT_NAME} main.cpp definitions.cpp ../../../utils.cpp)
#set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})
target_link_libraries(${PROJECT_NAME} ${HERMES2D})

