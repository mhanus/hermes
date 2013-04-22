project(stankovski)
add_executable(${PROJECT_NAME} main.cpp ../../utils.cpp)
#add_executable(${PROJECT_NAME}-kelly main-kelly.cpp ../../utils.cpp)

#set_common_target_properties(${PROJECT_NAME} "HERMES2D")
#set_common_target_properties(${PROJECT_NAME}-kelly "HERMES2D")

set(FLAGS "${FLAGS} -O2")

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})
#set_property(TARGET ${PROJECT_NAME}-kelly PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})
#target_link_libraries(${PROJECT_NAME}-kelly ${HERMES2D})
