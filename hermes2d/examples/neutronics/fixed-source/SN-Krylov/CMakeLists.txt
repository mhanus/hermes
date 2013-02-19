project(SN-Krylov)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp ../../utils.cpp)

set (FLAGS "${FLAGS} -g")
#set (FLAGS "${FLAGS} -O2")
set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})
