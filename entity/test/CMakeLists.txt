message("Test build...")

set(ENTITY_INCLUDE "entity/")

set(ENTITY_TEST_HEADERS "entity/test/include/")
set(ENTITY_TEST_SOURCES "entity/test/lib/")

set(CXX_HEADERS )
list(TRANSFORM CXX_HEADERS PREPEND ${ENTITY_TEST_HEADERS})

set(CXX_SOURCES )
list(TRANSFORM CXX_SOURCES PREPEND ${ENTITY_TEST_SOURCES})

add_executable(entity_test entity/test/main.cpp ${CXX_HEADERS} ${CXX_SOURCES})

target_include_directories(entity_test PUBLIC ${ENTITY_INCLUDE})
target_link_libraries(entity_test entity)
