message("Source build")

set(ENTITY_HEADERS "entity/src/include/")
set(ENTITY_SOURCES "entity/src/lib/")

set(CXX_HEADERS random_engine.h default_random_engine.h modular_ops.h ModRing128.h)
list(TRANSFORM CXX_HEADERS PREPEND ${ENTITY_HEADERS})

set(CXX_SOURCES default_random_engine.cpp)
list(TRANSFORM CXX_SOURCES PREPEND ${ENTITY_SOURCES})

message("HERE ${CXX_SOURCES}")

add_library(entity SHARED ${CXX_HEADERS} ${CXX_SOURCES})
target_include_directories(entity PUBLIC ${ENTITY_HEADERS})
