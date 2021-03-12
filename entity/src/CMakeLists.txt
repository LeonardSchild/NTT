set(ENTITY_HEADERS "entity/src/include/")
set(ENTITY_SOURCES "entity/src/lib/")

set(CXX_HEADERS modular_ops.h mod_ring_128.h utils.h dummy.h)
list(TRANSFORM CXX_HEADERS PREPEND ${ENTITY_HEADERS})

set(CXX_SOURCES dummy.cpp)
list(TRANSFORM CXX_SOURCES PREPEND ${ENTITY_SOURCES})

add_library(entity SHARED ${CXX_HEADERS} ${CXX_SOURCES})
target_include_directories(entity PUBLIC ${ENTITY_HEADERS})
