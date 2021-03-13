set(ENTITY_HEADERS "entity/src/include/")
set(ENTITY_SOURCES "entity/src/lib/")

set(CXX_HEADERS basic_finite_fields.h FF128.h utils.h ntt.h interfaces.h entity.h)
list(TRANSFORM CXX_HEADERS PREPEND ${ENTITY_HEADERS})

set(CXX_SOURCES ntt.cpp)
list(TRANSFORM CXX_SOURCES PREPEND ${ENTITY_SOURCES})

add_library(entity SHARED ${CXX_HEADERS} ${CXX_SOURCES})
target_include_directories(entity PUBLIC ${ENTITY_HEADERS})
