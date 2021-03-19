set(ENTITY_HEADERS "entity/src/include/")
set(ENTITY_SOURCES "entity/src/lib/")

set(CXX_HEADERS basic_finite_fields.h FF128.h utils.h basic_ntt_engine.h interfaces.h entity.h fast_ntt_engine.h cu_finite_field.cuh)
list(TRANSFORM CXX_HEADERS PREPEND ${ENTITY_HEADERS})

set(CXX_SOURCES basic_ntt_engine.cpp)
list(TRANSFORM CXX_SOURCES PREPEND ${ENTITY_SOURCES})

add_library(entity SHARED ${CXX_HEADERS} ${CXX_SOURCES})
target_include_directories(entity PUBLIC ${ENTITY_HEADERS})