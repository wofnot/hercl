
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11 -Wall")

project( hercl CXX )
find_package( OpenCV REQUIRED )
set(BIN "../bin")
set(INTER "../src/inter")
set(METRICS "../src/metrics")
set(BASE "../src/base")
set(SRC "../src")

set(needsCV art-school debug_art-school runMetrics)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN})

SET_SOURCE_FILES_PROPERTIES( ${SRC}/hashmap/hashmap.c ${BASE}/intersearch.c ${BASE}/interact.c ${BASE}/step.c ${BASE}/gradient.c ${BASE}/scan_print.c ${BASE}/point.c ${BASE}/cross_mutate.c ${BASE}/eval.c ${BASE}/message.c ${BASE}/debug_interact.c PROPERTIES LANGUAGE CXX )

add_executable( art-school ${SRC}/hashmap/hashmap.c ${INTER}/art-school.cpp ${INTER}/OCVCanvas.cpp ${BASE}/intersearch.c ${BASE}/interact.c ${BASE}/step.c ${BASE}/gradient.c ${BASE}/scan_print.c ${BASE}/point.c ${BASE}/cross_mutate.c ${BASE}/eval.c ${BASE}/message.c ${METRICS}/metrics.cpp)

add_executable( hm ${SRC}/hashmap/main.c ${SRC}/hashmap/hashmap.c)

add_executable( debug_art-school ${SRC}/hashmap/hashmap.c ${INTER}/art-school.cpp ${INTER}/OCVCanvas.cpp ${BASE}/debug_interact.c ${BASE}/interact.c ${BASE}/step.c ${BASE}/gradient.c ${BASE}/scan_print.c ${BASE}/point.c ${BASE}/cross_mutate.c ${BASE}/eval.c ${BASE}/message.c ../src/metrics/metrics.cpp )


add_executable( debug_mimic ${INTER}/mimic.c ${BASE}/debug_interact.c ${BASE}/interact.c ${BASE}/step.c ${BASE}/gradient.c ${BASE}/scan_print.c ${BASE}/point.c ${BASE}/cross_mutate.c ${BASE}/eval.c ${BASE}/message.c )

SET_TARGET_PROPERTIES(debug_art-school debug_mimic PROPERTIES COMPILE_FLAGS -DDEBUG)

add_executable( runMetrics ${METRICS}/runMetrics.cpp ${METRICS}/metrics.cpp  )

include_directories(${OpenCV_INCLUDE_DIRS})

foreach(exe ${needsCV})
  target_link_libraries( ${exe} ${OpenCV_LIBS} )
endforeach(exe)
