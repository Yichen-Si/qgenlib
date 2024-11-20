cmake_minimum_required(VERSION 2.8...3.23)
project(qgenlib)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread -std=c++11")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")


find_path(HTS_INCLUDE_DIRS htslib/bgzf.h HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ ../htslib/)
if(NOT HTS_INCLUDE_DIRS )
    message(FATAL_ERROR "libhts HTS_INCLUDE_DIRS not found")
endif()

message(STATUS ${HTS_INCLUDE_DIRS})

find_library(HTS_LIBRARIES libhts.a hts HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ ../htslib/)
if(NOT HTS_LIBRARIES)
    message(FATAL_ERROR "libhts HTS_LIBRARIES not found")
endif()

message(STATUS ${HTS_LIBRARIES})

include_directories(${HTS_INCLUDE_DIRS})

set(SOURCE_FILES
    genome_interval.cpp
    gtf.cpp
    hts_utils.cpp
    params.cpp
    phred_helper.cpp
    qgen_error.cpp
    qgen_utils.cpp
    tsv_reader.cpp
    commands.cpp
    dataframe.cpp
    seq_utils.cpp
    qgenlib/dataframe.h
    qgenlib/commands.h
    qgenlib/genome_interval.h
    qgenlib/gtf.h
    qgenlib/hts_utils.h
    qgenlib/phred_helper.h
    qgenlib/qgen_error.h
    qgenlib/qgen_utils.h
    qgenlib/genome_loci.h
    qgenlib/pos_loci.h
    qgenlib/gtf_interval_tree.h
    qgenlib/params.h
    qgenlib/qgen_constant.h
    qgenlib/qgen_except.h
    qgenlib/tsv_reader.h
    qgenlib/seq_utils.h)

ADD_LIBRARY(qgen STATIC ${SOURCE_FILES})

SET(APP_EXE testqgen)

ADD_EXECUTABLE(${APP_EXE} test_qgen.cpp)

find_library(ZLIB z HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT ZLIB)
    message(FATAL_ERROR "libz library not found")
endif()

find_library(CURLLIB curl HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CURLLIB)
    message(FATAL_ERROR "libcurl library not found")
endif()

find_library(DEFLATELIB libdeflate.a deflate HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ ../libdeflate/)
if(NOT DEFLATELIB)
    message(FATAL_ERROR "libdeflate library not found")
endif()

find_package (BZip2)
if (NOT BZIP2_FOUND)
    message(FATAL_ERROR "libbz2 library not found")
else()
    include_directories(${BZIP2_INCLUDE_DIRS})
    set(BZIP2 ${BZIP2_LIBRARIES})
endif()

find_library(LZMA lzma HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT LZMA)
    message(WARNING "liblzma library not found, if you specified --disable-lzma when compiling libhts, please ignore this warning")
    set(LZMA "")
else()
    set(LZMA ${LZMA_LIBRARIES})
endif()

find_library(CRYPTOLIB crypto HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CRYPTOLIB)
    message(FATAL_ERROR "libcrypto library not found")
endif()

TARGET_LINK_LIBRARIES(${APP_EXE} qgen ${HTS_LIBRARIES} ${ZLIB} ${ZLIB} ${BZIP2} ${LZMA} ${CURLLIB} ${DEFLATELIB} ${CRYPTOLIB})

install(TARGETS ${APP_EXE} RUNTIME DESTINATION bin)
