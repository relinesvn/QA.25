cmake_minimum_required(VERSION 3.20)
project(ConsoleApplication1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- locate model.cpp/model.h (can be in repo root OR in ConsoleApplication1/) ---
set(MODEL_CPP "")
set(MODEL_H   "")

if (EXISTS "${CMAKE_SOURCE_DIR}/model.cpp")
    set(MODEL_CPP "${CMAKE_SOURCE_DIR}/model.cpp")
    set(MODEL_H   "${CMAKE_SOURCE_DIR}/model.h")
elseif (EXISTS "${CMAKE_SOURCE_DIR}/ConsoleApplication1/model.cpp")
    set(MODEL_CPP "${CMAKE_SOURCE_DIR}/ConsoleApplication1/model.cpp")
    set(MODEL_H   "${CMAKE_SOURCE_DIR}/ConsoleApplication1/model.h")
else()
    message(FATAL_ERROR
        "model.cpp not found. Put model.cpp/model.h either in project root or in ConsoleApplication1/ folder.")
endif()

# --- core library with business logic ---
add_library(core STATIC
    ${MODEL_CPP}
    ${MODEL_H}
)

target_include_directories(core PUBLIC
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/ConsoleApplication1
)

# --- app executable ---
add_executable(ConsoleApplication1
    ConsoleApplication1/ConsoleApplication1.cpp
)

target_link_libraries(ConsoleApplication1 PRIVATE core)
