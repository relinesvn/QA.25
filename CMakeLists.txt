cmake_minimum_required(VERSION 3.20)
project(ConsoleApplication1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1) Core-модуль з бізнес-логікою (testable)
add_library(core
    model.cpp
)

target_include_directories(core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# 2) Основний застосунок
add_executable(ConsoleApplication1
    ConsoleApplication1/ConsoleApplication1.cpp
)

target_link_libraries(ConsoleApplication1 PRIVATE core)
