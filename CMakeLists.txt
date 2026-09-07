cmake_minimum_required(VERSION 3.10)
project(ternaryCHNS_fig22D VERSION 1.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

function(add_fig22d_case target source output_directory)
    add_executable(${target} ${source})
    target_compile_definitions(
        ${target}
        PRIVATE OUTPUT_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/${output_directory}/\"
    )
    if(MSVC)
        target_compile_definitions(${target} PRIVATE _CRT_SECURE_NO_WARNINGS)
        target_compile_options(${target} PRIVATE /utf-8)
    else()
        target_link_libraries(${target} PRIVATE m)
    endif()
endfunction()

add_fig22d_case(ex2_theta30  chnsangle5b.c data5b)
add_fig22d_case(ex2_theta60  chnsangle6b.c data6b)
add_fig22d_case(ex2_theta90  chnsangle7b.c data7b)
add_fig22d_case(ex2_theta120 chnsangle8b.c data8b)
