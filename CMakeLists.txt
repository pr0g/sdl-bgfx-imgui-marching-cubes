cmake_minimum_required(VERSION 3.15)

project(sdl-bgfx-imgui-marching-cubes LANGUAGES CXX)

find_package(SDL2 REQUIRED CONFIG)
find_package(bgfx REQUIRED CONFIG)
find_package(imgui REQUIRED CONFIG)

include(FetchContent)
FetchContent_Declare(
    as-camera-input
    GIT_REPOSITORY https://github.com/pr0g/as-camera-input.git
    GIT_TAG        777d7895874115bb81933f8e117414bf8efc9d29)
FetchContent_Declare(
    thh-bgfx-debug
    GIT_REPOSITORY https://github.com/pr0g/thh-bgfx-debug
    GIT_TAG        b0eb4ddeb9cfe19ee18986444eb437c08fb323fe)
FetchContent_MakeAvailable(as-camera-input thh-bgfx-debug)

add_executable(${PROJECT_NAME})
target_sources(
    ${PROJECT_NAME} PRIVATE
    main.cpp
    sdl-imgui/imgui_impl_sdl.cpp
    bgfx-imgui/imgui_impl_bgfx.cpp
    marching-cubes/marching-cubes.cpp)

target_link_libraries(
    ${PROJECT_NAME} PRIVATE
    SDL2::SDL2
    SDL2::SDL2main
    bgfx::bgfx
    imgui::imgui
    as-camera-input
    thh-bgfx-debug)

# required by as and as-camera-input library
target_compile_definitions(${PROJECT_NAME}
    PRIVATE $<$<BOOL:${AS_PRECISION_FLOAT}>:AS_PRECISION_FLOAT>
            $<$<BOOL:${AS_PRECISION_DOUBLE}>:AS_PRECISION_DOUBLE>
            $<$<BOOL:${AS_COL_MAJOR}>:AS_COL_MAJOR>
            $<$<BOOL:${AS_ROW_MAJOR}>:AS_ROW_MAJOR>)

target_compile_features(
    ${PROJECT_NAME} PRIVATE cxx_std_17)

# required for project when using visual studio generator
set_target_properties(
    ${PROJECT_NAME} PROPERTIES
    VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

if (WIN32)
# copy the SDL2.dll to the same folder as the executable
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    $<TARGET_FILE:SDL2::SDL2>
    $<TARGET_FILE_DIR:${PROJECT_NAME}>
    VERBATIM)
endif()
