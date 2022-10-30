cmake_minimum_required(VERSION 3.17)

project(fluidsim VERSION 1.0
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

############################################## OWN LIBRARIES ###########################################################


set(SOURCES
        "src/fieldsvis.cpp"
        "src/fluidsim.cpp"
        "src/matrixtexture.cpp"
        "src/ui.cpp"
        "src/uidata.cpp"
        )

set(LIBRARY_NAME engine)
add_library(${LIBRARY_NAME} OBJECT ${SOURCES})

target_include_directories(${LIBRARY_NAME} PUBLIC "${PROJECT_SOURCE_DIR}/include")


################################################ EXECUTABLE ############################################################


add_executable(main src/main.cpp)

target_link_libraries(main PRIVATE ${LIBRARY_NAME})


############################################ EXTERNAL LIBRARIES ########################################################

find_package(SFML COMPONENTS graphics window system REQUIRED)
target_link_libraries(${LIBRARY_NAME} PRIVATE
        sfml-graphics
        sfml-window
        sfml-system
        )

find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(main PRIVATE Eigen3::Eigen)
target_link_libraries(${LIBRARY_NAME} PRIVATE Eigen3::Eigen)

find_package(imgui CONFIG REQUIRED)
target_link_libraries(main PRIVATE imgui::imgui)
target_link_libraries(${LIBRARY_NAME} PRIVATE imgui::imgui)

find_package(ImGui-SFML CONFIG REQUIRED)
target_link_libraries(main PRIVATE ImGui-SFML::ImGui-SFML)
target_link_libraries(${LIBRARY_NAME} PRIVATE ImGui-SFML::ImGui-SFML)


find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(${LIBRARY_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()