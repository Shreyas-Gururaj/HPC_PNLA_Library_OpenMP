
# To add an executable named dev made out of develop.cpp and my_utilities.cpp
add_executable(dev develop.cpp my_utilities.cpp)

# Tells where headers of target are to be found
target_include_directories(dev PUBLIC {CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(dev PUBLIC ../include)

# Links first_app to library pnla
target_link_libraries(dev pnla)


