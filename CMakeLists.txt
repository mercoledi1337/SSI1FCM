cmake_minimum_required(VERSION 4.0)
project(SSI)

set(CMAKE_CXX_STANDARD 20)

# Vcpkg toolchain
set(CMAKE_TOOLCHAIN_FILE "C:/Users/mrser/Documents/VCPKG/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "" FORCE)

# Wyszukiwanie Boost
find_package(Boost REQUIRED CONFIG COMPONENTS iostreams filesystem system)

# Tworzymy exe
add_executable(SSI3 main.cpp gnuplot_iostream.h)

# Linkowanie bibliotek Boost
target_link_libraries(SSI3 PRIVATE Boost::iostreams Boost::filesystem Boost::system)
