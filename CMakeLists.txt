cmake_minimum_required(VERSION 3.10)

project(rtc VERSION 0.0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(rtc_VERSION_MAJOR 0)
set(rtc_VERSION_MINOR 0)

add_library(librtc SHARED
	Curve.hpp
	Curve.cpp
)

set_target_properties(librtc PROPERTIES
	OUTPUT_NAME rtc
	VERSION 0.0.0
	SOVERSION 0
)

add_executable(rtc
	Curve.hpp
	Curve_test.cpp
)

target_link_libraries(rtc
	librtc
)
