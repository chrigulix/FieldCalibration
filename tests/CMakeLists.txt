# Testing stuff. I use gtest for testing, for instruction on how to install or use it
# go here https://github.com/google/googletest and click around. A test template file
# is supplied too. Just add your new test to the TEST_FILES Variable.

include_directories(${GTEST_INCLUDE_DIRS} ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include)

# Add any further test files here
set(TEST_FILES
        test_utilities.cpp
        )

add_executable(RunTests ${TEST_FILES})
target_link_libraries(RunTests gtest fieldcal ${ROOT_LIBRARIES} -lpthread -lm)