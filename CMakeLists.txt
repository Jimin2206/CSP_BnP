cmake_minimum_required(VERSION 3.15)
project(HighsExample)

set(CMAKE_CXX_STANDARD 17)

# HiGHS 설치 경로
set(HIGHS_DIR "C:/highs_test/HiGHS")

# 헤더 경로 추가
include_directories(
	${HIGHS_DIR}/highs
	${HIGHS_DIR}/build
)

# 라이브러리 경로 추가
link_directories(${HIGHS_DIR}/build/Debug/bin)

# 실행 파일 생성
#add_executable(HighsExample call_highs_from_cpp.cpp)
add_executable(CSP_BnP main.cpp "BnP.cpp" "BnP.h" "ProblemData.h" "ProblemData.cpp")

# HiGHS 라이브러리 연결
target_link_libraries(CSP_BnP highs)