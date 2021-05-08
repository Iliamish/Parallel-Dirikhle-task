#include "Dirikhle.h"
#include "Markdown.h"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

//#define DEBUG

TEST_CASE("HEADER", "[classic]") {
	writeColumnsNames();
	REQUIRE("0" == "0");
}

#ifdef DEBUG

#else

TEST_CASE("Simple", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(501);
	clock_t begin_time = clock();
	simple_version(vec, 500, 500, 0, 1, 0, 1, 1000, 0.000001);
	std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Simple", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}

TEST_CASE("Sequential STL container std::sort 12", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(501);
	clock_t begin_time = clock();
	simple_tbb_version(vec, 500, 500, 0, 1, 0, 1, 1000, 0.000001);
	std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Sequential", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}

TEST_CASE("Sequential STL container std::sort 10", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(101);
	clock_t begin_time = clock();
	tbb_version_2(vec, 100, 100, 0, 1, 0, 1, 10000, 0.000001);
	std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Sequential", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}

TEST_CASE("Sequential STL container std::sort 121", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(101);
	clock_t begin_time = clock();
	tbb_version_3(vec, 100, 100, 0, 1, 0, 1, 10000, 0.000001);
	std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Sequential", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}

TEST_CASE("Sequential STL container std::sort 1221", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(101);
	clock_t begin_time = clock();
	tbb_version_4(vec, 100, 100, 0, 1, 0, 1, 10000, 0.000001);
	std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Sequential", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}


TEST_CASE("Sequential STL container std::sort 12211", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(101);
	clock_t begin_time = clock();
	tbb_version_5(vec, 100, 100, 0, 1, 0, 1, 10000, 0.000001);
	std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Sequential", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}

TEST_CASE("Sequential STL container std::sort 12d211", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(101);
	clock_t begin_time = clock();
	tbb_version_6(vec, 100, 100, 0, 1, 0, 1, 10000, 0.000001);
	std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Sequential", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}
#endif // DEBUG