#include "Dirikhle.h"
#include "Markdown.h"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#define N_dim 500

//#define DEBUG

TEST_CASE("HEADER", "[classic]") {
	writeColumnsNames();
	REQUIRE("0" == "0");
}

#ifdef DEBUG

#else

TEST_CASE("Simple", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(N_dim + 1);
	clock_t begin_time = clock();
	int S = simple_version(vec, N_dim, N_dim, 0, 1, 0, 1, 10000, 0.0001);
	std::cout << "S: " << S << " time: " << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Simple", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}

TEST_CASE("Sequential STL container std::sort 12", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(N_dim + 1);
	clock_t begin_time = clock();
	int S = simple_tbb_version(vec, N_dim, N_dim, 0, 1, 0, 1, 10000, 0.0001);
	std::cout << "S: " << S << " time: " << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Sequential", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}

TEST_CASE("Sequential STL container std::sort 10", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(N_dim + 1);
	clock_t begin_time = clock();
	int S = tbb_version_2(vec, N_dim , N_dim , 0, 1, 0, 1, 10000, 0.0001);
	std::cout << "S: " << S << " time: " << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Sequential", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}

TEST_CASE("Sequential STL container std::sort 10as", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(N_dim + 1);
	clock_t begin_time = clock();
	int S = tbb_version_2_1(vec, N_dim, N_dim, 0, 1, 0, 1, 10000, 0.0001);
	std::cout << "S: " << S << " time: " << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Sequential", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}

TEST_CASE("Sequential STL container std::sort 121", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(N_dim + 1);
	clock_t begin_time = clock();
	int S = tbb_version_3(vec, N_dim , N_dim , 0, 1, 0, 1, 10000, 0.0001);
	std::cout << "S: " << S << " time: " << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Sequential", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}

TEST_CASE("Sequential STL container std::sort 1221", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(N_dim + 1);
	clock_t begin_time = clock();
	int S = tbb_version_4(vec, N_dim , N_dim , 0, 1, 0, 1, 10000, 0.0001);
	std::cout << "S: " << S << " time: " << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Sequential", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}


TEST_CASE("Sequential STL container std::sort 12211", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(N_dim + 1);
	clock_t begin_time = clock();
	int S = tbb_version_5(vec, N_dim , N_dim , 0, 1, 0, 1, 10000, 0.0001);
	std::cout << "S: " << S << " time: " << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Sequential", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}

TEST_CASE("Sequential STL container std::sort 122qq11", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(11);
	clock_t begin_time = clock();
	int S = tbb_version_6(vec, 10, 10, 0, 1, 0, 1, 10000, 0.0001);
	std::cout << "S: " << S << " time: " << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Sequential", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}


TEST_CASE("Sequential STL container std::sort 1ss2211", "[classic]") {
	const size_t N = 10000000;

	std::vector<std::vector<double>> vec(11);
	clock_t begin_time = clock();
	int S = tbb_version_8(vec, 10, 10, 0, 1, 0, 1, 10000, 0.0001);
	std::cout << "S: " << S << " time: " << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
	writeResult("std::sort", "Sequential", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);

	REQUIRE("0" == "0");
}

//TEST_CASE("Sequential STL container std::sort 12d211", "[classic]") {
//	const size_t N = 10000000;
//
//	std::vector<std::vector<double>> vec(N_dim + 1);
//	clock_t begin_time = clock();
//	int S = tbb_version_6(vec, N_dim , N_dim , 0, 1, 0, 1, 10000, 0.0001);
//	std::cout << "S: " << S << " time: " << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
//	writeResult("std::sort", "Sequential", "STL vector", N, 5, float(clock() - begin_time) / CLOCKS_PER_SEC);
//
//	REQUIRE("0" == "0");
//}
#endif // DEBUG