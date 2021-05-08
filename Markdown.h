#pragma once
#include <fstream>
#include <string>

void writeColumnsNames() {
	std::ofstream outfile("Readme.md");
	outfile << "# Radix sort\n\n" << std::endl;
	outfile << "## Machine:\n" << std::endl;
	outfile << "	Processor: AMD Ryzen 7 2700 Eight-Core Processor 3.20 Gz\n" << std::endl;
	outfile << "	Memory: 32 Gb\n" << std::endl;
	outfile << "	System: Windows 10\n\n" << std::endl;

	outfile << "## Test results\n" << std::endl;

	outfile << "| Algorithm | Version | Container type | Number of elements | Digits in number | Sort time in seconds |" << std::endl;
	outfile << "| --- | --- | --- | --- | --- | --- |" << std::endl;
	outfile.close();
}

void writeResult(std::string testName, std::string testVersion, std::string containerName, int N, int digits, float elapsedTime) {
	std::ofstream outfile("Readme.md", std::ios::app);
	outfile << "| " << testName << " | " << testVersion << " | " << containerName << " | " << N << " | " << digits << " | " << elapsedTime << " | " << std::endl;
	outfile.close();
}