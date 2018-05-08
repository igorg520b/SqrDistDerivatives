#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <sstream>
#include "point_triangle_distance.h"
using namespace std;
int main()
{
	// number of test cases to run
	const int nCases = 10;

	ofstream file_output("results.csv", ios::out);
	// prepare header line
	for (int k = 0; k < 12; k++) file_output << "x_" << k << ",";
	file_output << "sqDist,";
	for (int k = 0; k < 12; k++) file_output << "fd_" << k << ",";
	for (int k = 0; k < 12; k++) 
		for (int m = 0; m < 12; m++)
		file_output << "sd_" << k << "_" << m << ",";
	file_output << "zeta2,zeta3" << endl;
	file_output << std::scientific << std::setprecision(17);


	for (int i = 0; i < nCases; i++) {
		double x[12]; // p0, p1, p2, p3
		for (int k = 0; k < 12; k++) x[k] = 2 * (double)rand() / RAND_MAX - 1;

		// apply the function
		double fd[12] = {};
		double sd[12][12] = {};
		double s, t;
		int branch;

		// perform the calculation
		double dsq = pt(x, fd, sd, s, t, branch);

		// save to file
		for (int k = 0; k < 12; k++) file_output << x[k] << ",";
		file_output << dsq << ",";
		for (int k = 0; k < 12; k++) file_output << fd[k] << ",";
		for (int k = 0; k < 12; k++)
			for (int m = 0; m < 12; m++)
				file_output << sd[k][m] << ",";
		file_output << s << "," << t << endl;

	}
	file_output.close();

	return 0;
}

