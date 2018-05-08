#pragma once
#include <cmath>
#include <algorithm>
#include <iostream>


// POINT-LINE
//  dot product of the form (p1-p0)(p2-p1)
double sp_dot3(double(&x)[9],	// input: coords of p0,p1,p2
	double(&fd)[9],				// output: first derivatives; must be cleared
	double(&sd)[9][9])			// output: second derivatives; must be cleared
{
	double x0 = x[0];
	double y0 = x[1];
	double z0 = x[2];

	double x1 = x[3];
	double y1 = x[4];
	double z1 = x[5];

	double x2 = x[6];
	double y2 = x[7];
	double z2 = x[8];

	fd[0] = x1 - x2;
	fd[1] = y1 - y2;
	fd[2] = z1 - z2;
	fd[3] = x0 - 2 * x1 + x2;
	fd[4] = y0 - 2 * y1 + y2;
	fd[5] = z0 - 2 * z1 + z2;
	fd[6] = x1 - x0;
	fd[7] = y1 - y0;
	fd[8] = z1 - z0;

	// second derivs
	for (int k = 0; k < 3; k++) {
		sd[k + 3][k + 3] = -2;
		sd[k][k + 6] = sd[k + 6][k] = -1;
		sd[k][k + 3] = sd[k + 3][k] = sd[k + 6][k + 3] = sd[k + 3][k + 6] = 1;
	}
	return (x1 - x0) * (x2 - x1) + (y1 - y0) * (y2 - y1) + (z1 - z0)* (z2 - z1);
}

// calculate the first and the second derivatives of f^2, given that f' and f'' are known
double function_squared(double f, // input: value at point
	double(&fd)[9],				// input: first and second derivatives
	double(&sd)[9][9],
	double(&fdOut)[9],				// output: first and second derivatives
	double(&sdOut)[9][9])
{
	for (int i = 0; i < 9; i++) {
		fdOut[i] = 2 * f*fd[i];
		for (int j = 0; j < 9; j++) sdOut[i][j] = 2 * (fd[i] * fd[j] + f*sd[i][j]);
	}
	return f*f;
}

double sp_dot3_squared(double(&x)[9],	// input: coords
	double(&fd)[9],				// output: first and second derivatives
	double(&sd)[9][9])
{
	double sp_dot3_fd[9] = {};
	double sp_dot3_sd[9][9] = {};

	double sp_dot3_value = sp_dot3(x, sp_dot3_fd, sp_dot3_sd);

	double result = function_squared(sp_dot3_value, sp_dot3_fd, sp_dot3_sd, fd, sd);
	return result;
}


// value, 1st and 2nd derivatives of the squared distance between points selected by idx1 and idx2
double vertex_vertex_distance_and_derivs(int idx1, int idx2,
	double(&x)[9],					// input: coords
	double(&sdd)[9],				// output: first derivatives
	double(&sdd2)[9][9]) {			// output: second derivatives

	int ix0 = idx1 * 3 + 0;
	int iy0 = idx1 * 3 + 1;
	int iz0 = idx1 * 3 + 2;

	int ix1 = idx2 * 3 + 0;
	int iy1 = idx2 * 3 + 1;
	int iz1 = idx2 * 3 + 2;

	double x0 = x[ix0];
	double y0 = x[iy0];
	double z0 = x[iz0];

	double x1 = x[ix1];
	double y1 = x[iy1];
	double z1 = x[iz1];

	sdd[ix0] = 2 * (x0 - x1);
	sdd[iy0] = 2 * (y0 - y1);
	sdd[iz0] = 2 * (z0 - z1);

	sdd[ix1] = -sdd[ix0];
	sdd[iy1] = -sdd[iy0];
	sdd[iz1] = -sdd[iz0];

	sdd2[ix0][ix0] = sdd2[iy0][iy0] = sdd2[iz0][iz0] = 2;
	sdd2[ix1][ix1] = sdd2[iy1][iy1] = sdd2[iz1][iz1] = 2;
	sdd2[ix0][ix1] = sdd2[iy0][iy1] = sdd2[iz0][iz1] = -2;
	sdd2[ix1][ix0] = sdd2[iy1][iy0] = sdd2[iz1][iz0] = -2;

	return (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1) + (z0 - z1)*(z0 - z1);
}

// value, 1st and 2nd derivatives of the squared distance between points selected by idx1 and idx2
double vertex_vertex_distance_and_derivs_12(int idx1, int idx2,
	double(&x)[12],					// input: coords
	double(&fd)[12],				// output: first derivatives
	double(&sd)[12][12]) {			// output: second derivatives

	int ix0 = idx1 * 3 + 0;
	int iy0 = idx1 * 3 + 1;
	int iz0 = idx1 * 3 + 2;

	int ix1 = idx2 * 3 + 0;
	int iy1 = idx2 * 3 + 1;
	int iz1 = idx2 * 3 + 2;

	double x0 = x[ix0];
	double y0 = x[iy0];
	double z0 = x[iz0];

	double x1 = x[ix1];
	double y1 = x[iy1];
	double z1 = x[iz1];

	fd[ix0] = 2 * (x0 - x1);
	fd[iy0] = 2 * (y0 - y1);
	fd[iz0] = 2 * (z0 - z1);

	fd[ix1] = -fd[ix0];
	fd[iy1] = -fd[iy0];
	fd[iz1] = -fd[iz0];

	sd[ix0][ix0] = sd[iy0][iy0] = sd[iz0][iz0] = 2;
	sd[ix1][ix1] = sd[iy1][iy1] = sd[iz1][iz1] = 2;
	sd[ix0][ix1] = sd[iy0][iy1] = sd[iz0][iz1] = -2;
	sd[ix1][ix0] = sd[iy1][iy0] = sd[iz1][iz0] = -2;

	return (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1) + (z0 - z1)*(z0 - z1);
}

double vertex_edge_distance_and_derivs(double(&x)[9],	// input coords; p0, line: (p1, p2)
	double(&sdd)[9],				// output: first derivatives
	double(&sdd2)[9][9]) {		// output: second derivatives

								// determine t, i.e. pc = p1(1-t)+p2(t)
	double x0 = x[0];
	double x1 = x[1];
	double x2 = x[2];
	double x3 = x[3];
	double x4 = x[4];
	double x5 = x[5];
	double x6 = x[6];
	double x7 = x[7];
	double x8 = x[8];

	double edge_length_sq = (-x3 + x6)*(-x3 + x6) + (-x4 + x7)*(-x4 + x7) + (-x5 + x8)*(-x5 + x8);
	double t = -((-x0 + x3)*(-x3 + x6) + (-x1 + x4)*(-x4 + x7) +
		(-x2 + x5)*(-x5 + x8)) / edge_length_sq;

	double u = 1 - t;
	double sqrDist = (x0 - u*x3 - t*x6)*(x0 - u*x3 - t*x6) +
		(x1 - u*x4 - t*x7)*(x1 - u*x4 - t*x7) +
		(x2 - u*x5 - t*x8)*(x2 - u*x5 - t*x8);
	double edge_length = sqrt(edge_length_sq);
	double dist = sqrt(sqrDist);
	double ratio = sqrDist / edge_length;


	/*
//	if (t < 0.01 && ratio < 0.001) {
		if (false) {
			// approximation 
		double fd_[9] = { 2 * (x0 - u*x3 - t*x6),2 * (x1 - u*x4 - t*x7),2 * (x2 - u*x5 - t*x8),-2 * u*(x0 - u*x3 - t*x6),-2 * u*(x1 - u*x4 - t*x7),-2 * u*(x2 - u*x5 - t*x8),-2 * t*(x0 - u*x3 - t*x6),-2 * t*(x1 - u*x4 - t*x7),-2 * t*(x2 - u*x5 - t*x8) };
		double sd_[9][9] = { { 2,0,0,-2 * u,0,0,-2 * t,0,0 },{ 0,2,0,0,-2 * u,0,0,-2 * t,0 },{ 0,0,2,0,0,-2 * u,0,0,-2 * t },{ -2 * u,0,0,2 * (u*u),0,0,2 * t*u,0,0 },{ 0,-2 * u,0,0,2 * (u*u),0,0,2 * t*u,0 },{ 0,0,-2 * u,0,0,2 * (u*u),0,0,2 * t*u },{ -2 * t,0,0,2 * t*u,0,0,2 * (t*t),0,0 },{ 0,-2 * t,0,0,2 * t*u,0,0,2 * (t*t),0 },{ 0,0,-2 * t,0,0,2 * t*u,0,0,2 * (t*t) } };

		for (int i = 0; i < 9; i++) {
			sdd[i] = fd_[i];
			for (int j = 0; j < 9; j++) {
				sdd2[i][j] = sd_[i][j];
			}
		}
	}
	*/

		double g;
		double g_fd[9] = {};
		double g_sd[9][9] = {};

		// |(x1-x0)(x2-x1)|^2 and its derivatives
		g = sp_dot3_squared(x, g_fd, g_sd);

		// f1 = |x1-x0|^2
		double f1;
		double f1fd[9] = {};
		double f1sd[9][9] = {};
		f1 = vertex_vertex_distance_and_derivs(1, 0, x, f1fd, f1sd);

		// f2 = |x2-x1|^2
		double f2;
		double f2fd[9] = {};
		double f2sd[9][9] = {};
		f2 = vertex_vertex_distance_and_derivs(2, 1, x, f2fd, f2sd);

		// combine together
		double f2sq = f2*f2;
		double f2cube = f2sq*f2;
		for (int i = 0; i < 9; i++) {
			sdd[i] = f1fd[i] + (g*f2fd[i] - f2*g_fd[i]) / f2sq;

			for (int j = 0; j < 9; j++) {
				double term1 = -2 * g*f2fd[i] * f2fd[j] / f2cube;
				double term2 = (g_fd[i] * f2fd[j] + g_fd[j] * f2fd[i]) / f2sq;
				double term3 = f1sd[i][j];
				double term4 = g*f2sd[i][j] / f2sq;
				double term5 = -g_sd[i][j] / f2;
				sdd2[i][j] = term1 + term2 + term3 + term4 + term5;
			}
		}

	return sqrDist;
}



// version for arrays with 12 elements
double vertex_edge_distance_and_derivs_12(double(&x)[12],	// input coords; p0, line: (p1, p2)
	int idx1, int idx2, // input indices for points p1 and p2
	double(&fd)[12],				// output: first derivatives
	double(&sd)[12][12]) {


	double _x[9];
	_x[0] = x[0];
	_x[1] = x[1];
	_x[2] = x[2];

	idx1 *= 3; idx2 *= 3;

	double p01 = (x[0] - x[0 + idx1])*(x[0] - x[0 + idx1]) + (x[1] - x[1 + idx1])*(x[1] - x[1 + idx1]) + (x[2] - x[2 + idx1])*(x[2] - x[2 + idx1]);
	double p02 = (x[0] - x[0 + idx2])*(x[0] - x[0 + idx2]) + (x[1] - x[1 + idx2])*(x[1] - x[1 + idx2]) + (x[2] - x[2 + idx2])*(x[2] - x[2 + idx2]);

	if (p01 > p02) {
		// swap indices
		int tmp_idx = idx1;
		idx1 = idx2;
		idx2 = tmp_idx;
	}

	_x[3] = x[0 + idx1];
	_x[4] = x[1 + idx1];
	_x[5] = x[2 + idx1];

	_x[6] = x[0 + idx2];
	_x[7] = x[1 + idx2];
	_x[8] = x[2 + idx2];

	double _fd[9] = {};
	double _sd[9][9] = {};

	double result = vertex_edge_distance_and_derivs(_x, _fd, _sd);

	// distribute _fd and _sd

	for (int i = 0; i < 3; i++) {
		fd[i] = _fd[i];
		fd[idx1 + i] = _fd[3 + i];
		fd[idx2 + i] = _fd[6 + i];

		for (int j = 0; j < 3; j++) {
			sd[i][j] = _sd[i][j];
			sd[i + idx1][j] = _sd[3 + i][j];
			sd[i][j + idx1] = _sd[i][3 + j];
			sd[i + idx1][j + idx1] = _sd[i + 3][j + 3];
			sd[i + idx1][j + idx2] = _sd[i + 3][j + 6];
			sd[i + idx2][j + idx1] = _sd[i + 6][j + 3];
			sd[i][j + idx2] = _sd[i][j + 6];
			sd[i + idx2][j] = _sd[i + 6][j];
			sd[i + idx2][j + idx2] = _sd[i + 6][j + 6];
		}
	}

	return result;
}