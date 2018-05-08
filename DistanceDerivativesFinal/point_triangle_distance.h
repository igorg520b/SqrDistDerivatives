#pragma once
#include <cmath>
#include <algorithm>
#include <iostream>
#include"point_line.h"

// POINT-PLANE (intended to use in the interior of the triangle)
bool make_pause = false;
double dist_side_ratio;
// second derivatives of a
double a2[12][12] = {
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0 },
	{ 0, 0, 0, -2, 0, 0, 2, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, -2, 0, 0, 2, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, -2, 0, 0, 2, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };

// second derivatives of b
double b2[12][12] = {
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 2, 0, 0, -1, 0, 0, -1, 0, 0 },
	{ 0, 0, 0, 0, 2, 0, 0, -1, 0, 0, -1, 0 },
	{ 0, 0, 0, 0, 0, 2, 0, 0, -1, 0, 0, -1 },
	{ 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0 },
	{ 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0 },
	{ 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1 },
	{ 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0 } };

// second derivatives of c
double c2[12][12] = {
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 2, 0, 0, 0, 0, 0, -2, 0, 0 },
	{ 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, -2, 0 },
	{ 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, -2 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, -2, 0, 0, 0, 0, 0, 2, 0, 0 },
	{ 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 2, 0 },
	{ 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 2 } };

// second derivatives of d
double d2[12][12] = {
	{ 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0 },
	{ 1, 0, 0, -2, 0, 0, 1, 0, 0, 0, 0, 0 },
	{ 0, 1, 0, 0, -2, 0, 0, 1, 0, 0, 0, 0 },
	{ 0, 0, 1, 0, 0, -2, 0, 0, 1, 0, 0, 0 },
	{ -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };

// second derivatives of e
double e2[12][12] = {
	{ 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0 },
	{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0 },
	{ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1 },
	{ 1, 0, 0, -2, 0, 0, 0, 0, 0, 1, 0, 0 },
	{ 0, 1, 0, 0, -2, 0, 0, 0, 0, 0, 1, 0 },
	{ 0, 0, 1, 0, 0, -2, 0, 0, 0, 0, 0, 1 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0 } };

// second derivatives of f
double f2[12][12] = {
	{ 2, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 2, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 2, 0, 0, -2, 0, 0, 0, 0, 0, 0 },
	{ -2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, -2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, -2, 0, 0, 2, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };


// helper functions

void cross(double x0, double x1, double x2, double x3, double x4, double x5,
	double &out_v0, double &out_v1, double &out_v2) {
	out_v0 = -x2 * x4 + x1 * x5;
	out_v1 = x2 * x3 - x0 * x5;
	out_v2 = -x1 * x3 + x0 * x4;
}

void normalize(double &x, double &y, double &z) {
	double mag = sqrt(x*x + y*y + z*z);
	x /= mag;
	y /= mag;
	z /= mag;
}

void plane_normal(double x0, double x1, double x2,
	double x3, double x4, double x5,
	double x6, double x7, double x8,
	double &out_v0, double &out_v1, double &out_v2) {

	cross(x3 - x0, x4 - x1, x5 - x2,
		x6 - x0, x7 - x1, x8 - x2,
		out_v0, out_v1, out_v2);
	normalize(out_v0, out_v1, out_v2);
}

double dot(double x0, double x1, double x2,
	double x3, double x4, double x5) {
	return x0*x3 + x1*x4 + x2*x5;
}

// squared distance between two points
double vertex_vertex_distance(int idx1, int idx2, double(&x)[12]) {

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

	return (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1) + (z0 - z1)*(z0 - z1);
}

// POINT-PLANE universal (automatic selection between degenerate and non-degenerate)
bool point_plane_approximation_used;

double xd(int idx1, int idx2) {
	return idx1 == idx2 ? 1. : 0;
}

double point_plane_distance(double(&x)[12],
	double(&fd)[12], double(&sd)[12][12])
{
	double output_s, output_t; // for testing
	double x0 = x[0];
	double x1 = x[1];
	double x2 = x[2];
	double x3 = x[3];
	double x4 = x[4];
	double x5 = x[5];
	double x6 = x[6];
	double x7 = x[7];
	double x8 = x[8];
	double x9 = x[9];
	double x10 = x[10];
	double x11 = x[11];

	double abcdef[6] = { (-x3 + x6)*(-x3 + x6) + (-x4 + x7)*(-x4 + x7) + (-x5 + x8)*(-x5 + x8),(x10 - x4)*(-x4 + x7) + (x11 - x5)*(-x5 + x8) + (-x3 + x6)*(-x3 + x9),(x10 - x4)*(x10 - x4) + (x11 - x5)*(x11 - x5) + (-x3 + x9)*(-x3 + x9),(-x0 + x3)*(-x3 + x6) + (-x1 + x4)*(-x4 + x7) + (-x2 + x5)*(-x5 + x8),(x10 - x4)*(-x1 + x4) + (x11 - x5)*(-x2 + x5) + (-x0 + x3)*(-x3 + x9),(-x0 + x3)*(-x0 + x3) + (-x1 + x4)*(-x1 + x4) + (-x2 + x5)*(-x2 + x5) };

	double a = abcdef[0];
	double b = abcdef[1];
	double c = abcdef[2];
	double d = abcdef[3];
	double e = abcdef[4];
	double f = abcdef[5];

	double det = a*c - b*b;
	double detsq = det * det;
	double detcube = detsq * det;

	double s = b*e - c*d;
	double t = b*d - a*e;

	double invDet = 1. / det;
	s *= invDet;
	t *= invDet;
	output_s = s;
	output_t = t;
	double u = 1 - (s + t);

	double sqrDistance = (-x0 + x6*s + x9*t + x3*u)*(-x0 + x6*s + x9*t + x3*u) +
		(-x1 + x7*s + x10*t + x4*u)*(-x1 + x7*s + x10*t + x4*u) +
		(-x2 + x8*s + x11*t + x5*u)*(-x2 + x8*s + x11*t + x5*u);
	double dist = sqrt(sqrDistance);

	// select either normal case or degenerate case

	double p12 = sqrt(vertex_vertex_distance(1, 2, x));
	double p13 = sqrt(vertex_vertex_distance(1, 3, x));
	double p23 = sqrt(vertex_vertex_distance(2, 3, x));
	double smallest_side = std::min(std::min(p12, p13), p23);
	double max_ust = std::max(std::max(u, s), t);
	bool use_approximation = (dist / smallest_side < 1e-4 && max_ust > 0.999);
	use_approximation = false;

	dist_side_ratio = dist / smallest_side;

	double s2[12][12], t2[12][12], det2[12][12];

	double u1[12], u2[12][12]; // experimental/testing

	// derivatives of s and t
	// first derivatives of the above quantities
	double a1[12] = { 0,0,0,-2 * (-x3 + x6),-2 * (-x4 + x7),-2 * (-x5 + x8),2 * (-x3 + x6),2 * (-x4 + x7),2 * (-x5 + x8),0,0,0 };
	double b1[12] = { 0,0,0,2 * x3 - x6 - x9,-x10 + 2 * x4 - x7,-x11 + 2 * x5 - x8,-x3 + x9,x10 - x4,x11 - x5,-x3 + x6,-x4 + x7,-x5 + x8 };
	double c1[12] = { 0,0,0,-2 * (-x3 + x9),-2 * (x10 - x4),-2 * (x11 - x5),0,0,0,2 * (-x3 + x9),2 * (x10 - x4),2 * (x11 - x5) };
	double d1[12] = { x3 - x6, x4 - x7, x5 - x8, x0 - 2 * x3 + x6, x1 - 2 * x4 + x7, x2 - 2 * x5 + x8, -x0 + x3, -x1 + x4, -x2 + x5, 0, 0, 0 };
	double e1[12] = { x3 - x9, -x10 + x4, -x11 + x5, x0 - 2 * x3 + x9, x1 + x10 - 2 * x4, x11 + x2 - 2 * x5, 0, 0, 0, -x0 + x3, -x1 + x4, -x2 + x5 };
	double f1[12] = { -2 * (-x0 + x3),-2 * (-x1 + x4),-2 * (-x2 + x5),2 * (-x0 + x3),2 * (-x1 + x4),2 * (-x2 + x5),0,0,0,0,0,0 };
	double s1[12], t1[12], det1[12];

	// first derivatives 
	for (int i = 0; i < 12; i++)
	{
		det1[i] = c*a1[i] + a*c1[i] - 2 * b*b1[i];
		s1[i] = ((c*d - b*e)*det1[i]) / detsq + ((e*b1[i] + b*e1[i]) - (d*c1[i] + c*d1[i])) / det;
		t1[i] = ((a*e - b*d)*det1[i]) / detsq + ((d*b1[i] + b*d1[i]) - (a*e1[i] + e*a1[i])) / det;
		u1[i] = -(s1[i] + t1[i]);

		fd[i] = -2 * (x0 - x6*s - x9*t - x3*u)*
			(x6*s1[i] + x9*t1[i] + x3*u1[i] - xd(0, i) + u*xd(3, i) +
				s*xd(6, i) + t*xd(9, i)) -
			2 * (x1 - x7*s - x10*t - x4*u)*
			(x7*s1[i] + x10*t1[i] + x4*u1[i] - xd(1, i) + u*xd(4, i) +
				s*xd(7, i) + t*xd(10, i)) -
			2 * (x2 - x8*s - x11*t - x5*u)*
			(x8*s1[i] + x11*t1[i] + x5*u1[i] - xd(2, i) + u*xd(5, i) +
				s*xd(8, i) + t*xd(11, i));
	}

	/*
	if (use_approximation) {
		point_plane_approximation_used = true;
		double fd_[12] = { 2 * (x0 - u*x3 - s*x6 - t*x9),2 * (x1 - t*x10 - u*x4 - s*x7),2 * (-(t*x11) + x2 - u*x5 - s*x8),-2 * u*(x0 - u*x3 - s*x6 - t*x9),-2 * u*(x1 - t*x10 - u*x4 - s*x7),-2 * u*(-(t*x11) + x2 - u*x5 - s*x8),-2 * s*(x0 - u*x3 - s*x6 - t*x9),-2 * s*(x1 - t*x10 - u*x4 - s*x7),-2 * s*(-(t*x11) + x2 - u*x5 - s*x8),-2 * t*(x0 - u*x3 - s*x6 - t*x9),-2 * t*(x1 - t*x10 - u*x4 - s*x7),-2 * t*(-(t*x11) + x2 - u*x5 - s*x8) };
		for (int i = 0; i < 12; i++) fd[i] = fd_[i];
	}
	*/

	for (int i = 0; i < 12; i++)
		for (int j = 0; j < 12; j++)
		{
			det2[i][j] = -2 * b1[i] * b1[j] + a1[j] * c1[i] + a1[i] * c1[j] + c*a2[i][j] - 2 * b*b2[i][j] + a*c2[i][j];

			s2[i][j] =
				+(-(c1[j] * d1[i]) - c1[i] * d1[j] + b1[j] * e1[i] + b1[i] * e1[j] + e*b2[i][j] - d*c2[i][j] - c*d2[i][j] + b*e2[i][j]) / det
				- ((det1[j] * (e*b1[i] - d*c1[i] - c*d1[i] + b*e1[i])) + (det1[i] * (e*b1[j] - d*c1[j] - c*d1[j] + b*e1[j])) + ((-(c*d) + b*e)*det2[i][j])) / detsq
				+ (2 * (-(c*d) + b*e)*det1[i] * det1[j]) / detcube;

			t2[i][j] =
				+(b1[j] * d1[i] + b1[i] * d1[j] - a1[j] * e1[i] - a1[i] * e1[j] - e*a2[i][j] + d*b2[i][j] + b*d2[i][j] - a*e2[i][j]) / det
				- ((det1[j] * (-(e*a1[i]) + d*b1[i] + b*d1[i] - a*e1[i])) + (det1[i] * (-(e*a1[j]) + d*b1[j] + b*d1[j] - a*e1[j])) + ((b*d - a*e)*det2[i][j])) / detsq
				+ (2 * (b*d - a*e)*det1[i] * det1[j]) / detcube;

			u2[i][j] = -(s2[i][j] + t2[i][j]);

			sd[i][j] = 2 * ((x6*s1[i] + x9*t1[i] + x3*u1[i] - xd(0, i) + u*xd(3, i) +
				s*xd(6, i) + t*xd(9, i))*
				(x6*s1[j] + x9*t1[j] + x3*u1[j] - xd(0, j) + u*xd(3, j) +
					s*xd(6, j) + t*xd(9, j)) -
					(x0 - x6*s - x9*t - x3*u)*
				(-0 + 0 * s + 0 * t + 0 * u + x6*s2[i][j] + x9*t2[i][j] +
					x3*u2[i][j] + u1[j] * xd(3, i) + u1[i] * xd(3, j) + s1[j] * xd(6, i) +
					s1[i] * xd(6, j) + t1[j] * xd(9, i) + t1[i] * xd(9, j)) +
					(x7*s1[i] + x10*t1[i] + x4*u1[i] - xd(1, i) + u*xd(4, i) +
						s*xd(7, i) + t*xd(10, i))*
						(x7*s1[j] + x10*t1[j] + x4*u1[j] - xd(1, j) + u*xd(4, j) +
							s*xd(7, j) + t*xd(10, j)) -
							(x1 - x7*s - x10*t - x4*u)*
				(-0 + 0 * s + 0 * t + 0 * u + x7*s2[i][j] + x10*t2[i][j] +
					x4*u2[i][j] + u1[j] * xd(4, i) + u1[i] * xd(4, j) + s1[j] * xd(7, i) +
					s1[i] * xd(7, j) + t1[j] * xd(10, i) + t1[i] * xd(10, j)) +
					(x8*s1[i] + x11*t1[i] + x5*u1[i] - xd(2, i) + u*xd(5, i) +
						s*xd(8, i) + t*xd(11, i))*
						(x8*s1[j] + x11*t1[j] + x5*u1[j] - xd(2, j) + u*xd(5, j) +
							s*xd(8, j) + t*xd(11, j)) -
							(x2 - x8*s - x11*t - x5*u)*
				(-0 + 0 * s + 0 * t + 0 * u + x8*s2[i][j] + x11*t2[i][j] +
					x5*u2[i][j] + u1[j] * xd(5, i) + u1[i] * xd(5, j) + s1[j] * xd(8, i) +
					s1[i] * xd(8, j) + t1[j] * xd(11, i) + t1[i] * xd(11, j)));
		}

	/*
	if (make_pause) {
		std::cout << "*";
	}
	*/
	return sqrDistance;
}



// POINT-TRIANGE (ALL CASES COMBINED)

// point-triangle derivatives; returns squared distance
double pt(double(&x)[12],
	double(&fd)[12], double(&sd)[12][12], double &output_s, double &output_t, int &branch)
{
	double x0 = x[0];
	double x1 = x[1];
	double x2 = x[2];
	double x3 = x[3];
	double x4 = x[4];
	double x5 = x[5];
	double x6 = x[6];
	double x7 = x[7];
	double x8 = x[8];
	double x9 = x[9];
	double x10 = x[10];
	double x11 = x[11];

	double abcdef[6] = { (-x3 + x6)*(-x3 + x6) + (-x4 + x7)*(-x4 + x7) + (-x5 + x8)*(-x5 + x8),
		(x10 - x4)*(-x4 + x7) + (x11 - x5)*(-x5 + x8) + (-x3 + x6)*(-x3 + x9),
		(x10 - x4)*(x10 - x4) + (x11 - x5)*(x11 - x5) + (-x3 + x9)*(-x3 + x9),
		(-x0 + x3)*(-x3 + x6) + (-x1 + x4)*(-x4 + x7) + (-x2 + x5)*(-x5 + x8),
		(x10 - x4)*(-x1 + x4) + (x11 - x5)*(-x2 + x5) + (-x0 + x3)*(-x3 + x9),
		(-x0 + x3)*(-x0 + x3) + (-x1 + x4)*(-x1 + x4) + (-x2 + x5)*(-x2 + x5) };

	double a = abcdef[0];
	double b = abcdef[1];
	double c = abcdef[2];
	double d = abcdef[3];
	double e = abcdef[4];
	double f = abcdef[5];

	double det = a*c - b*b;
//	double detsq = det * det;
//	double detcube = detsq * det;

	double s = b*e - c*d;
	double t = b*d - a*e;

	branch = -1;
	double sqrDistance;


	if (s + t <= det)
	{
		if (s < 0)
		{
			if (t < 0)  // region 4
			{
				if (d < 0)
				{
					t = 0;
					if (-d >= a)
					{
						branch = 1;
						s = 1;
						sqrDistance = vertex_vertex_distance_and_derivs_12(0, 2, x, fd, sd);
					}
					else {
						branch = 5;
						s = -d / a;
						sqrDistance = vertex_edge_distance_and_derivs_12(x, 2, 1, fd, sd);
					}
				}
				else {
					s = 0;
					if (e >= 0)
					{
						branch = 3;
						t = 0;
						sqrDistance = vertex_vertex_distance_and_derivs_12(0, 1, x, fd, sd);
					}
					else if (-e >= c)
					{
						branch = 2;
						t = 1;
						sqrDistance = vertex_vertex_distance_and_derivs_12(0, 3, x, fd, sd);
					}
					else {
						branch = 4;
						t = -e / c;
						sqrDistance = vertex_edge_distance_and_derivs_12(x, 3, 1, fd, sd);
					}
				}
			}
			else  // region 3
			{
				s = 0;
				if (e >= 0)
				{
					branch = 3;
					t = 0;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 1, x, fd, sd);
				}
				else if (-e >= c)
				{
					branch = 2;
					t = 1;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 3, x, fd, sd);
				}
				else {
					branch = 4;
					t = -e / c;
					sqrDistance = vertex_edge_distance_and_derivs_12(x, 3, 1, fd, sd);
				}
			}
		}
		else if (t < 0)  // region 5
		{
			t = 0;
			if (d >= 0)
			{
				branch = 3;
				s = 0;
				sqrDistance = vertex_vertex_distance_and_derivs_12(0, 1, x, fd, sd);
			}
			else if (-d >= a)
			{
				branch = 1;
				s = 1;
				sqrDistance = vertex_vertex_distance_and_derivs_12(0, 2, x, fd, sd);
			}
			else {
				branch = 5;
				s = -d / a;
				sqrDistance = vertex_edge_distance_and_derivs_12(x, 1, 2, fd, sd);
			}
		}
		else  // region 0
		{
			branch = 0; // interior point
			double invDet = (1) / det;
			s *= invDet;
			t *= invDet;
			sqrDistance = point_plane_distance(x, fd, sd);
		}
	}
	else {
		double tmp0, tmp1, numer, denom;

		if (s < 0)  // region 2
		{
			tmp0 = b + d;
			tmp1 = c + e;
			if (tmp1 > tmp0)
			{
				numer = tmp1 - tmp0;
				denom = a - 2*b + c;
				if (numer >= denom)
				{
					branch = 1;
					s = 1;
					t = 0;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 2, x, fd, sd);
				}
				else {
					branch = 6;
					s = numer / denom;
					t = 1 - s;
					sqrDistance = vertex_edge_distance_and_derivs_12(x, 2, 3, fd, sd);
				}
			}
			else {
				s = 0;
				if (tmp1 <= 0)
				{
					branch = 2;
					t = 1;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 3, x, fd, sd);
				}
				else if (e >= 0)
				{
					branch = 3;
					t = 0;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 1, x, fd, sd);
				}
				else {
					branch = 4;
					t = -e / c;
					sqrDistance = vertex_edge_distance_and_derivs_12(x, 1, 3, fd, sd);
				}
			}
		}
		else if (t < 0)  // region 6
		{
			tmp0 = b + e;
			tmp1 = a + d;
			if (tmp1 > tmp0)
			{
				numer = tmp1 - tmp0;
				denom = a - 2* b + c;
				if (numer >= denom)
				{
					branch = 2;
					t = 1;
					s = 0;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 3, x, fd, sd);
				}
				else {
					branch = 7;
					t = numer / denom;
					s = 1 - t;
					sqrDistance = vertex_edge_distance_and_derivs_12(x, 2, 3, fd, sd);
				}
			}
			else {
				t = 0;
				if (tmp1 <= 0)
				{
					branch = 1;
					s = 1;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 2, x, fd, sd);
				}
				else if (d >= 0)
				{
					branch = 3;
					s = 0;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 1, x, fd, sd);
				}
				else {
					branch = 5;
					s = -d / a;
					sqrDistance = vertex_edge_distance_and_derivs_12(x, 1, 2, fd, sd);
				}
			}
		}
		else  // region 1
		{
			numer = c + e - b - d;
			if (numer <= 0)
			{
				branch = 2;
				s = 0;
				t = 1;
				sqrDistance = vertex_vertex_distance_and_derivs_12(0, 3, x, fd, sd);
			}
			else {
				denom = a - 2 * b + c;
				if (numer >= denom)
				{
					branch = 1;
					s = 1;
					t = 0;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 2, x, fd, sd);
				}
				else {
					branch = 6;
					s = numer / denom;
					t = 1 - s;
					sqrDistance = vertex_edge_distance_and_derivs_12(x, 2, 3, fd, sd);
				}
			}
		}
	}

	// testing (remove later)
	if (branch == -1)  std::cout << "-1 branch" << std::endl;
	if (sqrDistance < 0) std::cout << "negative distance" << std::endl;

	output_s = s;
	output_t = t;
	return sqrDistance;
}