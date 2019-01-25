#pragma once
#include <math.h>
#include <stdio.h>

#define THETA 1.1
int sign(float x) {
	if (x < 0) return -1;
	else return 1;
}

float minmod(float x, float y) {
	return ( sign(x) + sign(y) ) * fmin(fabs(x), fabs(y)) / 2.0;
}

float minmod3(float x, float y, float z) {
   return minmod( x, minmod(y,z) );
}

float approximateDerivative(float x, float y, float z) {
	float l = THETA * (y - x);
	float c = (z - x) / 2.0;
	float r = THETA * (z - y);
	return minmod3(l, c, r);
}
