#pragma once
#include "Complex.h"
#include <vector>

using namespace std;

void fir1(int n, type_data fc, type_data fs, type_data delay, vector<type_data>& h);
void idft_fir1(vector<Complex>& x_in, int num, vector<type_data>& x_out);
void hamming(int n, vector<type_data>& x);