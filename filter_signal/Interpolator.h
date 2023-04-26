#pragma once
#include "Complex.h"
#include <vector>
#include <iostream>

using namespace std;

void conv(vector<Complex>& x, int num1, vector<type_data>& h, vector<Complex>& svertka);
void fir1(int n, type_data fc, type_data fs, type_data delay, vector<type_data>& h);
void idft_fir1(vector<Complex>& x_in, int num, vector<type_data>& x_out);
void hamming(int n, vector<type_data>& x);
class Interpolator
{
private:
    vector <type_data> h;
    vector <Complex> ostatok;
    vector <Complex> x_add_zero;
    vector <Complex> for_conv_sig;
    int index_nachalo;
    int n;
    int length_conv;
    int piece_flag;
    type_data fs;
    int M;
    vector <Complex> result_conv;

public:

    Interpolator(int length_conv = 200, type_data fs = 1, int M = 1, int n = 50);
    ~Interpolator();

    void Upsample(vector<Complex>& signal, int koef);
    void filter(vector <Complex>& x_out, int& len_x_out, vector <Complex>& x_in, int& len_x_in, int& flag);
};

