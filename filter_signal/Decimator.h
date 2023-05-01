#pragma once
#include <vector>
#include "Complex.h"
#include "fir1.h"
using namespace std;

class Decimator 
{
private:

    vector <type_data> h;
    vector <Complex> ostatok;
    vector <Complex> svertka;
    vector <Complex> x_add_zero;
    int index_nachalo;
    int n;
    int len_conv;
    int piece_flag;
    type_data fs;
    int M;

public:

    Decimator(int len_conv = 200, type_data fs = 1, int M = 1, int n = 50);
    ~Decimator();

    void filter(vector <Complex>& x_out, int& len_x_out, vector <Complex>& x_in, int& len_x_in, int& flag);
};
