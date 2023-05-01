#pragma once
#include "Complex.h"
#include "fir1.h"
#include <vector>

using namespace std;

void poly_sort(vector<vector<type_data>>& poly_filt, vector<vector<type_data>>& polyphazes, int num, int m, type_data delay, int count_poly);

class Resampler
{
private:
    int piece_flag;
    int vichet;
    int poly_num;
    int i;
    int sdvig_increment;
    int index_vichet;
    int order_poly;
    int index;
    vector <Complex> prev_block;
    vector <Complex> x_add_ostatok;
    vector <int> sdvig;
    vector<vector <type_data>> polyphazes;
    int len_x_add_ostatok;
    int len_piece_signal;
    int L;
    int M;

public:
    Resampler(vector<vector <type_data>>* polyphazes = nullptr, vector <int>* sdvig = nullptr,
        int len_x_add_ostatok = 1, int L = 1, int M = 1, int order_poly = 50); 

    static void create_polyphazes(int& L, int& M, int& order_poly, vector<vector <type_data>>& polyphazes, vector <int>& sdvig); 
    ~Resampler();

    void resample(vector <Complex>& x_out, int& len_x_out, vector <Complex>& x_in, int& len_x_in, int& flag);
};

