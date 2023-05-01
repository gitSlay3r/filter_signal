#include "Resampler.h"

void poly_sort(vector<vector<type_data>>& poly_filt, vector<vector<type_data>>& polyphazes, int num, int m, type_data delay, int count_poly) {

    for (int i = 0, j = count_poly; i < count_poly; i++, j--) {
        if ((type_data)j / count_poly == delay)
            memcpy(&polyphazes[m][0], &poly_filt[i][0], sizeof(type_data) * num);
    }
}


Resampler::Resampler(vector<vector <type_data>>* polyphazes, vector <int>* sdvig,
    int len_x_add_ostatok, int L, int M, int order_poly) 
{

    if (polyphazes != nullptr && sdvig != nullptr) {
        this->polyphazes = *polyphazes;
        this->sdvig = *sdvig;
    }
    this->order_poly = order_poly;
    this->len_x_add_ostatok = len_x_add_ostatok;
    this->prev_block = vector <Complex>(order_poly);
    this->x_add_ostatok = vector <Complex>(len_x_add_ostatok);
    this->L = L;
    this->M = M;
    this->i = 0;
    this->index = 0;
    this->vichet = 0;
    this->poly_num = 0;
    this->piece_flag = 0;
    this->index_vichet = 0;
    this->sdvig_increment = 0;
    this->len_piece_signal = 0;


}
void Resampler::create_polyphazes(int& L, int& M, int& order_poly, vector<vector<type_data>>& polyphazes, vector<int>& sdvig)
{
    int count_poly = 100;
    vector <type_data> dt(L);
    for (int i = 0; i < L; i++)
        dt[i] = i * (type_data)M / L - i;

    for (int i = 1; i < L; i++) {

        if (trunc(abs(dt[i])) != 0)
            sdvig[i] = (int)(trunc(abs(dt[i])));

        if (dt[i] < 0)
            sdvig[i] = (-1) * sdvig[i];
    }

    vector <type_data> h(order_poly + 1);
    vector <vector<type_data>> poly_filt(count_poly, vector<type_data>(order_poly + 1));
    for (int j = 0, i = count_poly; j < count_poly; j++, i--) {

        type_data delay = (type_data)i / count_poly;
        if (i == count_poly)
            delay = 0;

        memset(&h[0], 0, sizeof(type_data) * (order_poly + 1));
        type_data koef = 0;
        if (M > L)
        {
            type_data koef = L / M;
        }
        else
            type_data koef = M / L;

        fir1(order_poly, koef, 0, delay, h);

        for (int z = 0; z < order_poly + 1; z++)
            poly_filt[j][z] = h[z];

    }

    for (int i = 0; i < L; i++) {
        type_data delay = abs(round((dt[i] - trunc(dt[i])) * 100) / 100);
        if (delay == 0)
            delay = 1;
        poly_sort(poly_filt, polyphazes, order_poly + 1, i, delay, count_poly);
    }

    if (L > M) {
        for (int i = 0; i < L; i++) {
            for (int low = 0, high = order_poly; low < high; low++, high--)
                swap(polyphazes[i][low], polyphazes[i][high]);
        }
    }
}

Resampler::~Resampler(){}

void Resampler::resample(vector<Complex>& x_out, int& len_x_out, vector<Complex>& x_in, int& len_x_in, int& flag)
{
    if (flag == 1) {
        memset(&x_add_ostatok[0], 0, sizeof(Complex) * size(x_add_ostatok));
    }
    if (piece_flag == 0) {

        memcpy(&x_add_ostatok[order_poly / 2], &x_in[0], sizeof(Complex) * len_x_in);

        len_piece_signal = len_x_in + order_poly / 2;
        index_vichet = len_piece_signal - order_poly;
    }

    if (piece_flag != 0) {

        memcpy(&x_add_ostatok[0], &prev_block[0], sizeof(Complex) * order_poly);
        memcpy(&x_add_ostatok[order_poly], &x_in[0], sizeof(Complex) * len_x_in);

        len_piece_signal = len_x_in + order_poly;
        index_vichet = len_x_in;
        if (flag == 1) {
            len_piece_signal += order_poly / 2;
        }
    }

    len_x_out = 0;
    memset(&x_out[0], 0, sizeof(Complex) * size(x_out));
    index = i + sdvig_increment + sdvig[poly_num] - vichet;
    while (index + order_poly < len_piece_signal) {

        for (int l = 0; l <= order_poly; l++) {
            x_out[len_x_out] += x_add_ostatok[index + l] * polyphazes[poly_num][l];
        }

        poly_num++; i++; len_x_out++;

        if (poly_num > L - 1) {
            sdvig_increment += M - L;
            poly_num = 0;
        }

        index = i + sdvig_increment + sdvig[poly_num] - vichet;
    }

    vichet += index_vichet;

    if (flag != 1) {
        memcpy(&prev_block[0], &x_in[len_x_in - order_poly], sizeof(Complex) * order_poly);
    }

    piece_flag++;
}
