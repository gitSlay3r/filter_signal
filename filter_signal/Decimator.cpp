#include "Decimator.h"

Decimator::Decimator(int len_conv, type_data fs, int M, int n)
{
    this->n = n;
    this->len_conv = len_conv;
    this->M = M;
    this->fs = fs;
    this->piece_flag = 0;
    this->index_nachalo = 0;
    this->svertka = vector <Complex>(len_conv);
    this->ostatok = vector <Complex>(n);
    this->h = vector <type_data>(n + 1);
    this->x_add_zero = vector <Complex>(len_conv + 2 * n);
    fir1(n, 1.f / M, fs, 0, h);
}

Decimator::~Decimator(){}

void Decimator::filter(vector<Complex>& x_out, int& len_x_out, vector<Complex>& x_in, int& len_x_in, int& flag)
{
    if (flag == 1) {
        memset(&x_add_zero[0], 0, sizeof(Complex) * size(x_add_zero));
    }
    memcpy(&x_add_zero[n], &x_in[0], sizeof(Complex) * len_x_in);
    conv(x_add_zero, len_x_in, h, svertka);

    for (int t = 0; t < n; t++)
        svertka[t] += ostatok[t];

    memcpy(&ostatok[0], &svertka[len_x_in], sizeof(Complex) * n);

    int nachalo = 0;
    int konec = len_x_in;

    if (piece_flag == 0)
        nachalo = n / 2;

    if (flag == 1)
        konec += n / 2;

    int m = 0;
    int index_konec = 0;

    for (int p = index_nachalo + nachalo; p < konec; p += M) {
        x_out[m] = svertka[p];
        index_konec = p - nachalo;
        m++;
    }
    len_x_out = m;

    index_nachalo = M - ((konec - nachalo) - index_konec);
    piece_flag++;
}
