#include "Interpolator.h"

Interpolator::Interpolator(int length_conv, type_data fs, int M, int n)
{
    this->n = n;
    vector <type_data> h(n + 1);
    fir1(n, 1. / M, fs, 0, h);
    this->h = h;
    this->length_conv = length_conv;
    this->M = M;
    this->fs = fs;
    this->piece_flag = 0;
    this->index_nachalo = 0;

    vector <Complex> for_conv_sig(length_conv - n, 0.0);
    vector <Complex> for_conv_h(length_conv, 0.0);
    vector <Complex> result_conv(length_conv, 0.0);

    this->ostatok = vector <Complex>(n);


    this->for_conv_sig = for_conv_sig;
    this->result_conv = result_conv;
    this->x_add_zero = vector <Complex>(length_conv + n);
}
Interpolator::~Interpolator() {}

void Interpolator::Upsample(vector<Complex>& signal, int koef)
{
    int original_length = signal.size();
    for (int i = 0; i < original_length; i++) {
        for_conv_sig[i * koef] = signal[i];
    }
}

void Interpolator::filter(vector<Complex>& x_out, int& len_x_out, vector<Complex>& x_in, int& len_x_in, int& flag)
{
    if (flag == 1) {
        memset(&x_add_zero[0], 0, sizeof(Complex) * size(x_add_zero));
    }

    Upsample(x_in, M);
    memcpy(&x_add_zero[n], &for_conv_sig[0], sizeof(Complex) * for_conv_sig.size());
    conv(x_add_zero, for_conv_sig.size(), h, result_conv); // свертка временная 

    for (int t = 0; t < n; t++)
        result_conv[t] += ostatok[t];

    memcpy(&ostatok[0], &result_conv[length_conv - n], sizeof(Complex) * n);

    int nachalo = 0;
    int konec = for_conv_sig.size();

    if (piece_flag == 0)
        nachalo = n / 2;

    if (flag == 1)
    {
        konec += n / 2;

    }

    int m = 0;
    int index_konec = 0;

    memcpy(&x_out[0], &result_conv[nachalo], sizeof(Complex) * (konec - nachalo));
    for (int kl = 0; kl < x_out.size(); kl++)
    {
        x_out[kl] = x_out[kl] * M;
    }
    len_x_out = konec - nachalo;

    piece_flag++;
}

