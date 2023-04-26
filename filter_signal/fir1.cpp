#include "fir1.h"

void idft_fir1(vector<Complex>& x_in, int num, vector<type_data>& x_out) {

    for (int n = 0; n < num; n++) {
        for (int k = 0; k < num; k++) {
            type_data phi = (type_data)((2 * PI * k * n) / num);
            Complex t = x_in[k] * Complex(cos(phi), sin(phi));
            x_out[n] += t.real();
        }
        x_out[n] = x_out[n] / num;
    }
}
void fir1(int n, type_data fc, type_data fs, type_data delay, vector<type_data>& h) {

    int N = n + 1;
    type_data K = (type_data)n / 2 * fc;
    vector <int> H(N);
    vector <type_data> f1(N);

    for (int i = -N / 2, k = 0; i <= N / 2 - 1; i++, k++)
        f1[k] = (type_data)i / N;

    vector <type_data> phi(N);
    if (N % 2 != 0) {
        for (int i = n / 2 - (int)K; i < n / 2 + 1 + (int)K; i++)
            H[i] = 1; // À×Õ ÔÍ× ôèëüòðà

        for (int i = 0; i < N; i++)
            phi[i] = (-(type_data)n / 2 - delay) * 2 * PI * f1[i]; // Ô×Õ ÔÍ× ôèëüòðà

        for (int i = 0; i < N; i++)
            phi[i] = phi[i] - phi[n / 2];
    }
    else {
        for (int i = N / 2 - (int)K - 1; i < N / 2 + (int)K; i++)
            H[i] = 1; // À×Õ ÔÍ× ôèëüòðà

        for (int i = 0; i < N; i++)
            phi[i] = (-(type_data)N / 2 + 0.5f - delay) * 2 * PI * f1[i]; // Ô×Õ ÔÍ× ôèëüòðà

        for (int i = 0; i < N; i++)
            phi[i] = phi[i] - phi[N / 2];
    }

    int m = (N - 1) / 2 + 1;
    vector <Complex> h_shift(N);
    for (int i = 0; i < N; i++) {
        h_shift[m] = (type_data)H[i] * Complex(cos(phi[i]), sin(phi[i]));;
        m++;
        if (m == N)
            m = 0;
    }

    idft_fir1(h_shift, N, h);
    hamming(n, h);
}
void hamming(int n, vector<type_data>& x)
{
    type_data a0 = 0.54f;
    for (int i = 0; i < n + 1; i++) {
        type_data phi = (2 * PI * (type_data)i / ((type_data)n + 1));
        x[i] *= a0 - (1 - a0) * cos(phi);
    }
}