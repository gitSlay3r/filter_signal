#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>
#include <vector>
#include <memory.h>
#include <ctime>
#include "Complex.h"
#include "Interpolator.h"
#include "Resampler.h"
using namespace std;
int len_piece = 98304;

int length_block = 8192;
int lim_block = 500000;
int L = 20; // целый интерполятор
int fs = 30e3;
int M2 = 3; // для полифазника вниз
int L2 = 4; // для полифазника вверх
type_data FS_L = fs * L;
type_data FS_res = FS_L * M2;
type_data fs_after = fs * L * L2 / M2;

int shift_f[14] = { 226420, 201808, 175924, 151435, 126585, 76113, 26325, 475, -48627, -73840, -124001, -173430, -199275, -223844 };

FILE* file;

int main()
{
    setlocale(LC_ALL, "ru");
    double A = 1;
    float T = 5;
    double dt = 1.0 / fs;
    float iter = fs * T;
   
    file = fopen("C:/Users/maksi/source/repos/peredis/peredis/chan_30kHz_float.pcm", "rb"); // путь пк
    //file = fopen("C:/Users/AhtemiichukMaxim/source/repos/NIR_peredis/chan_30kHz_float.pcm", "rb"); // путь ноут
    if (file == NULL) {
        cout << "error" << endl;
    }
    fseek(file, 0, SEEK_END);
    int len_file = ftell(file);
    rewind(file);
    int len_signal = len_file / sizeof(type_data);
    vector <type_data> data(len_signal);
    fread(&data[0], 1, len_file, file);
    fclose(file);

    int len_signal_elem = len_signal / 2;
    vector <Complex> to_up(len_piece, 0.0);
    int n = 32;
    int len_resample = (int)(ceil((type_data)to_up.size() * (L * L2) / M2));
    vector <Complex> result_block((length_block * L) + n);
    vector <Complex> result(len_resample);
    vector <Complex> block_signal(length_block);
    int order_poly = 32;

    int len_x_add_ostatok = size(result) + 3 * order_poly / 2;

    vector <Complex> signal_peredis((int)(ceil((type_data)(len_x_add_ostatok)*L2 / M2)));

    int numb_section = to_up.size() / length_block;
    int flag = 0;

    int length_conv = (length_block * L) + n;
    int len_result = 0;

    int dlina = 0;
    int dlina2 = 0;
    int len_piece_resample = 0;
    int count_channel = 14; // колличество каналов
    vector<int> sdvig(L2);

    vector<vector <type_data>> polyphazes(L2, vector<type_data>(order_poly + 1));
    Resampler::create_polyphazes(L2, M2, order_poly, polyphazes, sdvig);
   
    vector <Resampler> f3(count_channel); // объекты класса для полифазника
    for (int i = 0; i < count_channel; i++) {
        f3[i] = Resampler(&polyphazes, &sdvig, len_x_add_ostatok, L2, M2, order_poly);
    }

    vector <Interpolator> f2(count_channel); // объекты класса для целого 
    for (int i = 0; i < count_channel; i++) {
        f2[i] = Interpolator(length_conv, fs, L, n);
    }

    vector <Complex> result_fs(result.size(), 0.0); // сюда пойдет итог по каналам с переносом

    type_data t_up = 1 / fs_after; // для переноса

    int temp_len = 0;

    clock_t start = clock();

    for (int k = 0; k < count_channel; k++)
    {
        for (int count_up = 0, j = 0; count_up < len_piece; count_up++, j += 2) // запись re, im для 
        {
            to_up[count_up] = Complex(data[(j)+temp_len], data[(j + 1) + temp_len]);
        }

        for (int i = 0; i < numb_section; i++)
        {
            memcpy(&block_signal[0], &to_up[length_block * i], sizeof(Complex) * (length_block));
            if (i == numb_section - 1)
            {
                flag = 1;
            }
            f2[k].filter(result_block, len_result, block_signal, length_block, flag);
            f3[k].resample(signal_peredis, len_piece_resample, result_block, len_result, flag);

            memcpy(&result[(dlina)], &signal_peredis[0], sizeof(Complex) * (len_piece_resample));
            dlina2 += len_result;
            dlina += len_piece_resample;
        }

        type_data param = -2 * PI * shift_f[k];

        // перенос по частоте канала
        for (int l = 0; l < result.size(); l++)
        {
            result_fs[l] += result[l] * Complex(cos(param), sin(param)) * (t_up);
        }

        temp_len += len_piece;
        dlina = 0;
        dlina2 = 0;
    }

    clock_t end = clock();
    double elapsed_time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    cout << "Затраченное время: " << elapsed_time << " секунд." << endl;
    int size_out = result_fs.size() * 2;
    /*cout << "Объем данных типа float:" << size_out << " мб." << endl;*/
    type_data speed = (size_out / elapsed_time);
    cout << "Скорость обработки: " << speed / 1e6 << " мSamp/сек" << endl;

    ofstream toMATLAB;
    //string path = "Y:\\Documents\\MATLAB\\NIR\\cpp_tests\\test_sin.bin"; //путь ноут 
    string path = "C:\\Users\\maksi\\OneDrive\\Документы\\MATLAB\\vs19\\test_sin.bin"; //путь пк
    toMATLAB.open(path, ios::binary); //запись 

    type_data re = 0;
    type_data im = 0;

    for (int i = 0; i < result_fs.size(); i++)
    {
        re = result_fs[i].real();
        toMATLAB.write(reinterpret_cast<const char*>(&re), sizeof(type_data));
        im = result_fs[i].imag();
        toMATLAB.write(reinterpret_cast<const char*>(&im), sizeof(type_data));
    }

    toMATLAB.close();
    return 0;
}