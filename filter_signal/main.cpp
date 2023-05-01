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
#include "Decimator.h"

using namespace std;
int len_piece = 98304;

int length_block = 8192;
int lim_block = 500000;
int fs_down = 800e3;
int L = 20; // целый интерполятор
int M1 = 5; // целый дециматор
int M2 = 2;
int M3 = 2;
int fs_up = 30e3;
int poly_M = 3; // для полифазника вниз
int poly_L = 4; // для полифазника вверх
type_data  FS_M1 = (fs_down) / M1;
type_data FS_M2 = FS_M1 / M2;
type_data FS_M3 = FS_M2 / M3;
type_data FS_L = fs_up * L;
type_data FS_res = FS_L * poly_M;
type_data fs_after = fs_up * L * poly_L / poly_M;

int shift_f[14] = { 226420, 201808, 175924, 151435, 126585, 76113, 26325, 475, -48627, -73840, -124001, -173430, -199275, -223844 };

FILE* file;
FILE* file2;

void filter_down()
{
    // здесь вниз
    file = fopen("C://Users//maksi//source//repos//filter_signal//filter_signal//tetra_chans_0_000MHz_100KHz.pcm", "rb");
    if (file == NULL) {
        cout << "error" << endl;
    }
    fseek(file, 0, SEEK_END);
    int len_file = ftell(file);
    rewind(file);

    int len_signal = len_file / sizeof(short);
    vector <short> data(len_signal);
    fread(&data[0], 1, len_file, file);

    fclose(file);

    int len_signal_elem = len_signal / 2;
    int len_pieces = 8192;

    int n = 32;
    int order_poly = 32;

    //int count_channel = sizeof(shift_f) / sizeof(int);
    int count_channel = 14;
    int len_conv = len_pieces + n;
    int len_signal_shift = 0;
    int len_sig_decim1 = 0;
    int len_sig_decim2 = 0;

    vector <Complex> signal_shift(len_pieces);
    vector <Complex> signal_after_decimM1((int)(ceil(((type_data)(len_conv + n / 2) / M1))));
    vector <Complex> signal_after_decimM2((int)((type_data)(size(signal_after_decimM1) + n / 2) / M2));

    int len_piece_resample = 0;
    int len_x_add_ostatok = size(signal_after_decimM2) + 3 * order_poly / 2;
    vector <Complex> signal_peredis((int)(ceil((type_data)(len_x_add_ostatok)*L / M3)));

    int dlina1 = 0;
    int dlina2 = 0;
    int len_resample = (int)(ceil((type_data)len_signal_elem * L / M1 / M2 / M3));
    vector <vector<Complex>> matrix_chan(count_channel, vector<Complex>(len_resample));

    vector <Decimator> h1(count_channel);
    for (int i = 0; i < count_channel; i++) {
        h1[i] = Decimator(len_conv, fs_down, M1, n);
    }

    vector <Decimator> h2(count_channel);
    for (int i = 0; i < count_channel; i++) {
        h2[i] = Decimator(size(signal_after_decimM1) + n, FS_M1, M2, n);
    }

    vector<int> sdvig(L);
    vector<vector <type_data>> polyphazes(L, vector<type_data>(order_poly + 1));
    Resampler::create_polyphazes(L, M3, order_poly, polyphazes, sdvig);

    vector <Resampler> h3(count_channel);
    for (int i = 0; i < count_channel; i++) {
        h3[i] = Resampler(&polyphazes, &sdvig, len_x_add_ostatok, L, M3, order_poly);
    }

    int flag = 0;
    int dlina = 0;
    int piece = 0;
    int konec = 0;
    clock_t start = clock();
    while (flag == 0) {

        konec = (piece + 1) * len_pieces;
        if (konec >= len_signal_elem) {
            flag = 1;
            konec = len_signal_elem;
        }

        for (int f_number = 0; f_number < count_channel; f_number++) {

            //перенос по частоте         
            for (int h = len_pieces * piece, m = 0, j = 2 * h; h < konec; h++, m++, j += 2) {
                type_data param = 2 * PI * shift_f[f_number] * (type_data)h / fs_down;
                signal_shift[m] = Complex(data[j], data[j + 1]) * Complex(cos(param), sin(param));
                //signal_shift[m] = Complex(data[j], data[j + 1]);
                len_signal_shift = m + 1;
            }

            //прореживание
            //первый дециматор--------------------------------------------------
            h1[f_number].filter(signal_after_decimM1, len_sig_decim1, signal_shift, len_signal_shift, flag);

            //второй дециматор--------------------------------------------------           
            h2[f_number].filter(signal_after_decimM2, len_sig_decim2, signal_after_decimM1, len_sig_decim1, flag);

            //ресемплер
            h3[f_number].resample(signal_peredis, len_piece_resample, signal_after_decimM2, len_sig_decim2, flag);

            memcpy(&matrix_chan[f_number][dlina], &signal_peredis[0], sizeof(Complex) * len_piece_resample);
        }
        piece++;
        dlina1 += len_sig_decim1;
        dlina2 += len_sig_decim2;
        dlina += len_piece_resample;
    }
    clock_t end = clock();

    if (dlina < len_resample) {
        for (int i = 0; i < count_channel; i++)
            matrix_chan[i].erase(matrix_chan[i].begin() + dlina, matrix_chan[i].end());
    }

    int w = 0;
    vector <Complex> massiv_out(dlina * count_channel);
    for (int i = 0; i < count_channel; i++) {
        memcpy(&massiv_out[i * dlina], &matrix_chan[i][0], sizeof(Complex) * dlina);
    }

    /*file_out = fopen("C:/Users/79814/Desktop/zadanieNIR/chan 30kHz float.pcm", "wb");
    if (file_out == NULL) {
        cout << "error" << endl;
    }
    fwrite(&massiv_out[0], sizeof(type_data), 2 * dlina * count_channel, file_out);
    fclose(file_out);*/

    //Test(massiv_out, count_channel * dlina);

  

    cout << "Decimation + Resample" << endl;
    double elapsed_time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    cout << "Затраченное время: " << elapsed_time << " секунд." << endl;
    int size_out = len_signal_elem * 2;
    /*cout << "Объем данных типа float:" << size_out << " мб." << endl;*/
    type_data speed = (size_out / elapsed_time);
    cout << "Скорость обработки: " << speed / 1e6 << " мSamp/сек" << endl;

}

void filter_up()
{
    // это вверх 

    file2 = fopen("C:/Users/maksi/source/repos/peredis/peredis/chan_30kHz_float.pcm", "rb"); // путь пк
    //file = fopen("C:/Users/AhtemiichukMaxim/source/repos/NIR_peredis/chan_30kHz_float.pcm", "rb"); // путь ноут
    if (file2 == NULL) {
        cout << "error" << endl;
    }
    fseek(file2, 0, SEEK_END);
    int len_file = ftell(file2);
    rewind(file2);
    int len_signal = len_file / sizeof(type_data);
    vector <type_data> data(len_signal);
    fread(&data[0], 1, len_file, file2);
    fclose(file2);

    int len_signal_elem = len_signal / 2;
    vector <Complex> to_up(len_piece, 0.0);
    int n = 32;
    int len_resample = (int)(ceil((type_data)to_up.size() * (L * poly_L) / poly_M));
    vector <Complex> result_block((length_block * L) + n);
    vector <Complex> result(len_resample);
    vector <Complex> block_signal(length_block);
    int order_poly = 32;

    int len_x_add_ostatok = size(result) + 3 * order_poly / 2;

    vector <Complex> signal_peredis((int)(ceil((type_data)(len_x_add_ostatok)*poly_L / poly_M)));

    int numb_section = to_up.size() / length_block;
    int flag = 0;

    int length_conv = (length_block * L) + n;
    int len_result = 0;

    int dlina = 0;
    int dlina2 = 0;
    int len_piece_resample = 0;
    int count_channel = 14; // колличество каналов
    vector<int> sdvig(poly_L);

    vector<vector <type_data>> polyphazes(poly_L, vector<type_data>(order_poly + 1));
    Resampler::create_polyphazes(poly_L, M2, order_poly, polyphazes, sdvig);

    vector <Resampler> f3(count_channel); // объекты класса для полифазника
    for (int i = 0; i < count_channel; i++) {
        f3[i] = Resampler(&polyphazes, &sdvig, len_x_add_ostatok, poly_L, poly_M, order_poly);
    }

    vector <Interpolator> f2(count_channel); // объекты класса для целого 
    for (int i = 0; i < count_channel; i++) {
        f2[i] = Interpolator(length_conv, fs_up, L, n);
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

    cout << "Interpolation + Resample" << endl;
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
}
int main()
{
    setlocale(LC_ALL, "ru");
    
    filter_down();
    filter_up();

}