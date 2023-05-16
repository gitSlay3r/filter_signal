#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>
#include <vector>
#include <memory.h>
#include <ctime>
#include <stdio.h>
#include <ipp.h>
#include "Complex.h"
#include "Interpolator.h"
#include "Resampler.h"
#include "Decimator.h"

using namespace std;


int lim_block = 500000;
int fs_down = 800e3;

//параметры вниз
vector <int> M{ 2 , 2 , 2}; // коэфицент для стадии
int Sum_M = 8;
int poly_M_down = 4; // для полифазника вниз
int poly_L_down = 3; // для полифазника вверх

//параметры вверх 
int L = 8; // целый интерполятор
int poly_L_up = 4;
int poly_M_up = 3;

//type_data FS_M1 = (fs_down) / M1;
//type_data FS_M2 = FS_M1 / M2;
//type_data FS_L = fs_up * L;
//type_data FS_res = FS_L * poly_M;

int shift_f[14] = { 226420, 201808, 175924, 151435, 126585, 76113, 26325, 475, -48627, -73840, -124001, -173430, -199275, -223844 };

FILE* file;
FILE* file2;

// частота дискретизации для каждой стадии 
void fs_stage_M(int& fs, vector <int>& fs_stage, vector <int> M)
{
    int temp = fs;
    for (int i = 0; i < M.size(); i++)
    {
        fs_stage[i] = temp / M[i];
        temp = fs_stage[i];
    }
}

void Stage_filter(vector <vector<Decimator>>& decim, int& f_number, vector <Complex>& x_in,int& len_signal_buff,vector<int>& M, int& fs, int& len_pieces, int& n, vector <Complex>& signal_buff_M, int& flag)
{
    int len_conv = len_pieces + n; // длина свертки 
    int fs_stage = fs; // частота на каждом этапе  
    int original_size = len_conv;
    //входной буфер
    vector <Complex> signal_buff = x_in;
    int temp = 0;
    for (int i = 0; i < M.size(); i++)
    {
        int len_signal_buff_M = 0;
        // выходной буффер
        signal_buff_M.resize((int)(ceil(((type_data)(original_size + n/2 - temp)/ M[i]))));
        // прореживание
        decim[i][f_number].filter(signal_buff_M, len_signal_buff_M, signal_buff, len_signal_buff, flag);
        // для следующей итерации(другой коэффициент)
        signal_buff = signal_buff_M;
        fs_stage = fs_stage / M[i];
        original_size += n / 2;
        original_size = original_size / M[i] + n;
        len_signal_buff = len_signal_buff_M;
        temp += n;
    }   
}
vector<vector<Decimator>> Make_class_dec(vector<int>& M, int& count_channel, int& len_pieces, int& n,int& fs)
{
    // первое значание- класс для опредленного коэффициента 
    // второе значение- класс для каннала
    int len_conv = len_pieces + n; // длина свертки 
    int fs_stage = fs; // частота на каждом этапе  
    int original_size = len_conv;
    vector <vector<Decimator>> decim(M.size(), vector<Decimator> (count_channel));
    for (int k = 0; k < M.size(); k++)
    {
        for (int i = 0; i < count_channel; i++)
        {
            decim[k][i] = Decimator(original_size, fs_stage, M[k], n);
        }
        fs_stage = fs_stage / M[k];
        original_size += n / 2;
        original_size = original_size / M[k] + n;
    }
    return decim;
}

void filter_down()
{
    // здесь вниз
    //file = fopen("C://Users//maksi//source//repos//filter_signal//filter_signal//tetra_chans_0_000MHz_100KHz.pcm", "rb"); // путь пк
    file = fopen("C://Users//AhtemiichukMaxim//source//repos//filter_signal//filter_signal//tetra_chans_0_000MHz_100KHz.pcm", "rb"); // путь ноут
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
    
    int len_piece_resample = 0;
   

    int dlina1 = 0;
    int dlina2 = 0;
    int len_resample = (int)(ceil((type_data)len_signal_elem * poly_L_down / poly_M_down));

    for (int i = 0, temp = len_resample; i < M.size(); i++)
    {
        len_resample = (int)(ceil((type_data)temp / M[i]));
        temp = len_resample;
    }

    vector <vector<Complex>> matrix_chan(count_channel, vector<Complex>(len_resample));

 /*   vector <Decimator> h1(count_channel);
    for (int i = 0; i < count_channel; i++) {
        h1[i] = Decimator(len_conv, fs_down, M1, n);
    }

    vector <Decimator> h2(count_channel);
    for (int i = 0; i < count_channel; i++) {
        h2[i] = Decimator(size(signal_after_decimM1) + n, FS_M1, M2, n);
    }*/

    vector <int> fs_M(M.size());
    fs_stage_M(fs_down, fs_M, M);

    vector<int> sdvig(poly_L_down);
    vector<vector <type_data>> polyphazes(poly_L_down, vector<type_data>(order_poly + 1));
    type_data fs_resampler = fs_M[M.size()]; // частоты дискретизации перед ресемплером 
    Resampler::create_polyphazes(poly_L_down, poly_M_down, order_poly, polyphazes, sdvig, fs_resampler);

    vector <Resampler> h3(count_channel);
    /*for (int i = 0; i < count_channel; i++) {
        h3[i] = Resampler(&polyphazes, &sdvig, len_x_add_ostatok, poly_L, poly_M, order_poly);
    }*/

    int flag = 0;
    int dlina = 0;
    int piece = 0;
    int konec = 0;
    clock_t start = clock();

    
  
    // вычисление размера вектора для ресемплера
    int len_signal_after_decimM = 0;
    for (int i = 0, temp = n/2, temp2 = len_conv; i < M.size(); i++, temp = n/2 )
    {
        len_signal_after_decimM = ceil(temp2 + temp)/ M[i];
        temp2 = len_signal_after_decimM;
    }
    vector <Complex> signal_after_decimM(len_signal_after_decimM);
    // создание объектов класса дециматор для каждого канала 
    vector <vector<Decimator>> decim = Make_class_dec(M, count_channel, len_pieces, n, fs_down);
    //размер и сам буффер для ресемплера
    int res_ostatok = size(signal_after_decimM) + 3 * order_poly / 2;
    vector <Complex> signal_peredis((int)(ceil((type_data)(res_ostatok)*poly_L_down / poly_M_down)));
    // создание объектов класса ресемплер для каждого канала
    for (int i = 0; i < count_channel; i++)
    {
        h3[i] = Resampler(&polyphazes, &sdvig, res_ostatok, poly_L_down, poly_M_down, order_poly);
    }
    while (flag == 0) 
    {

        konec = (piece + 1) * len_pieces;
        if (konec >= len_signal_elem) 
        {
            flag = 1;
            konec = len_signal_elem;
        }

        for (int f_number = 0; f_number < count_channel; f_number++) 
        {

            //перенос по частоте         
            for (int h = len_pieces * piece, m = 0, j = 2 * h; h < konec; h++, m++, j += 2) {
                type_data param = 2 * PI * shift_f[f_number] * (type_data)h / fs_down;
                signal_shift[m] = Complex(data[j], data[j + 1]) * Complex(cos(param), sin(param));
                //signal_shift[m] = Complex(data[j], data[j + 1]);
                len_signal_shift = m + 1;
            }

            //прореживание
            Stage_filter(decim, f_number, signal_shift, len_signal_shift, M, fs_down, len_pieces, n, signal_after_decimM, flag);
            
            //ресемплер
            h3[f_number].resample(signal_peredis, len_piece_resample, signal_after_decimM, len_signal_shift, flag);
            memcpy(&matrix_chan[f_number][dlina], &signal_peredis[0], sizeof(Complex) * len_piece_resample);
        }
        piece++;
        dlina1 += len_signal_shift;
        dlina += len_piece_resample;
    }
    clock_t end = clock();

    if (dlina < len_resample) 
    {
        for (int i = 0; i < count_channel; i++)
            matrix_chan[i].erase(matrix_chan[i].begin() + dlina, matrix_chan[i].end());
    }

    int w = 0;
    vector <Complex> massiv_out(dlina * count_channel);
    for (int i = 0; i < count_channel; i++) {
        memcpy(&massiv_out[i * dlina], &matrix_chan[i][0], sizeof(Complex) * dlina);
    }

    FILE* file_out;
    //file_out = fopen("C://Users//maksi//source//repos//filter_signal//filter_signal//chan_30kHz_float.pcm", "wb"); // путь пк
    file_out = fopen("Y:\\Documents\\MATLAB\\NIR\\cpp_tests\\chan_30kHz_float_stage.pcm", "wb"); // путь ноут
    if (file_out == NULL) 
    {
        cout << "error" << endl;
    }
    fwrite(&massiv_out[0], sizeof(type_data), 2 * dlina * count_channel, file_out);
    fclose(file_out);

    //Test(massiv_out, count_channel * dlina);


    cout << "Decimation + Resample" << endl;
    double elapsed_time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    cout << "Затраченное время: " << elapsed_time << " секунд." << endl;
    int size_out = len_signal_elem * 2;
    /*cout << "Объем данных типа float:" << size_out << " мб." << endl;*/
    type_data speed = (size_out / elapsed_time);
    cout << "Скорость обработки: " << speed / 1e6 << " мSamp/сек" << endl;

}

int length_block = 8192;


// Функция для выполнения переноса по частоте сигнала после интерполяции
void frequencyShift(vector<Complex>& signal, double shiftFrequency, double sampleRate) {
    int N = signal.size();
    double dt = 1.0 / sampleRate;

    // Выполнение преобразования Фурье
    vector<Complex> spectrum(N);
    for (int k = 0; k < N; ++k) {
        for (int n = 0; n < N; ++n) {
            double theta = 2.0 * PI * k * n / N;
            spectrum[k] += signal[n] * Complex(cos(-theta), sin(-theta));
        }
    }

    // Выполнение переноса частоты путем умножения на комплексный экспонент
    for (int k = 0; k < N; ++k) {
        double freq = k / (N * dt);
        double phaseShift = 2.0 * PI * shiftFrequency * freq;
        spectrum[k] *= Complex(cos(phaseShift), sin(phaseShift));
    }

    // Выполнение обратного преобразования Фурье
    for (int n = 0; n < N; ++n) {
        signal[n] = Complex(0.0, 0.0);
        for (int k = 0; k < N; ++k) {
            double theta = 2.0 * PI * k * n / N;
            signal[n] += spectrum[k] * Complex(cos(theta), sin(theta));
        }
        signal[n] = signal[n]/ N;
    }
}

// частота дискретизации для каждой стадии 
// 
void fs_stage_L(int& fs, vector <int>& fs_stage, vector <int>& L)
{
    int temp = fs;
    for (int i = 0; i < L.size(); i++)
    {
        fs_stage[i] = temp * L[i];
        temp = fs_stage[i];
    }
}

void filter_up()
{
    int fs_up = fs_down * poly_L_down / Sum_M * poly_M_down;
   
    //file2 = fopen("C://Users//maksi//source//repos//peredis//peredis//chan_30kHz_float.pcm", "rb"); // путь пк
    //file2 = fopen("C://Users//maksi//source//repos//filter_signal//filter_signal//chan_30kHz_float.pcm", "rb");
    file2 = fopen("Y:\\Documents\\MATLAB\\NIR\\cpp_tests\\chan_30kHz_float_stage.pcm", "rb"); // путь ноут
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

    type_data fs_after = fs_up * L * poly_L_up / poly_M_up;

    int count_channel = 14; // колличество каналов

    int len_piece = (data.size() / (2 * count_channel)); // длина входного сигнала
    int len_signal_elem = len_signal / 2;
    vector <Complex> to_up(len_piece, 0.0);
    int n = 32;
    int len_resample = (int)(ceil((type_data)to_up.size() * (L * poly_L_up) / poly_M_up));
    vector <Complex> result_block((length_block * L) + n);
    vector <Complex> result(len_resample);
    vector <Complex> block_signal(length_block);
    int order_poly = 32;

    int len_x_add_ostatok = size(result) + 3 * order_poly / 2;

    vector <Complex> signal_peredis((int)(ceil((type_data)(len_x_add_ostatok)*poly_L_up / poly_M_up)));

    int numb_section = to_up.size() / length_block;
    int flag = 0;

    int length_conv = (length_block * L) + n;
    int len_result = 0;

    int dlina = 0;
    int dlina2 = 0;
    int len_piece_resample = 0;
  
   
    vector<int> sdvig(poly_L_up);

    vector<vector <type_data>> polyphazes(poly_L_up, vector<type_data>(order_poly + 1));
    type_data fs_to_res = fs_up * L;
    Resampler::create_polyphazes(poly_L_up, poly_M_up, order_poly, polyphazes, sdvig, fs_to_res);
 
    vector <Resampler> f3(count_channel); // объекты класса для полифазника
    for (int i = 0; i < count_channel; i++) {
        f3[i] = Resampler(&polyphazes, &sdvig, len_x_add_ostatok, poly_L_up, poly_M_up, order_poly);
    }
    
    vector <Interpolator> f2(count_channel); // объекты класса для целого 
    for (int i = 0; i < count_channel; i++) {
        f2[i] = Interpolator(length_conv, fs_up, L, n);
    }

    vector <Complex> result_fs(result.size(), 0.0); // сюда пойдет итог по каналам с переносом

    type_data t_up = 1 / fs_after; // для переноса

    int temp_len = 0;

    clock_t start = clock();
    int kl = 0;
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

       
        type_data param = -2 * PI * shift_f[k] * t_up;
        //frequencyShift(result, shift_f[k], fs_after);
        // перенос по частоте канала
         // Выполнение переноса частоты путем умножения на комплексный экспонент

        int size = result.size();
        for (int n = 0; n < size; ++n) {
            type_data freq = n / (size * fs_down);
            type_data phaseShift = -2.0 * PI * shift_f[k] * freq;
            result_fs[n] += result[n] * Complex(cos(phaseShift), sin(phaseShift));
        }

        /*for (int l = 0; l < result.size(); l++)
        {
            result_fs[l] += result[l];
        }*/

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
    //string path = "C:\\Users\\maksi\\OneDrive\\Документы\\MATLAB\\vs19\\out_filter.bin"; //путь пк
    //toMATLAB.open(path, ios::binary); //запись 

    type_data re = 0;
    type_data im = 0;
    FILE *file_out2;
    //file_out2 = fopen("C:\\Users\\maksi\\OneDrive\\Документы\\MATLAB\\vs19\\out_filter.bin", "wb"); // путь пк
    file_out2 = fopen("Y:\\Documents\\MATLAB\\NIR\\cpp_tests\\test_sin.bin", "wb"); // путь ноут 
   /* for (int i = 0; i < result_fs.size(); i++)
    {
        re = result_fs[i].real();
        toMATLAB.write(reinterpret_cast<const char*>(&re), sizeof(type_data));
        im = result_fs[i].imag();
        toMATLAB.write(reinterpret_cast<const char*>(&im), sizeof(type_data));
    }*/
    fwrite(&result_fs[0], sizeof(type_data),2 * result_fs.size(), file_out2);
    fclose(file_out2);
    //toMATLAB.close();
}

int main()
{
    setlocale(LC_ALL, "ru");
    filter_down();
    filter_up();
   //// Входные сигналы
   // Ipp32f signal1[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
   // Ipp32f signal2[8] = { 1, 2, 1, 2, 1, 2, 1, 2 };
   // //Ipp32f vector[8] = {1, 2, 3, 4, 5, 6, 7, 8 };
   // //Ipp32f vector[8] = {1, 2, 1, 2, 1, 2, 1, 2 };
   // // Размер сигналов
   // int signalSize = 8;
 
   // // Выделение памяти для результатов
   // Ipp32f* fftResult1 = ippsMalloc_32f(signalSize);
   // Ipp32f* fftResult2 = ippsMalloc_32f(signalSize);
   // Ipp32f* convResult = ippsMalloc_32f(2 * signalSize - 1);

   // // Создание объектов FFT
   // IppsFFTSpec_R_32f* fftSpec = nullptr;
   // int specSize, initSize, bufferSize;
   // Ipp8u* buffer;

   // ippsFFTGetSize_R_32f(signalSize, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, &specSize, &initSize, &bufferSize);
   // fftSpec = (IppsFFTSpec_R_32f*)ippsMalloc_8u(specSize);
   // buffer = ippsMalloc_8u(bufferSize);
   // // Инициализация объекта FFT
   // ippsFFTGetSize_R_32f(signalSize, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, &specSize, &initSize, &bufferSize);
   //
   // // Прямое преобразование Фурье для обоих сигналов
   // ippsFFTFwd_RToCCS_32f(signal1, fftResult1, fftSpec, buffer);
   // ippsFFTFwd_RToCCS_32f(signal2, fftResult2, fftSpec, buffer);
   // 

   // // Поэлементное умножение в частотной области
   // ippsMul_32f(fftResult1, fftResult2, convResult, signalSize);

   // // Обратное преобразование Фурье для получения свертки
   // ippsFFTInv_CCSToR_32f(convResult, convResult, fftSpec, buffer);

   // FILE* file_out_test;
   // //file_out = fopen("C://Users//maksi//source//repos//filter_signal//filter_signal//chan_30kHz_float.pcm", "wb"); // путь пк
   // file_out_test = fopen("Y:\\Documents\\MATLAB\\NIR\\cpp_tests\\test_ipp.pcm", "wb"); // путь ноут
   // if (file_out_test == NULL)
   // {
   //     cout << "error" << endl;
   // }
   // fwrite(&convResult[0], sizeof(type_data), 2 * signalSize - 1, file_out_test);
   // fclose(file_out_test);

   // // Освобождение ресурсов
   // ippsFree(fftResult1);
   // ippsFree(fftResult2);
   // ippsFree(convResult);
   // ippsFree(buffer);
   // ippsFree(fftSpec);

    return 0;
}