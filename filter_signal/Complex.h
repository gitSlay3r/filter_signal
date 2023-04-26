#pragma once
#include <iostream>
#include <cmath>
typedef float type_data;
const type_data PI = (type_data)(atan(1)) * 4;

class Complex {
private:
    type_data re, im;
public:
    Complex(type_data r = 0, type_data i = 0);
    Complex(const Complex& c);
    type_data real() const;
    type_data imag() const;
    ~Complex();
    type_data abs();
    Complex& operator +=(Complex a);
    Complex& operator += (type_data a);
    Complex& operator *= (Complex a);
    Complex& operator *= (type_data a);
    Complex& operator = (const Complex& c);
    Complex operator + (const Complex& c);
    Complex operator - (const Complex& c);
    Complex operator * (const Complex& c);
    Complex conj();
    Complex operator / (const Complex& c);
    friend std::ostream& operator<< (std::ostream&, const Complex&);
    friend std::istream& operator>> (std::istream&, Complex&);
};

Complex operator + (Complex a, type_data b);
Complex operator + (type_data a, Complex b);
Complex operator * (Complex a, type_data b);
Complex operator * (type_data a, Complex b);


