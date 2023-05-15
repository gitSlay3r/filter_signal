#include "complex.h"
#include <cmath>

Complex::Complex(type_data r, type_data i) : re(r), im(i) {}

Complex::Complex(const Complex& c) : re(c.re), im(c.im) {}

type_data Complex::real() const {
    return re;
}

type_data Complex::imag() const {
    return im;
}

Complex::~Complex() {}

type_data Complex::abs()
{
    return sqrt(re * re + im * im);
}

type_data Complex::arg()
{
    return atan2(re, im);
}

Complex& Complex::operator +=(Complex a)
{
    re += a.re;
    im += a.im;
    return *this;
}

Complex& Complex::operator += (type_data a)
{
    re += a;
    return *this;
}

Complex& Complex::operator *= (Complex a)
{
    re = re * a.re - im * a.im;
    im = re * a.im + im * a.re;
    return *this;
}

Complex& Complex::operator *= (type_data a)
{
    re *= a;
    im *= a;
    return *this;
}

Complex& Complex::operator = (const Complex& c)
{
    if (this != &c) {
        re = c.re;
        im = c.im;
    }
    return *this;
}

Complex Complex::operator + (const Complex& c)
{
    return Complex(re + c.re, im + c.im);
}

Complex Complex::operator - (const Complex& c)
{
    return Complex(re - c.re, im - c.im);
}

Complex Complex::operator * (const Complex& c)
{
    return Complex(re * c.re - im * c.im, re * c.im + im * c.re);
}

Complex Complex::conj()
{
    return Complex(re, -im);
}

Complex Complex::operator / (const Complex& c)
{
    Complex temp;
    type_data r = c.re * c.re + c.im * c.im;
    temp.re = (re * c.re + im * c.im) / r;
    temp.im = (im * c.re - re * c.im) / r;
    return temp;
}

//перегрузка оператора <<
std::ostream& operator<< (std::ostream& out, const Complex& c)
{
    out << "(" << c.re << ", " << c.im << ")";
    return out;
}

// перегрузка оператора >>
std::istream& operator>> (std::istream& in, Complex& c)
{
    in >> c.re >> c.im;
    return in;
}

Complex operator + (Complex a, type_data b)
{
    Complex r = a;
    return r += b;
}

Complex operator + (type_data a, Complex b)
{
    Complex r = b;
    return r += a;
}

Complex operator * (Complex a, type_data b)
{
    Complex r = a;
    return r *= b;
}

Complex operator * (type_data a, Complex b)
{
    Complex r = b;
    return r *= a;
}

inline bool operator==(Complex a, Complex b)
{
    return a.real() == b.real() && a.imag() == b.imag();
}

inline bool operator==(Complex a, type_data b)
{
    return a.real() == b && a.imag() == 0.0;
}

inline bool operator!=(Complex a, type_data b)
{
    return a.real() != b && a.imag() != 0.0;
}

inline bool operator<(Complex a, type_data b)
{
    return (a.real() < b && a.imag() < b) || ((a.real() == b && a.imag() == b) && (a.real() > b && a.imag() > b));
}
