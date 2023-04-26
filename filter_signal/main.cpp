#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>
#include <vector>
#include <memory.h>
#include <ctime>
#include "Complex.h"
using namespace std;

int main()
{
	Complex a(2, 5);
	Complex b(3, 4);
	Complex c = a + b;
	cout << c << endl;

}