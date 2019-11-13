#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <tchar.h>
#include <time.h>
#include <ctime>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <vector>
#include <string>
#include <chrono>
using namespace std;

/*const for compare with null*/
const double almost_zero = 1e-16;

/*const precision of floating-point numbers compare*/
const double CONST_COMPARE_PRECISION = 1e-12;
/*const precision of floating-point numbers display*/
const streamsize CONST_DISPLAY_PRECISION = cout.precision();
/*current precision of floating-point numbers compare*/
extern double CURRENT_COMPARE_PRECISION;
/*current precision of floating-point numpers display*/
extern streamsize CURRENT_DISPLAY_PRECISION;

/*set new precision of floating-point numbers compare*/
void set_compare_precision(double precision);
/*get current precision of floating-point numbers compare*/
double get_current_compare_precision();
/*set typical precision of floating-point numbers compare*/
void reset_compare_precision();

/*set new precision of floating-point numbers display*/
void set_display_precision(streamsize precision);
/*get current precision of floating-point numbers display*/
streamsize get_current_display_precision();
/*set typical precision of floating-point numbers display*/
void reset_display_precision();

template <typename T>
void string_to_number(string s, T &n)
{
	stringstream str;
	str << s;
	str >> n;
}

template <typename T>
void number_to_string(string &s, T n)
{
	stringstream str;
	str << n;
	str >> s;
}

template <typename T>
void array_shift(vector<T> &array, T shift)
{
	for(size_t i = 0; i < array.size(); i++)
		array[i] += shift;
}
