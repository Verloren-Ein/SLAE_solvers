#include "functions.h"

complex<double> conjugate_prod(const vector<complex<double>> &a, const vector<complex<double>> &b)
{
	if(a.size() != b.size())
	{
		cout << "Error: conjugate_prod: vectors of different sizes obtained" << endl;
	}

	double re_result = 0.;
	double im_result = 0.;
	for(size_t i = 0; i < a.size(); i++)
	{
		double re_a = a[i].real();
		double im_a = a[i].imag();
		double re_b = b[i].real();
		double im_b = b[i].imag();

		re_result += re_a * re_b - im_a * im_b;
		im_result += re_a * im_b + im_a * re_b;
	}
	return complex<double>(re_result, im_result);
}

complex<double> dot_prod(const vector<complex<double>> &a, const vector<complex<double>> &b)
{
	if(a.size() != b.size())
	{
		cout << "Error: dot_prod: vectors of different sizes obtained" << endl;
		return false;
	}
	double re_result = 0.;
	double im_result = 0.;
	for(size_t i = 0; i < a.size(); i++)
	{
		double re_a = a[i].real();
		double im_a = a[i].imag();
		double re_b = b[i].real();
		double im_b = b[i].imag();

		re_result += re_a * re_b + im_a * im_b;
		im_result += re_a * im_b - im_a * re_b;
	}
	return complex<double>(re_result, im_result);
}

double dot_prod_2(const vector<complex<double>>& a, const vector<complex<double>>& b)
{
	if(a.size() != b.size())
	{
		cout << "Error: dot_prod: vectors of different sizes obtained" << endl;
		return false;
	}
	double result = 0.0;
	for(size_t i = 0; i < a.size(); i++)
	{
		result += a[i].real()*b[i].real() + a[i].imag()*b[i].imag();
	}
	return result;
}

double dot_prod(const vector<complex<double>> &a)
{
	double result = 0.0;
	for(size_t i = 0; i < a.size(); i++)
	{
		double re = a[i].real();
		double im = a[i].imag();

		result += re * re + im * im;
	}
	return result;
}
