#pragma once
#include "common.h"

template <typename T>
class IO
{
public:
	IO() {}
	~IO() {}

	// reads array from txt file
	bool read_txt(const char *file_name, vector<T> &array);
	// reads array from binary file
	bool read_bin(const char *file_name, vector<T> &array);
	bool write_txt(const char *file_name, const vector<T> &array);
	bool write_bin(const char *file_name, const vector<T> &array);

	// gets file name, reads one first number
	bool read_txt(const char *file_name, T *number);
	// gets stream, reads one number
	bool read_txt(ifstream &istr, T *number);

	// gets file name, overwrites only one number
	bool write_txt(const char *file_name, T number);
	// gets stream, appends only one number to the file
	bool write_txt(ofstream &ostr, T number);
};

template<typename T>
inline bool IO<T>::read_txt(const char * file_name, vector<T>& array)
{
	string line;		//line read from file
	string num;			//string containing one number
	stringstream nums;	//string containing all numbers from line
	ifstream in;
	in.open(file_name, ifstream::in);

	if(!in.is_open())
	{
		cout << "Error: Cannot open file " << file_name << " for reading" << endl;
		return false;
	}

	for(size_t i = 0; i < array.size(); )
	{
		std::getline(in, line);
		while(line == "")	//skip empty lines
			std::getline(in, line);
		nums.str(line);
		while(nums >> num)
		{
			T temp;
			string_to_number(num, temp);
			array[i] = temp;
			i++;
		}
		nums.clear();
		nums.str(string());
	}
	in.close();

	return true;
}

template<typename T>
inline bool IO<T>::read_bin(const char * file_name, vector<T>& array)
{
	T temp;
	FILE *fp;

	// open file for reading
	fopen_s(&fp, file_name, "rb");
	if(fp == 0)
	{
		cout << "Error: Cannot open file " << file_name << " for reading" << endl;
		return false;
	}

	for(size_t i = 0; i < array.size(); i++)
	{
		size_t rdb = fread(&temp, sizeof(temp), 1, fp);
		if(rdb == 0) throw;
		array[i] = temp;
	}
	fclose(fp);

	return true;
}

template<typename T>
inline bool IO<T>::write_txt(const char *file_name, const vector<T> &array)
{
	ofstream ostr;
	ostr.open(file_name);
	ostr << fixed;
	ostr.precision(CURRENT_DISPLAY_PRECISION);

	for(size_t i = 0; i < array.size(); i++)
	{
		ostr << array[i] << endl;
	}

	return true;
}

template<typename T>
inline bool IO<T>::read_txt(const char * file_name, T * number)
{
	T temp;

	string num; //string containing one number

	ifstream in;
	in.open(file_name, ifstream::in);

	if(!in.is_open())
	{
		cout << "Error: Cannot open file " << file_name << " for reading" << endl;
		return false;
	}
	std::getline(in, num);

	while(num == "")	//skip empty lines
		std::getline(in, num);
	string_to_number(num, temp);

	in.close();
	*number = temp;

	return true;
}

template<typename T>
inline bool IO<T>::read_txt(ifstream &istr, T * number)
{
	T temp;
	string num; //string containing one number

	std::getline(istr, num);
	while(num == "")	//skip empty lines
		std::getline(istr, num);
	string_to_number(num, temp);

	*number = temp;

	return true;
}

template<typename T>
inline bool IO<T>::write_txt(const char * file_name, T number)
{
	ofstream ostr;
	ostr.open(file_name);
	ostr << fixed;
	ostr.precision(CURRENT_DISPLAY_PRECISION);
	ostr << number << endl;
	ostr.close();

	return true;
}

template<typename T>
inline bool IO<T>::write_txt(ofstream & ostr, T number)
{
	ostr << fixed;
	ostr.precision(CURRENT_DISPLAY_PRECISION);
	ostr << number << endl;

	return true;
}

template<typename T1, typename T2>
bool read_n_eps_maxiter(const char *file_name, T1 *n, T2 *eps, T1 *maxiter)
{
	IO<T1> in_T1;
	IO<T2> in_T2;

	ifstream istr;
	istr.open(file_name, ifstream::in);

	if(!istr.is_open())
	{
		cout << "Error: Cannot open file " << file_name << " for reading" << endl;
		return false;
	}
	in_T1.read_txt(istr, n);
	in_T2.read_txt(istr, eps);
	in_T1.read_txt(istr, maxiter);

	istr.close();

	return true;
}
 