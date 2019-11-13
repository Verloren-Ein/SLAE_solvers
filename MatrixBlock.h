#pragma once
#include "common.h"
#include "IO.h"

class MatrixBlock
{
protected:
	bool size_set;

	//convert double matrix representation to complex
	void double_to_complex(const vector<size_t> &index_array, const size_t end,
		const vector<double> &dmatrix_values, vector<complex<double>> &cmatrix_values);

public:
	size_t n;
	size_t nb;

	vector<size_t> ig;
	vector<size_t> jg;
	vector<size_t> idi;
	vector<size_t> ijg;

	vector<double> di;
	vector<double> rp;

	vector<complex<double>> complex_di;
	vector<complex<double>> complex_rp;

	virtual void clear();
	virtual void values_to_complex() = 0;

	// reads right part from file
	bool read_rp_txt(const char *f_rp);
	bool read_rp_bin(const char *f_rp);

	MatrixBlock()
		:size_set(false)
	{}
	MatrixBlock(size_t n)
		:n(n)
		, nb(n / 2)
		, size_set(true)
	{}
	virtual ~MatrixBlock() = 0
	{}
};

class MatrixBlock_Symmetric : public MatrixBlock
{
public:
	vector<double> gg;
	vector<complex<double>> complex_gg;

	// reads matrix from files
	bool read_matrix_txt(const char *f_ig, const char *f_jg, const char *f_idi,
		const char *f_ijg, const char *f_di, const char *f_gg);
	bool read_matrix_bin(const char *f_ig, const char *f_jg, const char *f_idi,
		const char *f_ijg, const char *f_di, const char *f_gg);

	// reads slae from files
	bool read_slae_txt(const char *f_ig, const char *f_jg, const char *f_idi, 
		const char *f_ijg, const char *f_di, const char *f_gg, const char *f_rp);
	bool read_slae_bin(const char *f_ig, const char *f_jg, const char *f_idi,
		const char *f_ijg, const char *f_di, const char *f_gg, const char *f_rp);

	void clear() override;
	void values_to_complex() override;

	MatrixBlock_Symmetric() : MatrixBlock(){}
	MatrixBlock_Symmetric(size_t n) : MatrixBlock(n){}
	~MatrixBlock_Symmetric() {}
};

class MatrixBlock_Asymmetric : public MatrixBlock
{
public:
	vector<double> ggu;
	vector<double> ggl;
	vector<complex<double>> complex_ggu;
	vector<complex<double>> complex_ggl;

	// reads matrix from files
	bool read_matrix_txt(const char *f_ig, const char *f_jg, const char *f_idi,
		const char *f_ijg, const char *f_di, const char *f_ggu, const char *f_ggl);
	bool read_matrix_bin(const char *f_ig, const char *f_jg, const char *f_idi,
		const char *f_ijg, const char *f_di, const char *f_ggu, const char *f_ggl);

	// reads slae from files
	bool read_slae_txt(const char *f_ig, const char *f_jg, const char *f_idi,
		const char *f_ijg, const char *f_di, const char *f_ggu, const char *f_ggl, const char *f_rp);
	bool read_slae_bin(const char *f_ig, const char *f_jg, const char *f_idi,
		const char *f_ijg, const char *f_di, const char *f_ggu, const char *f_ggl, const char *f_rp);

	void clear() override;
	void values_to_complex() override;

	MatrixBlock_Asymmetric() : MatrixBlock() {}
	MatrixBlock_Asymmetric(size_t n) : MatrixBlock(n) {}
	~MatrixBlock_Asymmetric(){}
};


