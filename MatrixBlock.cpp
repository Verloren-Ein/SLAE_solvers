#include "MatrixBlock.h"

void MatrixBlock::double_to_complex(const vector<size_t> &index_array, const size_t end, const vector<double> &dmatrix_values, vector<complex<double>>& cmatrix_values)
{
	cmatrix_values.clear();
	for(size_t i = 0, count = 0; count < end; )
	{
		cmatrix_values.push_back(complex<double>(dmatrix_values[count], 0));
		count++;
		if(index_array[i + 1] - index_array[i] == 2)
		{
			cmatrix_values[i].imag(dmatrix_values[count]);
			count++;
		}
		i++;
	}
}

void MatrixBlock::clear()
{
	ig.clear();
	jg.clear();
	idi.clear();
	ijg.clear();
	di.clear();
	rp.clear();
	complex_di.clear();
	complex_rp.clear();
}

bool MatrixBlock::read_rp_txt(const char * f_rp)
{
	if(!size_set)
	{
		cout << "Error: Size of block matrix is not set " << endl;
		return false;
	}

	IO<double> io_double;
	rp.resize(n);
	io_double.read_txt(f_rp, rp);

	return true;
}

bool MatrixBlock::read_rp_bin(const char * f_rp)
{
	if(!size_set)
	{
		cout << "Error: Size of block matrix is not set " << endl;
		return false;
	}

	IO<double> io_double;
	rp.resize(n);
	io_double.read_bin(f_rp, rp);

	return true;
}


bool MatrixBlock_Symmetric::read_matrix_txt(const char * f_ig, const char * f_jg, const char * f_idi, const char * f_ijg, const char * f_di, const char * f_gg)
{
	if(!size_set)
	{
		cout << "Error: Size of block matrix is not set " << endl;
		return false;
	}

	/*CLASSES FOR INPUT/OUTPUT*/
	IO<double> io_double;
	IO<size_t> io_size_t;

	/*reading idi*/
	idi.resize(nb + 1);
	io_size_t.read_txt(f_idi, idi);
	if(idi[0] != 0)
		array_shift(idi, (size_t)-1);

	/*reading di*/
	di.resize(idi[nb]);
	io_double.read_txt(f_di, di);

	/*reading ig*/
	ig.resize(nb + 1);
	io_size_t.read_txt(f_ig, ig);
	if(ig[0] != 0)
		array_shift(ig, (size_t)-1);

	/*reading jg*/
	jg.resize(ig[nb]);
	io_size_t.read_txt(f_jg, jg);
	if(jg[0] != 0)
		array_shift(jg, (size_t)-1);

	/*reading ijg*/
	ijg.resize(ig[nb] + 1);
	io_size_t.read_txt(f_ijg, ijg);
	if(ijg[0] != 0)
		array_shift(ijg, (size_t)-1);

	/*reading gg*/
	gg.resize(ijg[ig[nb]]);
	io_double.read_txt(f_gg, gg);

	return true;
}

bool MatrixBlock_Symmetric::read_matrix_bin(const char * f_ig, const char * f_jg, const char * f_idi, const char * f_ijg, const char * f_di, const char * f_gg)
{
	if(!size_set)
	{
		cout << "Error: Size of block matrix is not set " << endl;
		return false;
	}

	/*CLASSES FOR INPUT/OUTPUT*/
	IO<double> io_double;
	IO<size_t> io_size_t;

	/*reading idi*/
	idi.resize(nb + 1);
	io_size_t.read_bin(f_idi, idi);
	if(idi[0] != 0)
		array_shift(idi, (size_t)-1);

	/*reading di*/
	di.resize(idi[nb]);
	io_double.read_bin(f_di, di);

	/*reading ig*/
	ig.resize(nb + 1);
	io_size_t.read_bin(f_ig, ig);
	if(ig[0] != 0)
		array_shift(ig, (size_t)-1);

	/*reading jg*/
	jg.resize(ig[nb]);
	io_size_t.read_bin(f_jg, jg);
	if(jg[0] != 0)
		array_shift(jg, (size_t)-1);

	/*reading ijg*/
	ijg.resize(ig[nb] + 1);
	io_size_t.read_bin(f_ijg, ijg);
	if(ijg[0] != 0)
		array_shift(ijg, (size_t)-1);

	/*reading gg*/
	vector<double> _gg;	//aux vector
	gg.resize(ijg[ig[nb]]);
	io_double.read_bin(f_gg, gg);

	return true;
}

bool MatrixBlock_Symmetric::read_slae_txt(const char * f_ig, const char * f_jg, const char * f_idi, const char * f_ijg, const char * f_di, const char * f_gg, const char * f_rp)
{
	read_matrix_txt(f_ig, f_jg, f_idi, f_ijg, f_di, f_gg);
	read_rp_txt(f_rp);

	return true;
}

bool MatrixBlock_Symmetric::read_slae_bin(const char * f_ig, const char * f_jg, const char * f_idi, const char * f_ijg, const char * f_di, const char * f_gg, const char * f_rp)
{
	read_matrix_bin(f_ig, f_jg, f_idi, f_ijg, f_di, f_gg);
	read_rp_bin(f_rp);

	return true;
}

void MatrixBlock_Symmetric::values_to_complex()
{
	complex_di.clear();
	complex_gg.clear();
	complex_rp.clear();

	double_to_complex(idi, idi[nb], di, complex_di);
	double_to_complex(ijg, ijg[ig[nb]], gg, complex_gg);
	for(size_t i = 0; i < nb; i++)
	{
		complex_rp.push_back(complex<double>(rp[2 * i], rp[2 * i + 1]));
	}
}

void MatrixBlock_Symmetric::clear()
{
	gg.clear();
	complex_gg.clear();
}


bool MatrixBlock_Asymmetric::read_matrix_txt(const char * f_ig, const char * f_jg, const char * f_idi, const char * f_ijg, const char * f_di, const char * f_ggu, const char * f_ggl)
{
	if(!size_set)
	{
		cout << "Error: Size of block matrix is not set " << endl;
		return false;
	}

	/*CLASSES FOR INPUT/OUTPUT*/
	IO<double> io_double;
	IO<size_t> io_size_t;

	/*reading idi*/
	idi.resize(nb + 1);
	io_size_t.read_txt(f_idi, idi);
	if(idi[0] != 0)
		array_shift(idi, (size_t)-1);

	/*reading di*/
	di.resize(idi[nb]);
	io_double.read_txt(f_di, di);

	/*reading ig*/
	ig.resize(nb + 1);
	io_size_t.read_txt(f_ig, ig);
	if(ig[0] != 0)
		array_shift(ig, (size_t)-1);

	/*reading jg*/
	jg.resize(ig[nb]);
	io_size_t.read_txt(f_jg, jg);
	if(jg[0] != 0)
		array_shift(jg, (size_t)-1);

	/*reading ijg*/
	ijg.resize(ig[nb] + 1);
	io_size_t.read_txt(f_ijg, ijg);
	if(ijg[0] != 0)
		array_shift(ijg, (size_t)-1);

	/*reading ggu*/
	ggu.resize(ijg[ig[nb]]);
	io_double.read_txt(f_ggu, ggu);

	/*reading ggl*/
	ggl.resize(ijg[ig[nb]]);
	io_double.read_txt(f_ggl, ggl);

	return true;
}

bool MatrixBlock_Asymmetric::read_matrix_bin(const char * f_ig, const char * f_jg, const char * f_idi, const char * f_ijg, const char * f_di, const char * f_ggu, const char * f_ggl)
{
	if(!size_set)
	{
		cout << "Error: Size of block matrix is not set " << endl;
		return false;
	}

	/*CLASSES FOR INPUT/OUTPUT*/
	IO<double> io_double;
	IO<size_t> io_size_t;

	/*reading idi*/
	idi.resize(nb + 1);
	io_size_t.read_bin(f_idi, idi);
	if(idi[0] != 0)
		array_shift(idi, (size_t)-1);

	/*reading di*/
	di.resize(idi[nb]);
	io_double.read_bin(f_di, di);

	/*reading ig*/
	ig.resize(nb + 1);
	io_size_t.read_bin(f_ig, ig);
	if(ig[0] != 0)
		array_shift(ig, (size_t)-1);

	/*reading jg*/
	jg.resize(ig[nb]);
	io_size_t.read_bin(f_jg, jg);
	if(jg[0] != 0)
		array_shift(jg, (size_t)-1);

	/*reading ijg*/
	ijg.resize(ig[nb] + 1);
	io_size_t.read_bin(f_ijg, ijg);
	if(ijg[0] != 0)
		array_shift(ijg, (size_t)-1);

	/*reading ggu*/
	ggu.resize(ijg[ig[nb]]);
	io_double.read_bin(f_ggu, ggu);

	/*reading ggl*/
	ggl.resize(ijg[ig[nb]]);
	io_double.read_bin(f_ggl, ggl);

	return true;
}

bool MatrixBlock_Asymmetric::read_slae_txt(const char * f_ig, const char * f_jg, const char * f_idi, const char * f_ijg, const char * f_di, const char * f_ggu, const char * f_ggl, const char * f_rp)
{
	read_matrix_txt(f_ig, f_jg, f_idi, f_ijg, f_di, f_ggu, f_ggl);
	read_rp_txt(f_rp);

	return true;
}

bool MatrixBlock_Asymmetric::read_slae_bin(const char * f_ig, const char * f_jg, const char * f_idi, const char * f_ijg, const char * f_di, const char * f_ggu, const char * f_ggl, const char * f_rp)
{
	read_matrix_bin(f_ig, f_jg, f_idi, f_ijg, f_di, f_ggu, f_ggl);
	read_rp_bin(f_rp);
	return false;
}

void MatrixBlock_Asymmetric::values_to_complex()
{
	complex_di.clear();
	complex_ggu.clear();
	complex_ggl.clear();
	complex_rp.clear();

	double_to_complex(idi, idi[nb], di, complex_di);
	double_to_complex(ijg, ijg[ig[nb]], ggu, complex_ggu);
	double_to_complex(ijg, ijg[ig[nb]], ggl, complex_ggl);
	for(size_t i = 0; i < nb; i++)
	{
		complex_rp.push_back(complex<double>(rp[2 * i], rp[2 * i + 1]));
	}
}

void MatrixBlock_Asymmetric::clear()
{
	MatrixBlock::clear();
	ggu.clear();
	ggl.clear();
	complex_ggu.clear();
	complex_ggl.clear();
}