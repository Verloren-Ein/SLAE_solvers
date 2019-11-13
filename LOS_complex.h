#pragma once
#include "solver.h"
#include "common.h"
#include "functions.h"

// for symmetric complex block matrix
class LOS_Complex : public SYMMETRIC_COMPLEX_BLOCK_SOLVER
{
private:
	// residual smoothing
	bool smooth;

	// aux vectors
	vector<complex<double>> x;
	vector<complex<double>> r;
	vector<complex<double>> z;
	vector<complex<double>> p;
	vector<complex<double>> s;
	vector<complex<double>> a;
	vector<complex<double>> w;
	vector<complex<double>> aux;
	vector<complex<double>> sm_x;	//for smoothing
	vector<complex<double>> sm_r;	//for smoothing

public:
	void init(vector<size_t> _ig, vector<size_t> _jg, vector<complex<double>> _di,
		vector<complex<double>> _gg, vector<complex<double>> _rp, size_t _n);
	void solve(const vector<complex<double>> &initial_approximation, const size_t maxiter, const double eps, vector<complex<double>> &solution);
	void set_smooth(bool value);

	LOS_Complex() : SYMMETRIC_COMPLEX_BLOCK_SOLVER()
	{
		name = "LOS_Complex";
	}
	~LOS_Complex() {}

private:
	virtual string create_log_file_name();
	virtual void create_general_info();
};