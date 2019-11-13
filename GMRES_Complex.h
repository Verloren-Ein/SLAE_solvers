#pragma once
#include "solver.h"
#include "common.h"
#include "functions.h"
#include "Matrix.h"

#define DEFAULT_GMRES_DEPTH 5

// for symmetric complex block matrix
class GMRES_Complex : public SYMMETRIC_COMPLEX_BLOCK_SOLVER
{
private:
	// solver options
	size_t depth;

	// aux vectors
	vector<complex<double>> x;
	vector<complex<double>> r;
	vector<complex<double>> w;
	vector<complex<double>> d;
	vector<complex<double>> aux;

	// aux matrices
	Matrix<complex<double>> H;
	Matrix<complex<double>> V;

	void gmres_solve_WITHOUT(const vector<complex<double>> &initial_approximation, const size_t maxiter, const double eps, vector<complex<double>> &solution);
	void gmres_solve_D(const vector<complex<double>> &initial_approximation, const size_t maxiter, const double eps, vector<complex<double>> &solution);
	void gmres_solve_LDLT_LLT(const vector<complex<double>> &initial_approximation, const size_t maxiter, const double eps, vector<complex<double>> &solution);

	void H_to_triangular_view(size_t depth_current);
	void find_z(size_t depth_current);

public:
	void init(vector<size_t> _ig, vector<size_t> _jg, vector<complex<double>> _di,
		vector<complex<double>> _gg, vector<complex<double>> _rp, size_t _n);
	void solve(const vector<complex<double>> &initial_approximation, const size_t maxiter, const double eps, vector<complex<double>> &solution);
	void set_depth(size_t depth);

	GMRES_Complex() : SYMMETRIC_COMPLEX_BLOCK_SOLVER()
		,depth(DEFAULT_GMRES_DEPTH)
	{
		name = "GMRES_Complex";
	}
	~GMRES_Complex() {}

private:
	virtual string create_log_file_name();
	virtual void create_general_info();
};