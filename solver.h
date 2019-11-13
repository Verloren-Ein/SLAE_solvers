#pragma once
#include "common.h"
#include "functions.h"

/*preconditioner`s constants*/
#define WITHOUT 0
#define D		1
#define SSOR	2
#define LDLT	3
#define LLT		4

/*for log*/
#define GEN_INFO_LOG_FILE "gen_inf.txt"
#define DEFAULT_OUTPUT_ITERATIVE_SHIFT 1

/*for stagnation check*/
#define STAGNATION_CHECK true
#define MAX_STAGNATION_COUNT 20
#define STAGNATION_PRECISION 1e-7

//iterative solver for sparse symmetric complex block matrices
class SYMMETRIC_COMPLEX_BLOCK_SOLVER
{
protected:
	// solver options
	string name;	//method`s name
	size_t maxiter;
	double eps;

	// matrix data
	size_t n;
	vector<size_t> ig;
	vector<size_t> jg;

	vector<complex<double>> di;
	vector<complex<double>> gg;
	vector<complex<double>> rp;

	// aux vectors
	vector<complex<double>> L_di;	//for decomposition
	vector<complex<double>> L_gg;	//for decomposition
	
	void multiply_matrix_by_vector(const vector<complex<double>> &vec, vector<complex<double>> &result);
	void D_precond(const vector<complex<double>> &vec, vector<complex<double>> &result);
	void SSOR_precond(const vector<complex<double>> &vec, vector<complex<double>> &result);
	void LDLT_precond(const vector<complex<double>> &vec, vector<complex<double>> &result);
	void LLT_precond(const vector<complex<double>> &vec, vector<complex<double>> &result);
	bool LDLT_decomposition();
	bool LLT_decomposition();

	// S - lower; Q - upper
	void solve_S(const vector<complex<double>> &vec, vector<complex<double>> &result);
	void solve_Q(vector<complex<double>> &vec, vector<complex<double>> &result);
	void solve_QAS(const vector<complex<double>> &vec, vector<complex<double>> &result);
	void mult_Q(const vector<complex<double>> &vec, vector<complex<double>> &result);
	//void direct_move(const vector<complex<double>> &vec, vector<complex<double>> &result);
	//void reverse_move(const vector<complex<double>> &vec, vector<complex<double>> &result);

public:
	virtual void init(vector<size_t> _ig, vector<size_t> _jg, vector<complex<double>> _di,
		vector<complex<double>> _gg, vector<complex<double>> _rp, size_t _n) = 0;
	virtual void solve(const vector<complex<double>> &initial_approximation, 
		const size_t maxiter, const double eps, vector<complex<double>> &solution) = 0;
	void set_output_options(bool _log_console, bool _log_file, bool _gen_inf_console, bool _gen_inf_file, size_t _iter_shift);

	/*0 - without; 1 - diagonal; 2 - SSOR; 3 - LDLT; 4 - LLT;*/
	void set_preconditioner(size_t method);

	SYMMETRIC_COMPLEX_BLOCK_SOLVER()
		:log_console(false)
		, log_file(false)
		, gen_inf_console(false)
		, gen_inf_file(false)
		, precond_method(WITHOUT)
		, iter_shift(DEFAULT_OUTPUT_ITERATIVE_SHIFT)
	{}
	virtual ~SYMMETRIC_COMPLEX_BLOCK_SOLVER() = 0
	{}

protected:
	// print information about the residual to the console
	bool log_console;
	// print information about the residual to the log-file
	bool log_file;
	// print common information about solver and slae to the console
	bool gen_inf_console;
	// print common information about solver and slae to the file
	bool gen_inf_file;

	/*0 - without; 1 - diagonal; 2 - SSOR; 3 - LDLT; 4 - LLT;*/
	size_t precond_method;
	size_t iter_shift;

	string general_solver_information;
	string general_slae_information;

	void general_information_output();
	virtual string create_log_file_name();
	virtual void create_general_info();
};