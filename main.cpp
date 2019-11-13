#include "common.h"
#include "IO.h"
#include "MatrixBlock.h"
#include "COCG.h"
#include "COCR.h"
#include "GMRES_Complex.h"
#include "LOS_complex.h"

#define COCG_SOLVER		0
#define COCR_SOLVER		1
#define GMRES_SOLVER	2
#define LOS_SOLVER		3

int main()
{
	/*variables*/
	bool testing = false;
	string using_matrix = "/3";
	size_t using_solver = GMRES_SOLVER;
	size_t using_preconditioner = LLT;
	bool using_smooth = false;
	size_t gmres_depth = 10;

	MatrixBlock_Symmetric *A;
	GMRES_Complex *GMRES_solver;
	COCG *COCG_solver;
	COCR *COCR_solver;
	LOS_Complex *LOS_solver;

	string folder;

	size_t n;
	size_t maxiter;
	double eps;

	vector<complex<double>> solution;
	vector<complex<double>> initial_approximation;

	CURRENT_DISPLAY_PRECISION = 16;

	if(testing)
	{
		folder = "matrices/test";
		read_n_eps_maxiter(string(folder + "/kuslau.txt").c_str(), &n, &eps, &maxiter);

		A = new MatrixBlock_Symmetric(n);
		A->read_slae_txt(string(folder + "/ig.txt").c_str(), string(folder + "/jg.txt").c_str(), string(folder + "/idi.txt").c_str(),
			string(folder + "/ijg.txt").c_str(), string(folder + "/di.txt").c_str(), string(folder + "/gg.txt").c_str(), string(folder + "/pr.txt").c_str());
		A->values_to_complex();

		initial_approximation.resize(A->nb);
		solution.resize(A->nb);

		//solver = new COCR();
		//solver->set_output_options(true, true, true, true, 1);
		//solver->set_preconditioner(0);
		//solver->set_smooth(true);
		//solver->init(A->ig, A->jg, A->complex_di, A->complex_gg, A->complex_rp, A->nb);
		//solver->solve(initial_approximation, maxiter, eps, solution);

		//GMRES_solver = new GMRES_Complex();
		//GMRES_solver->set_depth(gmres_depth);
		//GMRES_solver->set_output_options(true, true, true, true, 1);
		//GMRES_solver->set_preconditioner(using_preconditioner);
		//GMRES_solver->init(A->ig, A->jg, A->complex_di, A->complex_gg, A->complex_rp, A->nb);
		//GMRES_solver->solve(initial_approximation, maxiter, eps, solution);

		LOS_solver = new LOS_Complex();
		LOS_solver->set_output_options(true, true, true, true, 1);
		LOS_solver->set_preconditioner(using_preconditioner);
		LOS_solver->init(A->ig, A->jg, A->complex_di, A->complex_gg, A->complex_rp, A->nb);
		LOS_solver->solve(initial_approximation, maxiter, eps, solution);

	}
	else
	{
		folder = "matrices" + using_matrix;
		read_n_eps_maxiter(string(folder + "/kuslau").c_str(), &n, &eps, &maxiter);

		A = new MatrixBlock_Symmetric(n);
		A->read_slae_bin(string(folder + "/ig").c_str(), string(folder + "/jg").c_str(), string(folder + "/idi").c_str(),
			string(folder + "/ijg").c_str(), string(folder + "/di").c_str(), string(folder + "/gg").c_str(), string(folder + "/pr").c_str());
		A->values_to_complex();

		initial_approximation.resize(A->nb);
		solution.resize(A->nb);

		switch(using_solver)
		{
		case COCG_SOLVER:
			COCG_solver = new COCG();
			COCG_solver->set_output_options(true, true, true, true, 1);
			COCG_solver->set_preconditioner(using_preconditioner);
			COCG_solver->set_smooth(using_smooth);
			COCG_solver->init(A->ig, A->jg, A->complex_di, A->complex_gg, A->complex_rp, A->nb);
			COCG_solver->solve(initial_approximation, maxiter, eps, solution);
			break;
		case COCR_SOLVER:
			COCR_solver = new COCR();
			COCR_solver->set_output_options(true, true, true, true, 1);
			COCR_solver->set_preconditioner(using_preconditioner);
			COCR_solver->set_smooth(using_smooth);
			COCR_solver->init(A->ig, A->jg, A->complex_di, A->complex_gg, A->complex_rp, A->nb);
			COCR_solver->solve(initial_approximation, maxiter, eps, solution);
			break;
		case GMRES_SOLVER:
			GMRES_solver = new GMRES_Complex();
			GMRES_solver->set_depth(gmres_depth);
			GMRES_solver->set_output_options(true, true, true, true, 1);
			GMRES_solver->set_preconditioner(using_preconditioner);
			GMRES_solver->init(A->ig, A->jg, A->complex_di, A->complex_gg, A->complex_rp, A->nb);
			GMRES_solver->solve(initial_approximation, maxiter, eps, solution);
			break;
		case LOS_SOLVER:
			LOS_solver = new LOS_Complex();
			LOS_solver->set_output_options(true, true, true, true, 1);
			LOS_solver->set_preconditioner(using_preconditioner);
			LOS_solver->set_smooth(using_smooth);
			LOS_solver->init(A->ig, A->jg, A->complex_di, A->complex_gg, A->complex_rp, A->nb);
			LOS_solver->solve(initial_approximation, maxiter, eps, solution);
		default:
			break;
		}
	}
	
	return 0;
}
