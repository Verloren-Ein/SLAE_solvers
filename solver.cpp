#include "solver.h"

void SYMMETRIC_COMPLEX_BLOCK_SOLVER::init(vector<size_t> _ig, vector<size_t> _jg, vector<complex<double>> _di, vector<complex<double>> _gg, vector<complex<double>> _rp, size_t _n)
{
	ig = _ig;
	jg = _jg;
	di = _di;
	gg = _gg;
	rp = _rp;
	n = _n;

	L_di.clear();
	L_gg.clear();

	L_di.resize(n);
	L_gg.resize(gg.size());

	for(size_t i = 0; i < n; i++)
		L_di[i] = di[i];
	for(size_t i = 0; i < gg.size(); i++)
		L_gg[i] = gg[i];
}

void SYMMETRIC_COMPLEX_BLOCK_SOLVER::multiply_matrix_by_vector(const vector<complex<double>> &vec, vector<complex<double>> &result)
{
	for(size_t i = 0; i < n; i++)
	{
		complex<double> v_el = vec[i];
		result[i] = di[i] * v_el;
		for(size_t k = ig[i], k1 = ig[i + 1]; k < k1; k++)
		{
			size_t j = jg[k];
			result[i] += gg[k] * vec[j];
			result[j] += gg[k] * v_el;
		}
	}
}

void SYMMETRIC_COMPLEX_BLOCK_SOLVER::D_precond(const vector<complex<double>>& vec, vector<complex<double>>& result)
{
	for(size_t i = 0; i < n; i++)
		result[i] = vec[i] / di[i];
}

void SYMMETRIC_COMPLEX_BLOCK_SOLVER::SSOR_precond(const vector<complex<double>>& vec, vector<complex<double>>& result)
{
	for(size_t k = 1; k <= n; k++)
	{
		complex<double> s = 0;
		for(size_t i = ig[k - 1]; i < ig[k]; i++)
			s += result[jg[i]] * L_gg[i];
		result[k - 1] = (vec[k - 1] - s) / L_di[k - 1];
	}

	for(size_t i = 0; i < n; i++)
		result[i] *= di[i];

	for(size_t k = n; k > 0; k--)
	{
		result[k - 1] = result[k - 1] / L_di[k - 1];
		for(size_t i = ig[k - 1]; i < ig[k]; i++)
			result[jg[i]] -= result[k - 1] * L_gg[i];
	}
}

void SYMMETRIC_COMPLEX_BLOCK_SOLVER::LDLT_precond(const vector<complex<double>>& vec, vector<complex<double>>& result)
{
	solve_S(vec, result);
	solve_Q(result, result);
}

void SYMMETRIC_COMPLEX_BLOCK_SOLVER::LLT_precond(const vector<complex<double>>& vec, vector<complex<double>>& result)
{
	solve_S(vec, result);
	solve_Q(result, result);
}

bool SYMMETRIC_COMPLEX_BLOCK_SOLVER::LDLT_decomposition()
{
	complex<double> sum_d, sum_l;

	for(size_t k = 0; k < n; k++)
	{
		sum_d = 0;
		size_t i_s = ig[k], i_e = ig[k + 1];

		for(size_t i = i_s; i < i_e; i++)
		{
			sum_l = 0;
			size_t j_s = ig[jg[i]], j_e = ig[jg[i] + 1];

			for(size_t m = i_s; m < i; m++)
			{
				for(size_t j = j_s; j < j_e; j++)
				{
					if(jg[m] == jg[j])
					{
						sum_l += L_gg[m] * L_gg[j] * L_di[jg[m]];
						j_s++;
					}
				}
			}
			L_gg[i] = (L_gg[i] - sum_l) / L_di[jg[i]];
			sum_d += L_gg[i] * L_gg[i] * L_di[jg[i]];
		}
		L_di[k] = L_di[k] - sum_d;
	}
	return true;
}

bool SYMMETRIC_COMPLEX_BLOCK_SOLVER::LLT_decomposition()
{
	complex<double> sum_d, sum_l;

	for(size_t k = 0; k < n; k++)
	{
		sum_d = 0;
		size_t i_s = ig[k], i_e = ig[k + 1];

		for(size_t i = i_s; i < i_e; i++)
		{
			sum_l = 0;
			size_t j_s = ig[jg[i]], j_e = ig[jg[i] + 1];

			for(size_t m = i_s; m < i; m++)
			{
				for(size_t j = j_s; j < j_e; j++)
				{
					if(jg[m] == jg[j])
					{
						sum_l += L_gg[m] * L_gg[j];
						j_s++;
					}
				}
			}
			L_gg[i] = (L_gg[i] - sum_l) / L_di[jg[i]];

			sum_d += L_gg[i] * L_gg[i];
		}
		L_di[k] = sqrt(L_di[k] - sum_d);
	}

	return true;
}

void SYMMETRIC_COMPLEX_BLOCK_SOLVER::solve_S(const vector<complex<double>>& vec, vector<complex<double>>& result)
{
	switch(precond_method)
	{
	case LDLT:
		for(size_t k = 1; k <= n; k++)
		{
			complex<double> s = 0;
			for(size_t i = ig[k - 1]; i < ig[k]; i++)
				s += result[jg[i]] * L_gg[i];
			result[k - 1] = (vec[k - 1] - s);
		}

		for(size_t i = 0; i < n; i++)
			result[i] /= L_di[i];
		break;
	case LLT:
		for(size_t k = 1; k <= n; k++)
		{
			complex<double> s = 0;
			for(size_t i = ig[k - 1]; i < ig[k]; i++)
				s += result[jg[i]] * L_gg[i];
			result[k - 1] = (vec[k - 1] - s) / L_di[k - 1];
		}
		break;
	default:
		break;
	}
}

void SYMMETRIC_COMPLEX_BLOCK_SOLVER::solve_Q(vector<complex<double>>& vec, vector<complex<double>>& result)
{
	switch(precond_method)
	{
	case LDLT:
		for(size_t k = n; k > 0; k--)
		{
			result[k - 1] = vec[k - 1];
			for(size_t i = ig[k - 1]; i < ig[k]; i++)
				vec[jg[i]] -= result[k - 1] * L_gg[i];
		}
		break;
	case LLT:
		for(size_t k = n; k > 0; k--)
		{
			result[k - 1] = vec[k - 1] / L_di[k - 1];
			for(size_t i = ig[k - 1]; i < ig[k]; i++)
				vec[jg[i]] -= result[k - 1] * L_gg[i];
		}
		break;
	default:
		break;
	}
}

void SYMMETRIC_COMPLEX_BLOCK_SOLVER::solve_QAS(const vector<complex<double>>& vec, vector<complex<double>>& result)
{
	vector<complex<double>> tmp = vector<complex<double>>(n);
	for(size_t i = 0; i < n; i++)
		result[i] = vec[i];
	solve_Q(result, result);
	multiply_matrix_by_vector(result, tmp);
	solve_S(tmp, result);
	tmp.clear();
}

void SYMMETRIC_COMPLEX_BLOCK_SOLVER::mult_Q(const vector<complex<double>>& vec, vector<complex<double>>& result)
{
	for(size_t i = 0; i < n; i++)
	{
		complex<double> v_el = vec[i];
		result[i] = L_di[i] * v_el;
		for(size_t k = ig[i], k1 = ig[i + 1]; k < k1; k++)
			result[jg[k]] += L_gg[k] * v_el;
	}
}

void SYMMETRIC_COMPLEX_BLOCK_SOLVER::set_output_options(bool _log_console, bool _log_file, bool _gen_inf_console, bool _gen_inf_file, size_t _iter_shift)
{
	log_console = _log_console;
	log_file = _log_file;
	gen_inf_console = _gen_inf_console;
	gen_inf_file = _gen_inf_file;
	iter_shift = _iter_shift;
}

void SYMMETRIC_COMPLEX_BLOCK_SOLVER::set_preconditioner(size_t method)
{
	precond_method = method;
}

void SYMMETRIC_COMPLEX_BLOCK_SOLVER::general_information_output()
{
	if(gen_inf_console)
	{
		cout << general_solver_information;
		cout << general_slae_information;
	}

	if(gen_inf_file)
	{
		ofstream ostr;
		ostr.open(GEN_INFO_LOG_FILE);
		ostr << general_solver_information;
		ostr << general_slae_information;
		ostr.close();
	}
}

string SYMMETRIC_COMPLEX_BLOCK_SOLVER::create_log_file_name()
{
	string log_file_name = "log_";
	log_file_name += name;

	switch(precond_method)
	{
	case WITHOUT:
		log_file_name += "_WITHOUT";
		break;
	case D:
		log_file_name += "_D";
		break;
	case SSOR:
		log_file_name += "_SSOR";
		break;
	case LDLT:
		log_file_name += "_LDLT";
		break;
	case LLT:
		log_file_name += "_LLT";
		break;
	default:
		break;
	}

	return log_file_name;
}

void SYMMETRIC_COMPLEX_BLOCK_SOLVER::create_general_info()
{
	general_solver_information = "________________________________\n";
	general_solver_information += "*** ";
	general_solver_information += name;
	general_solver_information += " general information ***\n";
	general_solver_information += "preconditioning method: ";

	switch(precond_method)
	{
	case WITHOUT:
		general_solver_information += "without preconditioning\n";
		break;
	case D:
		general_solver_information += "diagonal (Jacobi method)\n";
		break;
	case SSOR:
		general_solver_information += "SSOR (Gauss-Seidel method)\n";
		break;
	case LDLT:
		general_solver_information += "LDLT\n";
		break;
	case LLT:
		general_solver_information += "LLT\n";
		break;
	default:
		break;
	}
	general_solver_information += "maximum number of iterations: ";
	string num;
	number_to_string(num, maxiter);
	general_solver_information += num;
	general_solver_information += "\n";
	general_solver_information += "accuracy: ";
	number_to_string(num, eps);
	general_solver_information += num;
	general_solver_information += "\n";
	
	general_slae_information = "________________________________\n";
	general_slae_information += "*** SLAE information ***\n";
	general_slae_information += "size: ";
	number_to_string(num, n);
	general_slae_information += num;
	general_slae_information += "\n";
	general_slae_information += "number of blocks: ";
	number_to_string(num, 2 * gg.size() + di.size());
	general_slae_information += num;
	general_slae_information += "\n";
	general_slae_information += "________________________________\n";

}
