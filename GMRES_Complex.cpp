#include "GMRES_Complex.h"

void GMRES_Complex::gmres_solve_WITHOUT(const vector<complex<double>>& initial_approximation, const size_t maxiter, const double eps, vector<complex<double>>& solution)
{
	/*LOG-BLOCK*/
	/**/
	this->maxiter = maxiter;
	this->eps = eps;

	//displays information about slae and solver if required
	create_general_info();
	general_information_output();

	auto start = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(start);

	ofstream ostr;
	ostr.precision(CURRENT_DISPLAY_PRECISION);
	cout.precision(CURRENT_DISPLAY_PRECISION);
	if(log_file)
	{
		string log_file_name = create_log_file_name();
		ostr.open(log_file_name);
		ostr << "SLAE calculation started at " << std::ctime(&start_time) << "relative residual:" << endl;
	}
	if(log_console)
	{
		cout << "SLAE calculation started at " << std::ctime(&start_time) << "relative residual:" << endl;
	}

	// solver variables
	double right_part_norm;
	double residual_norm;
	bool found;
	size_t iteration;
	double relative_residual = 0.0;
	double relative_residual_old = 0.0;

	// x - initial approximation
	for(size_t i = 0; i < n; i++)
		x[i] = initial_approximation[i];

	// <vector> r = A*x0
	multiply_matrix_by_vector(x, r);

	// <vector> r = r0 = S^-1*(b - A*x0)
	for(size_t i = 0; i < n; i++)
		r[i] = rp[i] - r[i];

	// calculate right part norm
	right_part_norm = sqrt(dot_prod(rp));

	size_t depth_current = depth;
	size_t stagnation_count = 0;
	found = false;
	for(iteration = 0; iteration < maxiter && !found; iteration++)
	{
		// calculate residual norm
		residual_norm = sqrt(dot_prod(r));

		// calculate relative residual
		relative_residual_old = relative_residual;
		relative_residual = residual_norm / right_part_norm;

		if(log_console)
			if(iteration % iter_shift == 0)
				cout << iteration << " " << relative_residual << endl;
		if(log_file)
			if(iteration % iter_shift == 0)
				ostr << iteration << " " << relative_residual << endl;

		// stagnation check
		if(STAGNATION_CHECK)
			if(abs(relative_residual - relative_residual_old) <= STAGNATION_PRECISION)
			{
				stagnation_count++;
				if(stagnation_count > MAX_STAGNATION_COUNT)
				{
					found = true;
					if(log_console)
					{
						cout << "Warning: Stagnation detected! " << endl;
					}
					if(log_file)
					{
						ostr << "Warning: Stagnation detected! " << endl;
					}
				}
			}
			else
				stagnation_count = 0;
		if(relative_residual > eps && !found)
		{
			// <matrix> col V(0) = rk-1/|rk-1|
			for(size_t i = 0; i < n; i++)
				V(i, 0) = r[i] / residual_norm;

			// <vector> d = dk = (|rk-1|, 0, ..., 0)
			d[0] = residual_norm;
			for(size_t i = 1; i <= depth_current; i++)
				d[i] = 0.0;

			// <matrix> H = set to zero
			H.reset();

			for(size_t myu = 0; myu < depth_current; myu++)
			{
				// <vector> w = S^-1 * A * Q^-1 * V(myu)
				multiply_matrix_by_vector(V.get_col(myu), w);
				
				// <matrix> H(lambda, myu) = V(lambda)*w
				for(size_t lambda = 0; lambda <= myu; lambda++)
					H(lambda, myu) = dot_prod(V.get_col(lambda), w);

				// <vector> w = ~vmyu+1 = w - sum(H(lambda, myu)*V(lambda))
				for(size_t lambda = 0; lambda <= myu; lambda++)
					for(size_t i = 0; i < n; i++)
						w[i] -= H(lambda, myu) * V(i, lambda);

				// <matrix> H(myu+1, myu) = |~Vmyu+1|
				double v_norm = sqrt(dot_prod(w));
				H(myu + 1, myu) = v_norm;

				if(abs(v_norm) < 1e-15)
				{
					if(myu > 1)  depth_current = myu;
					else	depth_current = 1;
					break;
				}

				// <matrix> V(myu+1) = ~Vmyu+1 / H(mu+1, mu)
				if(myu + 1 < depth_current)
					for(size_t i = 0; i < n; i++)
						V(i, myu + 1) = w[i] / v_norm;
			}

			// H to triangular view
			H_to_triangular_view(depth_current);

			// <vector> d = zk
			find_z(depth_current);

			// <vector> x = xk = xk-1 + V*zk
			for(size_t i = 0; i < n; i++)
				for(size_t j = 0; j < depth_current; j++)
					x[i] += V(i, j) * d[j];

			// <vector> r = rk = S^-1*(b - A*Q^-1*xk)
			multiply_matrix_by_vector(x, r);
			for(size_t i = 0; i < n; i++)
				r[i] = rp[i] - r[i];
		}
		else
			found = true;

	}

	/*CALCULATE FINAL RESIDUAL*/
	/**/
	// <vector> s = A*xj
	multiply_matrix_by_vector(x, w);

	// <vector> r = A*x0
	multiply_matrix_by_vector(initial_approximation, r);

	// <vector> s = b - A*xj
	// <vector> r = r0 = b - A*x0
	for(size_t i = 0; i < n; i++)
	{
		w[i] = rp[i] - w[i];
		r[i] = rp[i] - r[i];
	}
	residual_norm = sqrt(dot_prod(w));
	right_part_norm = sqrt(dot_prod(r));

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	if(log_console)
	{
		cout << "finished computation at " << std::ctime(&end_time)
			<< "relative residual (b - Ax): " << residual_norm / right_part_norm << endl
			<< "elapsed time: " << elapsed_seconds.count() << "s\n";;
	}
	if(log_file)
	{
		ostr << "finished computation at " << std::ctime(&end_time)
			<< "relative residual (b - Ax): " << residual_norm / right_part_norm << endl
			<< "elapsed time: " << elapsed_seconds.count() << "s\n";;
		ostr.close();
	}

	for(size_t i = 0; i < n; i++)
		solution[i] = x[i];
}

void GMRES_Complex::gmres_solve_D(const vector<complex<double>>& initial_approximation, const size_t maxiter, const double eps, vector<complex<double>>& solution)
{
	/*LOG-BLOCK*/
	/**/
	this->maxiter = maxiter;
	this->eps = eps;

	//displays information about slae and solver if required
	create_general_info();
	general_information_output();

	auto start = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(start);

	ofstream ostr;
	ostr.precision(CURRENT_DISPLAY_PRECISION);
	cout.precision(CURRENT_DISPLAY_PRECISION);
	if(log_file)
	{
		string log_file_name = create_log_file_name();
		ostr.open(log_file_name);
		ostr << "SLAE calculation started at " << std::ctime(&start_time) << "relative residual:" << endl;
	}
	if(log_console)
	{
		cout << "SLAE calculation started at " << std::ctime(&start_time) << "relative residual:" << endl;
	}

	// solver variables
	double right_part_norm;
	double residual_norm;
	bool found;
	size_t iteration;
	double relative_residual = 0.0;
	double relative_residual_old = 0.0;

	// x - initial approximation
	for(size_t i = 0; i < n; i++)
		x[i] = initial_approximation[i];

	// <vector> r = A*x0
	multiply_matrix_by_vector(x, r);

	// <vector> r = r0 = S^-1*(b - A*x0)
	for(size_t i = 0; i < n; i++)
		r[i] = rp[i] - r[i];
	D_precond(r, r);
	//<vector> w = right part = S^-1*rp
	D_precond(rp, w);

	// calculate right part norm
	right_part_norm = sqrt(dot_prod(w));

	size_t depth_current = depth;
	size_t stagnation_count = 0;
	found = false;
	for(iteration = 0; iteration < maxiter && !found; iteration++)
	{
		// calculate residual norm
		residual_norm = sqrt(dot_prod(r));

		// calculate relative residual
		relative_residual_old = relative_residual;
		relative_residual = residual_norm / right_part_norm;

		if(log_console)
			if(iteration % iter_shift == 0)
				cout << iteration << " " << relative_residual << endl;
		if(log_file)
			if(iteration % iter_shift == 0)
				ostr << iteration << " " << relative_residual << endl;

		// stagnation check
		if(STAGNATION_CHECK)
			if(abs(relative_residual - relative_residual_old) <= STAGNATION_PRECISION)
			{
				stagnation_count++;
				if(stagnation_count > MAX_STAGNATION_COUNT)
				{
					found = true;
					if(log_console)
					{
						cout << "Warning: Stagnation detected! " << endl;
					}
					if(log_file)
					{
						ostr << "Warning: Stagnation detected! " << endl;
					}
				}
			}
			else
				stagnation_count = 0;
		if(relative_residual > eps && !found)
		{
			// <matrix> col V(0) = rk-1/|rk-1|
			for(size_t i = 0; i < n; i++)
				V(i, 0) = r[i] / residual_norm;

			// <vector> d = dk = (|rk-1|, 0, ..., 0)
			d[0] = residual_norm;
			for(size_t i = 1; i <= depth_current; i++)
				d[i] = 0.0;

			// <matrix> H = set to zero
			H.reset();

			for(size_t myu = 0; myu < depth_current; myu++)
			{
				// <vector> w = S^-1 * A * Q^-1 * V(myu)
				multiply_matrix_by_vector(V.get_col(myu), w);
				D_precond(w, w);

				// <matrix> H(lambda, myu) = V(lambda)*w
				for(size_t lambda = 0; lambda <= myu; lambda++)
					H(lambda, myu) = dot_prod(V.get_col(lambda), w);

				// <vector> w = ~vmyu+1 = w - sum(H(lambda, myu)*V(lambda))
				for(size_t lambda = 0; lambda <= myu; lambda++)
					for(size_t i = 0; i < n; i++)
						w[i] -= H(lambda, myu) * V(i, lambda);

				// <matrix> H(myu+1, myu) = |~Vmyu+1|
				double v_norm = sqrt(dot_prod(w));
				H(myu + 1, myu) = v_norm;

				if(abs(v_norm) < 1e-15)
				{
					depth_current = myu;
					break;
				}

				// <matrix> V(myu+1) = ~Vmyu+1 / H(mu+1, mu)
				if(myu + 1 < depth_current)
					for(size_t i = 0; i < n; i++)
						V(i, myu + 1) = w[i] / v_norm;
			}

			// H to triangular view
			H_to_triangular_view(depth_current);

			// <vector> d = zk
			find_z(depth_current);

			// <vector> x = xk = xk-1 + V*zk
			for(size_t i = 0; i < n; i++)
				for(size_t j = 0; j < depth_current; j++)
					x[i] += V(i, j) * d[j];

			// <vector> r = rk = S^-1*(b - A*Q^-1*xk)
			multiply_matrix_by_vector(x, r);
			for(size_t i = 0; i < n; i++)
				r[i] = rp[i] - r[i];
			D_precond(r, r);
		}
		else
			found = true;
	}

	/*CALCULATE FINAL RESIDUAL*/
	/**/
	// <vector> s = A*xj
	multiply_matrix_by_vector(x, w);

	// <vector> r = A*x0
	multiply_matrix_by_vector(initial_approximation, r);

	// <vector> s = b - A*xj
	// <vector> r = r0 = b - A*x0
	for(size_t i = 0; i < n; i++)
	{
		w[i] = rp[i] - w[i];
		r[i] = rp[i] - r[i];
	}
	residual_norm = sqrt(dot_prod(w));
	right_part_norm = sqrt(dot_prod(r));


	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	if(log_console)
	{
		cout << "finished computation at " << std::ctime(&end_time)
			<< "relative residual (b - Ax): " << residual_norm / right_part_norm << endl
			<< "elapsed time: " << elapsed_seconds.count() << "s\n";;
	}
	if(log_file)
	{
		ostr << "finished computation at " << std::ctime(&end_time)
			<< "relative residual (b - Ax): " << residual_norm / right_part_norm << endl
			<< "elapsed time: " << elapsed_seconds.count() << "s\n";;
		ostr.close();
	}

	for(size_t i = 0; i < n; i++)
		solution[i] = x[i];
}

void GMRES_Complex::gmres_solve_LDLT_LLT(const vector<complex<double>>& initial_approximation, const size_t maxiter, const double eps, vector<complex<double>>& solution)
{
	/*LOG-BLOCK*/
	/**/
	this->maxiter = maxiter;
	this->eps = eps;

	//displays information about slae and solver if required
	create_general_info();
	general_information_output();

	auto start = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(start);

	ofstream ostr;
	ostr.precision(CURRENT_DISPLAY_PRECISION);
	cout.precision(CURRENT_DISPLAY_PRECISION);
	if(log_file)
	{
		string log_file_name = create_log_file_name();
		ostr.open(log_file_name);
		ostr << "SLAE calculation started at " << std::ctime(&start_time) << "relative residual:" << endl;
	}
	if(log_console)
	{
		cout << "SLAE calculation started at " << std::ctime(&start_time) << "relative residual:" << endl;
	}

	// solver variables
	double right_part_norm;
	double residual_norm;
	bool found;
	size_t iteration;
	double relative_residual = 0.0;
	double relative_residual_old = 0.0;
	
	// <vector> x = x0 = Q*x0
	mult_Q(initial_approximation, x);
	// <vector> r = A*x0
	multiply_matrix_by_vector(x, r);
	// <vector> r = r0 = S^-1*(b - A*x0)
	for(size_t i = 0; i < n; i++)
		aux[i] = rp[i] - r[i];
	solve_S(aux, r);
	//<vector> w = right part = S^-1*rp
	solve_S(rp, w);

	// calculate right part norm
	right_part_norm = sqrt(dot_prod(w));

	size_t depth_current = depth;
	size_t stagnation_count = 0;
	found = false;
	for(iteration = 0; iteration < maxiter && !found; iteration++)
	{
		// calculate residual norm
		residual_norm = sqrt(dot_prod(r));

		// calculate relative residual
		relative_residual_old = relative_residual;
		relative_residual = residual_norm / right_part_norm;

		if(log_console)
			if(iteration % iter_shift == 0)
				cout << iteration << " " << relative_residual << endl;
		if(log_file)
			if(iteration % iter_shift == 0)
				ostr << iteration << " " << relative_residual << endl;

		// stagnation check
		if(STAGNATION_CHECK)
			if(abs(relative_residual - relative_residual_old) / relative_residual < STAGNATION_PRECISION)
			{
				stagnation_count++;
				if(stagnation_count > MAX_STAGNATION_COUNT)
				{
					found = true;
					if(log_console)
					{
						cout << "Warning: Stagnation detected! " << endl;
					}
					if(log_file)
					{
						ostr << "Warning: Stagnation detected! " << endl;
					}
				}
			}
			else
				stagnation_count = 0;
		if(relative_residual > eps && !found)
		{
			// <matrix> col V(0) = rk-1/|rk-1|
			for(size_t i = 0; i < n; i++)
				V(i, 0) = r[i] / residual_norm;

			// <vector> d = dk = (|rk-1|, 0, ..., 0)
			d[0] = residual_norm;
			for(size_t i = 1; i <= depth_current; i++)
				d[i] = 0.0;

			// <matrix> H = set to zero
			//H.reset();

			for(size_t myu = 0; myu < depth_current; myu++)
			{
				// <vector> w = S^-1 * A * Q^-1 * V(myu)
				solve_QAS(V.get_col(myu), w);

				// <matrix> H(lambda, myu) = V(lambda)*w
				for(size_t lambda = 0; lambda <= myu; lambda++)
					H(lambda, myu) = dot_prod(V.get_col(lambda), w);

				// <vector> w = ~vmyu+1 = w - sum(H(lambda, myu)*V(lambda))
				for(size_t lambda = 0; lambda <= myu; lambda++)
					for(size_t i = 0; i < n; i++)
						w[i] -= H(lambda, myu) * V(i, lambda);

				// <matrix> H(myu+1, myu) = |~Vmyu+1|
				double v_norm = sqrt(dot_prod(w));
				H(myu + 1, myu) = v_norm;

				if(abs(v_norm) < 1e-15)
				{
					depth_current = myu;
					break;
				}

				// <matrix> V(myu+1) = ~Vmyu+1 / H(mu+1, mu)
				if(myu + 1 < depth_current)
					for(size_t i = 0; i < n; i++)
						V(i, myu + 1) = w[i] / v_norm;
			}

			// H to triangular view
			H_to_triangular_view(depth_current);

			// <vector> d = zk
			find_z(depth_current);

			// <vector> x = xk = xk-1 + V*zk
			for(size_t i = 0; i < n; i++)
				for(size_t j = 0; j < depth_current; j++)
					x[i] += V(i, j) * d[j];

			// <vector> r = rk = S^-1*(b - A*Q^-1*xk)
			for(size_t i = 0; i < n; i++)
				aux[i] = x[i];
			solve_Q(aux, aux);
			multiply_matrix_by_vector(aux, r);
			for(size_t i = 0; i < n; i++)
				aux[i] = rp[i] - r[i];
			solve_S(aux, r);
		}
		else
			found = true;
	}

	// <vector> solution = Q^-1*x
	solve_Q(x, solution);

	/*CALCULATE FINAL RESIDUAL*/
	/**/
	// <vector> s = A*xj
	multiply_matrix_by_vector(solution, w);

	// <vector> r = A*x0
	multiply_matrix_by_vector(initial_approximation, r);

	// <vector> s = b - A*xj
	// <vector> r = r0 = b - A*x0
	for(size_t i = 0; i < n; i++)
	{
		w[i] = rp[i] - w[i];
		r[i] = rp[i] - r[i];
	}
	residual_norm = sqrt(dot_prod(w));
	right_part_norm = sqrt(dot_prod(r));


	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	if(log_console)
	{
		cout << "finished computation at " << std::ctime(&end_time)
			<< "relative residual (b - Ax): " << residual_norm / right_part_norm << endl
			<< "elapsed time: " << elapsed_seconds.count() << "s\n";;
	}
	if(log_file)
	{
		ostr << "finished computation at " << std::ctime(&end_time)
			<< "relative residual (b - Ax): " << residual_norm / right_part_norm << endl
			<< "elapsed time: " << elapsed_seconds.count() << "s\n";;
		ostr.close();
	}
}

void GMRES_Complex::H_to_triangular_view(size_t depth_current)
{
	for(size_t i = 0; i < depth_current; i++)
	{
		double val, val2, val3;
		val = H(i + 1, i).real();
		val2 = abs(H(i, i));
		val3 = 1.0 / sqrt(val * val + val2 * val2);
		double s = val * val3;
		complex<double> c = H(i, i) * val3;

		for(size_t j = i; j < depth_current; j++)
		{
			complex<double> H_ij = H(i, j);
			H(i, j) = H_ij * conj(c) + H(i + 1, j) * s;
			H(i + 1, j) = -H_ij * s + H(i + 1, j) * c;
		}
		complex<double> d_i = d[i];
		d[i] = d_i * conj(c) + d[i + 1] * s;
		d[i + 1] = -d_i * s + d[i + 1] * c;
	}
}

void GMRES_Complex::find_z(size_t depth_current)
{
	for(size_t i = 0; i < depth_current; i++)
	{
		d[i] = d[i] / H(i, i);
		for(size_t j = i + 1; j < depth_current; j++)
			H(i, j) /= H(i, i);
	}
	for(size_t i = depth_current - 1; i > 0; i--)
		for(size_t j = i - 1, ju = i; ju > 0; j--, ju--)
			d[j] -= H(j, i) * d[i];
}

void GMRES_Complex::init(vector<size_t> _ig, vector<size_t> _jg, vector<complex<double>> _di, vector<complex<double>> _gg, vector<complex<double>> _rp, size_t _n)
{
	ig = _ig;
	jg = _jg;
	di = _di;
	gg = _gg;
	rp = _rp;
	n = _n;

	x.clear();
	r.clear();
	w.clear();
	d.clear();
	aux.clear();
	L_di.clear();
	L_gg.clear();
	H.clear();
	V.clear();

	x.resize(n);
	r.resize(n);
	w.resize(n);
	d.resize(depth + 1);
	aux.resize(n);
	L_di.resize(n);
	L_gg.resize(gg.size());
	H = Matrix<complex<double>>(depth + 1, depth);
	V = Matrix<complex<double>>(n, depth);

	for(size_t i = 0; i < n; i++)
		L_di[i] = di[i];
	for(size_t i = 0; i < gg.size(); i++)
		L_gg[i] = gg[i];
}

void GMRES_Complex::solve(const vector<complex<double>>& initial_approximation, const size_t maxiter, const double eps, vector<complex<double>>& solution)
{
	switch(precond_method)
	{
	case WITHOUT:
		gmres_solve_WITHOUT(initial_approximation, maxiter, eps, solution);
		break;
	case D:
		gmres_solve_D(initial_approximation, maxiter, eps, solution);
		break;
	case LDLT:
		LDLT_decomposition();
		gmres_solve_LDLT_LLT(initial_approximation, maxiter, eps, solution);
		break;
	case LLT:
		LLT_decomposition();
		gmres_solve_LDLT_LLT(initial_approximation, maxiter, eps, solution);
		break;
	default:
		cout << "Error: preconditioner for number " << precond_method << " does not exist" << endl;
		break;
	}
}

void GMRES_Complex::set_depth(size_t depth)
{
	this->depth = depth;
}

string GMRES_Complex::create_log_file_name()
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

	string depth_str;
	number_to_string(depth_str, depth);

	log_file_name += "_depth=";
	log_file_name += depth_str;

	log_file_name += ".txt";

	return log_file_name;
}

void GMRES_Complex::create_general_info()
{
	SYMMETRIC_COMPLEX_BLOCK_SOLVER::create_general_info();
	general_solver_information += "depth = ";
	string depth_str;
	number_to_string(depth_str, depth);
	general_solver_information += depth_str;
	general_solver_information += "\n";
	general_solver_information += "________________________________\n";
}
