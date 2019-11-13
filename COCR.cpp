#include"COCR.h"

void COCR::init(vector<size_t> _ig, vector<size_t> _jg, vector<complex<double>> _di, vector<complex<double>> _gg, vector<complex<double>> _rp, size_t _n)
{
	ig = _ig;
	jg = _jg;
	di = _di;
	gg = _gg;
	rp = _rp;
	n = _n;

	x.clear();
	r.clear();
	p.clear();
	z.clear();
	s.clear();
	w.clear();
	a.clear();
	L_di.clear();
	L_gg.clear();
	sm_x.clear();
	sm_r.clear();

	x.resize(n);
	r.resize(n);
	p.resize(n);
	z.resize(n);
	s.resize(n);
	w.resize(n);
	a.resize(n);
	L_di.resize(n);
	L_gg.resize(gg.size());
	sm_x.resize(n);
	sm_r.resize(n);

	for(size_t i = 0; i < n; i++)
		L_di[i] = di[i];
	for(size_t i = 0; i < gg.size(); i++)
		L_gg[i] = gg[i];
}

void COCR::solve(const vector<complex<double>>& initial_approximation, const size_t maxiter, const double eps, vector<complex<double>>& solution)
{
	if(precond_method > 4)
	{
		cout << "Error: preconditioner for number " << precond_method << " does not exist" << endl;
		return;
	}

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
	complex<double> alpha;
	complex<double> beta;
	double eta;
	complex<double> product_1;
	complex<double> product_2;
	double right_part_norm;
	double residual_norm;
	bool found;
	size_t iteration;
	double relative_residual = 0.0;
	double relative_residual_old = 0.0;

	// x - initial approximation
	for(size_t i = 0; i < n; i++)
		sm_x[i] = x[i] = initial_approximation[i];
	// <vector> r = A*x0
	multiply_matrix_by_vector(x, r);

	// <vector> r = r0 = b - A*x0
	for(size_t i = 0; i < n; i++)
		sm_r[i] = r[i] = rp[i] - r[i];

	// <vector> p = p0 = s0 = M^-1*r0
	// <vector> a = a0 = z0 = A*s0
	// <vector> w = w0 = M^-1*z0
	switch(precond_method)
	{
	case WITHOUT:
		for(size_t i = 0; i < n; i++)
			s[i] = r[i];
		multiply_matrix_by_vector(s, z);
		for(size_t i = 0; i < n; i++)
			w[i] = z[i];
		break;
	case D:
		D_precond(r, s);
		multiply_matrix_by_vector(s, z);
		D_precond(z, w);
		break;
	case SSOR:
		SSOR_precond(r, s);
		multiply_matrix_by_vector(s, z);
		SSOR_precond(z, w);
		break;
	case LDLT:
		LDLT_decomposition();
		LDLT_precond(r, s);
		multiply_matrix_by_vector(s, z);
		LDLT_precond(z, w);
		break;
	case LLT:
		LLT_decomposition();
		LLT_precond(r, s);
		multiply_matrix_by_vector(s, z);
		LLT_precond(z, w);
		break;
	default:
		break;
	}

	for(size_t i = 0; i < n; i++)
	{
		p[i] = s[i];
		a[i] = z[i];
	}

	// calculate right part norm
	right_part_norm = sqrt(dot_prod(rp));

	// <double> product_1 = (`aj, sj)
	product_1 = conjugate_prod(a, s);

	size_t stagnation_count = 0;
	found = false;
	for(iteration = 0; iteration < maxiter && !found; iteration++)
	{
		// calculate residual norm
		if(smooth)
			residual_norm = sqrt(dot_prod(sm_r));
		else
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
			if(abs(relative_residual - relative_residual_old) / relative_residual <= STAGNATION_PRECISION)
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
			// alpha = alphaj = (`aj, sj)/(`zj, wj)
			alpha = product_1 / conjugate_prod(z, w);

			// <vector> x = xj+1 = xj + alphaj*pj
			// <vector> r = rj+1 = rj - alphaj*zj
			// <vector> s = sj+1 = sj - alphaj*wj
			for(size_t i = 0; i < n; i++)
			{
				x[i] += alpha * p[i];
				r[i] -= alpha * z[i];
				s[i] -= alpha * w[i];
				w[i] = r[i] - sm_r[i];	//for smoothing
			}

			if(smooth)
			{
				eta = -dot_prod_2(sm_r, w) / dot_prod(w);
				if(eta < 0.0)
					eta = 0.0;
				else if(eta > 1.0)
					eta = 1.0;
				for(size_t i = 0; i < n; i++)
				{
					sm_x[i] = (1 - eta)*sm_x[i] + eta * x[i];
					sm_r[i] = (1 - eta)*sm_r[i] + eta * r[i];
				}
			}

			// <vector> a = aj+1 = A*sj+1
			multiply_matrix_by_vector(s, a);

			// <double> product_2 = (`aj+1, sj+1)
			product_2 = conjugate_prod(a, s);

			// beta = betaj = (`aj+1, sj+1)/(`aj, sj)
			beta = product_2 / product_1;

			product_1 = product_2;

			// <vector> p = pj+1 = sj+1 + betaj*pj
			// <vector> z = zj+1 = aj+1 + betaj*zj
			for(size_t i = 0; i < n; i++)
			{
				p[i] = s[i] + beta * p[i];
				z[i] = a[i] + beta * z[i];
			}

			// <vector> w = wj+1 = M^-1*zj+1
			switch(precond_method)
			{
			case WITHOUT:
				for(size_t i = 0; i < n; i++)
					w[i] = z[i];
				break;
			case D:
				D_precond(z, w);
				break;
			case SSOR:
				SSOR_precond(z, w);
				break;
			case LDLT:
				LDLT_precond(z, w);
				break;
			case LLT:
				LLT_precond(z, w);
				break;
			default:
				break;
			}
		}
		else
			found = true;
	}

	/*CALCULATE FINAL RESIDUAL*/
	/**/
	// <vector> s = A*xj
	if(smooth)
		multiply_matrix_by_vector(sm_x, s);
	else
		multiply_matrix_by_vector(x, s);

	// <vector> r = A*x0
	multiply_matrix_by_vector(initial_approximation, r);

	// <vector> s = b - A*xj
	// <vector> r = r0 = b - A*x0
	for(size_t i = 0; i < n; i++)
	{
		s[i] = rp[i] - s[i];
		r[i] = rp[i] - r[i];
	}
	residual_norm = sqrt(dot_prod(s));
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

	if(smooth)
		for(size_t i = 0; i < n; i++)
			solution[i] = sm_x[i];
	else
		for(size_t i = 0; i < n; i++)
			solution[i] = x[i];

}

void COCR::set_smooth(bool value)
{
	smooth = value;
}

string COCR::create_log_file_name()
{
	string log_file_name = SYMMETRIC_COMPLEX_BLOCK_SOLVER::create_log_file_name();
	if(smooth)
		log_file_name += "_sm";
	log_file_name += ".txt";

	return log_file_name;
}

void COCR::create_general_info()
{
	SYMMETRIC_COMPLEX_BLOCK_SOLVER::create_general_info();

	general_solver_information += "residual smoothing:";
	if(smooth)
		general_solver_information += " yes\n";
	else
		general_solver_information += " no\n";
	general_solver_information += "________________________________\n";
}
