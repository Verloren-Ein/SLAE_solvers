#include"Gauss.h"

bool sol_gauss(double **A, double *b, double *x, int n)
{
	bool er;
	int i, j;
	double s;
	er = direct_st(A, b, n);
	if (!er) return er;
	for (i = n - 1; i >= 0; i--)
	{
		s = 0;
		for (j = n - 1; j > i; j--)
			s += A[i][j] * x[j];
		x[i] = (b[i] - s) / A[i][i];
	}
	return true;
}
bool direct_st(double **A, double *x, int n)
{
	int i, j, k;
	double koeff;

	for (i = 0; i < n; i++)
	{
		transform(A, x, i, n);
		if (fabs(A[i][i]) < 1E-20) return false;

		for (j = i + 1; j < n; j++)
		{
			koeff = -A[j][i] / A[i][i];
			for (k = i; k < n; k++)
				A[j][k] += A[i][k] * koeff;
			x[j] += x[i] * koeff;
		}
	}
	return true;
}
void transform(double **A, double *x, int i, int n)
{
	double max = A[i][i];
	int gl = i;
	int j;
	for (j = i + 1; j < n; j++)
		if (fabs(A[j][i]) > fabs(max))
		{
			max = A[j][i];
			gl = j;
		}
	if (gl != i) exchange(A, x, i, gl, n);
}
void exchange(double **A, double *x, int first, int second, int n)
{
	int i;
	double mid;
	mid = x[first];
	x[first] = x[second];
	x[second] = mid;
	for (i = 0; i < n; i++)
	{
		mid = A[first][i];
		A[first][i] = A[second][i];
		A[second][i] = mid;
	}
}

bool sol_gauss(double ** A, vector<double> b, vector<double> &x, int n)
{
	bool er;
	int i, j;
	double s;
	er = direct_st(A, b, n);
	if (!er) return er;
	for (i = n - 1; i >= 0; i--)
	{
		s = 0;
		for (j = n - 1; j > i; j--)
			s += A[i][j] * x[j];
		x[i] = (b[i] - s) / A[i][i];
	}
	return true;
}
bool direct_st(double ** A, vector<double> &x, int n)
{
	int i, j, k;
	double koeff;

	for (i = 0; i < n; i++)
	{
		transform(A, x, i, n);
		if (fabs(A[i][i]) < 1E-20) return false;

		for (j = i + 1; j < n; j++)
		{
			koeff = -A[j][i] / A[i][i];
			for (k = i; k < n; k++)
				A[j][k] += A[i][k] * koeff;
			x[j] += x[i] * koeff;
		}
	}
	return true;
}
void transform(double ** A, vector<double> &x, int i, int n)
{
	double max = A[i][i];
	int gl = i;
	int j;
	for (j = i + 1; j < n; j++)
		if (fabs(A[j][i]) > fabs(max))
		{
			max = A[j][i];
			gl = j;
		}
	if (gl != i) exchange(A, x, i, gl, n);
}
void exchange(double ** A, vector<double> &x, int first, int second, int n)
{
	int i;
	double mid;
	mid = x[first];
	x[first] = x[second];
	x[second] = mid;
	for (i = 0; i < n; i++)
	{
		mid = A[first][i];
		A[first][i] = A[second][i];
		A[second][i] = mid;
	}
}

bool sol_gauss(Matrix A, Vector b, Vector & x)
{
	bool er;
	int i, j;
	double s;
	int n = b.size();
	er = direct_st(A, b, n);
	if(!er) return er;
	for(i = n - 1; i >= 0; i--)
	{
		s = 0;
		for(j = n - 1; j > i; j--)
			s += A(i, j) * x[j];
		x[i] = (b[i] - s) / A(i, i);
	}
	return true;
}
bool direct_st(Matrix &A, Vector & x, int n)
{
	int i, j, k;
	double koeff;

	for(i = 0; i < n; i++)
	{
		transform(A, x, i, n);
		if(fabs(A(i, i)) < 1E-20) return false;

		for(j = i + 1; j < n; j++)
		{
			koeff = -A(j, i) / A(i, i);
			for(k = i; k < n; k++)
				A(j, k) += A(i, k) * koeff;
			x[j] += x[i] * koeff;
		}
	}
	return true;
}
void transform(Matrix &A, Vector & x, int i, int n)
{
	double max = A(i, i);
	int gl = i;
	int j;
	for(j = i + 1; j < n; j++)
		if(fabs(A(j, i)) > fabs(max))
		{
			max = A(j, i);
			gl = j;
		}
	if(gl != i) exchange(A, x, i, gl, n);
}
void exchange(Matrix &A, Vector & x, int first, int second, int n)
{
	int i;
	double mid;
	mid = x[first];
	x[first] = x[second];
	x[second] = mid;
	for(i = 0; i < n; i++)
	{
		mid = A(first, i);
		A(first, i) = A(second, i);
		A(second, i) = mid;
	}
}
