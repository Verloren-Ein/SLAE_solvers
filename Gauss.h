#pragma once
#include "common.h"
#include "Matrix.h"
#include "Vector.h"

bool sol_gauss(double **A, double *b, double *x, int n);
bool direct_st(double **A, double *x, int n);
void transform(double **A, double *x, int i, int n);
void exchange(double **A, double *x, int first, int second, int n);

bool sol_gauss(double **A, vector<double> b, vector<double> &x, int n);
bool direct_st(double **A, vector<double> &x, int n);
void transform(double **A, vector<double> &x, int i, int n);
void exchange(double **A, vector<double> &x, int first, int second, int n);

bool sol_gauss(Matrix A, Vector b, Vector &x);
bool direct_st(Matrix &A, Vector &x, int n);
void transform(Matrix &A, Vector &x, int i, int n);
void exchange(Matrix &A, Vector &x, int first, int second, int n);