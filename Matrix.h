#pragma once
#include "common.h"
//using primitive multiply algorithm

template<typename T>
class Matrix
{
private:
	size_t rows_, cols_;
	vector<T> elements_;

public:
	Matrix()
		:rows_()
		,cols_()
	{
		elements_.clear();
	}
	Matrix(size_t n)
		:rows_(n)
		, cols_(n)
	{
		elements_.clear();
		elements_.resize(n*n);
	}
	Matrix(size_t rows, size_t cols)
		:rows_(rows)
		,cols_(cols)
	{
		elements_.clear();
		elements_.resize(rows_*cols_);
	}
	~Matrix()
	{
		elements_.clear();
		elements_.~vector();
	}
	void clear();
	void reset();
	void display();
	void display(streamsize precision);
	void print(ofstream& ostr);
	void print(ofstream& ostr, streamsize precision);
	size_t rows();
	size_t cols();
	/*returns transposed matrix*/
	Matrix<T> transpose();
	/*returns inversed matrix*/
	Matrix<T> inverse();
	/*returns matrix determinant*/
	T determinant();
	/*strike out given row and col and returns new matrix*/
	Matrix<T> strike_out(size_t row, size_t col);

	/*sets current matrix with elements*/
	void set_matrix(T** elements);
	/*sets current matrix with element for all matrix elements*/
	void set_matrix(T element);
	/*sets current matrix to identity matrix*/
	void set_identity_matrix();
	/*sets current matrix to diagonal matrix with diagonal_element for all diagonal elements*/
	void set_diagonal_matrix(T diagonal_element);
	/*sets current matrix to diagonal matrix with diagonal_element for all diagonal elements*/
	/*diagonal_number values: 0 - main diagonal, >0 upper, <0 lower*/
	void set_subdiagonal_matrix(int diagonal_number, T diagonal_element);
	/*sets current matrix to diagonal matrix with diagonal_elements for diagonal*/
	void set_diagonal_matrix(vector<T> diagonal_elements);
	/*sets current matrix to diagonal matrix with diagonal_elements for diagonal*/
	/*diagonal_number values: 0 - main diagonal, >0 upper, <0 lower*/
	void set_subdiagonal_matrix(int diagonal_number, vector<T> diagonal_elements);
	/*returns vector contains diagonal_number elements*/
	/*diagonal_number values: 0 - main diagonal, >0 upper, <0 lower*/
	vector<T> get_diagonal(int diagonal_number);
	/*returns main diagonal elements*/
	vector<T> get_main_diagonal();
	/*returns vector contains row_number elements*/
	vector<T> get_row(size_t row_number);
	/*returns vector contains col_number elements*/
	vector<T> get_col(size_t col_number);


	/*operators overloading*/
	/*subscript operator*/
	T& operator() (size_t row, size_t col);
	T operator() (size_t row, size_t col) const;

	/*unary operators*/
	friend const Matrix<T>& operator+(const Matrix<T>& i)
	{
		return i;
	}
	friend const Matrix<T> operator-(const Matrix<T>& i)
	{
		Matrix<T> matr = Matrix<T>(i.rows_, i.cols_);
		for(size_t k = 0; k < i.elements_.size(); k++)
		{
			matr.elements_[k] = -i.elements_[k];
		}
		return matr;
	}

	friend Matrix<T>& operator+=(Matrix<T>& left, const Matrix<T>& right)
	{
		left.plus(right);
		return left;
	}
	friend Matrix<T>& operator-=(Matrix<T>& left, const Matrix<T>& right)
	{
		left.minus(right);
		return left;
	}
	friend Matrix<T>& operator*=(Matrix<T>& left, const Matrix<T>& right)
	{
		left.multiply(right);
		return left;
	}

	/*binary operators*/
	friend bool operator ==(Matrix<T> &A, Matrix<T> &B)
	{
		try
		{
			if(A.cols_ == B.cols_ && A.rows_ == B.rows_)
			{
				double precision = get_current_compare_precision();
				for(size_t i = 0; i < A.elements_.size(); i++)
				{
					if(abs(A.elements_[i] - B.elements_[i]) > precision)
					{
						return false;
					}
					else continue;
				}
				return true;
			}
			else
				throw "Error: attempt to compare matrix of different sizes";
		}
		catch(char *str)
		{
			cerr << str << endl;
		}
		return false;
	}
	friend bool operator !=(Matrix<T> &A, Matrix<T> &B)
	{
		try
		{
			if(A.cols_ == B.cols_ && A.rows_ == B.rows_)
			{
				double precision = get_current_compare_precision();
				for(size_t i = 0; i < A.elements_.size(); i++)
				{
					if(abs(A.elements_[i] - B.elements_[i]) > precision)
					{
						return true;
					}
					else continue;
				}
				return false;
			}
			else
				throw "Error: attempt to compare matrix of different sizes";
		}
		catch(char *str)
		{
			cerr << str << endl;
		}
		return false;
	}

	friend const Matrix<T> operator+(const Matrix<T>& left, const Matrix<T>& right)
	{
		Matrix<T> result = left;
		result.plus(right);
		return result;
	}
	friend const Matrix<T> operator-(const Matrix<T>& left, const Matrix<T>& right)
	{
		Matrix<T> result = left;
		result.minus(right);
		return result;
	}
	friend const Matrix<T> operator*(const Matrix<T>& left, const Matrix<T>& right)
	{
		Matrix<T> result = left;
		result.multiply(right);
		return result;
	}

protected:
	/*arithmetic operations*/
	virtual void plus(Matrix<T> arg);		// elements[i] + arg[i]
	virtual void minus(Matrix<T> arg);		// elements[i] - arg[i]
	virtual void multiply(Matrix<T> arg);	// elements[i] * arg[i]

private:
	/*sets current matrix with vector elements for all matrix elements : internal use only*/
	void set_matrix(vector<T> elements);
	/*recursive determinant calculation : internal use only*/
	double recursive_determinant(Matrix<T> matr);
	/*strike out given row and col and returns new matrix : internal use only*/
	Matrix<T> strike_out(Matrix<T> matrix, size_t row, size_t col);
};

template<typename T>
inline void Matrix<T>::clear()
{
	rows_ = size_t();
	cols_ = size_t();
	elements_.clear();
}

template<typename T>
inline void Matrix<T>::reset()
{
	set_matrix(T(0));
}

template<typename T>
inline void Matrix<T>::display()
{
	cout << fixed;
	cout.precision(CURRENT_DISPLAY_PRECISION);

	for(size_t i = 0; i < rows_; i++)
	{
		for(size_t j = 0; j < cols_; j++)
		{
			cout << this->operator()(i, j) << " ";
		}
		cout << endl;
	}
	cout << endl;
}

template<typename T>
inline void Matrix<T>::display(streamsize precision)
{
	cout << fixed;
	cout.precision(precision);

	for(size_t i = 0; i < rows_; i++)
	{
		for(size_t j = 0; j < cols_; j++)
		{
			cout << this->operator()(i, j) << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(CURRENT_DISPLAY_PRECISION);
}

template<typename T>
inline void Matrix<T>::print(ofstream & ostr)
{
	ostr << fixed;
	ostr.precision(CURRENT_DISPLAY_PRECISION);

	for(size_t i = 0; i < rows_; i++)
	{
		for(size_t j = 0; j < cols_; j++)
		{
			ostr << this->operator()(i, j) << " ";
		}
		ostr << endl;
	}
	ostr << endl;
}

template<typename T>
inline void Matrix<T>::print(ofstream & ostr, streamsize precision)
{
	ostr << fixed;
	ostr.precision(precision);

	for(size_t i = 0; i < rows_; i++)
	{
		for(size_t j = 0; j < cols_; j++)
		{
			ostr << this->operator()(i, j) << " ";
		}
		ostr << endl;
	}
	ostr << endl;
	ostr.precision(CURRENT_DISPLAY_PRECISION);
}

template<typename T>
inline size_t Matrix<T>::rows()
{
	return rows_;
}

template<typename T>
inline size_t Matrix<T>::cols()
{
	return cols_;
}

template<typename T>
inline Matrix<T> Matrix<T>::transpose()
{
	Matrix<T> transposed_matrix = Matrix<T>(this->cols(), this->rows());
	for(size_t i = 0; i < rows_; i++)
	{
		for(size_t j = 0; j < cols_; j++)
		{
			transposed_matrix(j, i) = this->operator()(i, j);
		}
	}

	return transposed_matrix;
}

template<typename T>
inline Matrix<T> Matrix<T>::inverse()
{
	Matrix<T> inverse_matrix;
	Matrix<T> new_matrix;

	/*calculate determinant*/
	T det = determinant();

	try
	{
		/*if determinant is zero, matrix can not be inversed*/
		if(abs(det) < almost_zero)
			throw "Error: matrix can not be inversed";
		else
		{
			inverse_matrix = Matrix<T>(rows_, cols_);
			for(size_t i = 0; i < rows_; i++)
			{
				for(size_t j = 0; j < cols_; j++)
				{
					new_matrix = strike_out(i, j);
					inverse_matrix(i, j) = pow(-1.0, i + j + 2)*recursive_determinant(new_matrix) / det;
					new_matrix.~Matrix();
				}
			}
			inverse_matrix = inverse_matrix.transpose();
		}
	}
	catch(char *str)
	{
		cerr << str << endl;
	}

	return inverse_matrix;
}

template<typename T>
inline T Matrix<T>::determinant()
{
	/*copy base matrix*/
	Matrix<T> new_matrix = Matrix<T>(this->rows_, this->cols_);
	new_matrix.set_matrix(this->elements_);
	T determinant = 0;

	try
	{
		if(rows_ != cols_)
			throw "Error: attempt to calculate a determinant of non-square matrix";
		else
		{
			determinant = recursive_determinant(new_matrix);
		}
	}
	catch(char *str)
	{
		cerr << str << endl;
	}

	new_matrix.~Matrix();
	return determinant;
}

template<typename T>
inline Matrix<T> Matrix<T>::strike_out(size_t row, size_t col)
{
	Matrix<T> new_matrix;
	try
	{
		if(row >= this->rows() || col >= this->cols())
			throw "Error: matrix subscript out of range";
		else
		{
			new_matrix = Matrix<T>(this->rows() - 1, this->cols() - 1);
			for(size_t i = 0, new_i = 0; i < this->rows(); i++)
			{
				if(i != row)
				{
					for(size_t j = 0, new_j = 0; j < this->cols(); j++)
					{
						if(j != col)
						{
							new_matrix(new_i, new_j) = this->operator()(i, j);
							new_j++;
						}
					}
					new_i++;
				}
			}
		}
	}
	catch(char *str)
	{
		cerr << str << endl;
	}
	return new_matrix;
}

template<typename T>
inline void Matrix<T>::set_matrix(T ** elements)
{
	for(size_t i = 0; i < rows_; i++)
	{
		for(size_t j = 0; j < cols_; j++)
		{
			this->operator()(i, j) = elements[i][j];
		}
	}
}

template<typename T>
inline void Matrix<T>::set_matrix(T element)
{
	for(size_t i = 0; i < rows_; i++)
	{
		for(size_t j = 0; j < cols_; j++)
		{
			this->operator()(i, j) = element;
		}
	}
}

template<typename T>
inline void Matrix<T>::set_identity_matrix()
{
	set_diagonal_matrix(T(1.));
}

template<typename T>
inline void Matrix<T>::set_diagonal_matrix(T diagonal_element)
{
	for(size_t i = 0; i < rows_; i++)
	{
		for(size_t j = 0; j < cols_; j++)
		{
			if(i == j)
				this->operator()(i, j) = diagonal_element;
			else
				this->operator()(i, j) = 0.;
		}
	}
}

template<typename T>
inline void Matrix<T>::set_subdiagonal_matrix(int diagonal_number, T diagonal_element)
{
	try
	{
		if(diagonal_number == 0)
			set_diagonal_matrix(diagonal_element);
		if(diagonal_number > 0)
		{
			if(diagonal_number > int(cols_ - 1))
				throw "Error: diagonal number in matrix out of range";
			else
			{
				reset();
				for(size_t i = 0, j = diagonal_number; i < rows_ && j < cols_; i++, j++)
				{
					this->operator()(i, j) = diagonal_element;
				}
			}
		}
		if(diagonal_number < 0)
		{
			if(abs(diagonal_number) > int(rows_ - 1))
				throw "Error: diagonal number in matrix out of range";
			else
			{
				reset();
				for(size_t i = 0, j = abs(diagonal_number); i < cols_ && j < rows_; i++, j++)
				{
					this->operator()(j, i) = diagonal_element;
				}
			}
		}
	}
	catch(char *str)
	{
		cerr << str << endl;
	}
}

template<typename T>
inline void Matrix<T>::set_diagonal_matrix(vector<T> diagonal_elements)
{
	for(size_t i = 0; i < rows_; i++)
	{
		for(size_t j = 0; j < cols_; j++)
		{
			if(i == j)
				this->operator()(i, j) = diagonal_elements[i];
			else
				this->operator()(i, j) = 0.;
		}
	}
}

template<typename T>
inline void Matrix<T>::set_subdiagonal_matrix(int diagonal_number, vector<T> diagonal_elements)
{
	//no check vector size
	try
	{
		if(diagonal_number == 0)
			set_diagonal_matrix(diagonal_elements);
		if(diagonal_number > 0)
		{
			if(diagonal_number > int(cols_ - 1))
				throw "Error: diagonal number in matrix out of range";
			else
			{
				reset();
				for(size_t i = 0, j = diagonal_number; i < rows_ && j < cols_; i++, j++)
				{
					this->operator()(i, j) = diagonal_elements[i];
				}
			}
		}
		if(diagonal_number < 0)
		{
			if(abs(diagonal_number) > int(rows_ - 1))
				throw "Error: diagonal number in matrix out of range";
			else
			{
				reset();
				for(size_t i = 0, j = abs(diagonal_number); i < cols_ && j < rows_; i++, j++)
				{
					this->operator()(j, i) = diagonal_elements[i];
				}
			}
		}
	}
	catch(char *str)
	{
		cerr << str << endl;
	}
}

template<typename T>
inline vector<T> Matrix<T>::get_diagonal(int diagonal_number)
{
	vector<T> diagonal;
	try
	{
		if(diagonal_number == 0)
			diagonal = get_main_diagonal();
		if(diagonal_number > 0)
		{
			if(diagonal_number > int(cols_ - 1))
				throw "Error: diagonal number in matrix out of range";
			else
			{
				for(size_t i = 0, j = diagonal_number; i < rows_ && j < cols_; i++, j++)
				{
					diagonal.push_back(this->operator()(i, j));
				}
			}
		}
		if(diagonal_number < 0)
		{
			if(abs(diagonal_number) > int(rows_ - 1))
				throw "Error: diagonal number in matrix out of range";
			else
			{
				for(size_t i = 0, j = abs(diagonal_number); i < cols_ && j < rows_; i++, j++)
				{
					diagonal.push_back(this->operator()(j, i));
				}
			}
		}
	}
	catch(char *str)
	{
		cerr << str << endl;
	}
	return diagonal;
}

template<typename T>
inline vector<T> Matrix<T>::get_main_diagonal()
{
	vector<T> diagonal;
	for(size_t i = 0; i < rows_; i++)
	{
		diagonal.push_back(this->operator()(i, i));
	}
	return diagonal;
}

template<typename T>
inline vector<T> Matrix<T>::get_row(size_t row_number)
{
	vector<T> row;
	try
	{
		if(row_number >= rows_)
			throw "Error: matrix subscript out of range";

		else
		{
			for(size_t i = 0; i < cols_; i++)
			{
				row.push_back(elements_[cols_*row_number + i]);
			}
		}
	}
	catch(char *str)
	{
		cerr << str << endl;
	}
	return row;
}

template<typename T>
inline vector<T> Matrix<T>::get_col(size_t col_number)
{
	vector<T> col;
	try
	{
		if(col_number >= cols_)
			throw "Error: matrix subscript out of range";

		else
		{
			for(size_t i = 0; i < rows_; i++)
			{
				col.push_back(elements_[col_number + i * cols_]);
			}
		}
	}
	catch(char *str)
	{
		cerr << str << endl;
	}
	return col;
}

template<typename T>
inline T & Matrix<T>::operator()(size_t row, size_t col)
{
	try
	{
		if(row >= rows_ || col >= cols_)
			throw "Error: matrix subscript out of range";

		else
			return elements_[cols_*row + col];
	}
	catch(char *str)
	{
		cerr << str << endl;
	}
}

template<typename T>
inline T Matrix<T>::operator()(size_t row, size_t col) const
{
	try
	{
		if(row >= rows_ || col >= cols_)
			throw "Error: matrix subscript out of range";

		else
			return elements_[cols_*row + col];
	}
	catch(char *str)
	{
		cerr << str << endl;
	}
}

template<typename T>
inline void Matrix<T>::plus(Matrix<T> arg)
{
	try
	{
		if(this->cols_ != arg.cols_ || this->rows_ != arg.rows_)
			throw "Error: matrix cannot be summarize";
		else
		{
			for(size_t i = 0; i < this->elements_.size(); i++)
			{
				elements_[i] += arg.elements_[i];
			}
		}
	}
	catch(char *str)
	{
		cerr << str << endl;
	}
}

template<typename T>
inline void Matrix<T>::minus(Matrix<T> arg)
{
	try
	{
		if(this->cols_ != arg.cols_ || this->rows_ != arg.rows_)
			throw "Error: matrix cannot be subtracted";
		else
		{
			for(size_t i = 0; i < this->elements_.size(); i++)
			{
				elements_[i] -= arg.elements_[i];
			}
		}
	}
	catch(char *str)
	{
		cerr << str << endl;
	}
}

template<typename T>
inline void Matrix<T>::multiply(Matrix<T> arg)
{
	Matrix matr = Matrix(rows_, arg.cols_);
	try
	{
		if(this->cols_ != arg.rows_)
			throw "Error: matrix cannot be multiplied";
		else
		{
			for(size_t row = 0; row < this->rows_; row++)
			{
				for(size_t col = 0; col < arg.cols_; col++)
				{
					for(size_t inner = 0; inner < this->cols_; inner++)
					{
						matr(row, col) += this->operator()(row, inner) * arg(inner, col);
					}
				}
			}
			this->elements_ = matr.elements_;
			this->rows_ = matr.rows_;
			this->cols_ = matr.cols_;
			matr.~Matrix();
		}
	}
	catch(char *str)
	{
		cerr << str << endl;
	}
}

template<typename T>
inline void Matrix<T>::set_matrix(vector<T> elements)
{
	for(size_t i = 0; i < elements_.size(); i++)
	{
		elements_[i] = elements[i];
	}
}

template<typename T>
inline double Matrix<T>::recursive_determinant(Matrix<T> matr)
{
	T current_determinant = 0;
	int dim = matr.cols();
	int k = 1; //power

	if(dim == 1)
	{
		current_determinant = matr(0, 0);
	}
	else
	{
		/*simple case handling*/
		if(dim == 2)
		{
			current_determinant = matr(0, 0)*matr(1, 1) - matr(1, 0)*matr(0, 1);
		}
		else
		{
			for(size_t i = 0; i < dim; i++)
			{
				Matrix<T> new_matr = strike_out(matr, 0, i);
				current_determinant = current_determinant + k * matr(0, i)*recursive_determinant(new_matr);
				k = -k;

				new_matr.~Matrix();
			}
		}
	}
	return current_determinant;
}

template<typename T>
inline Matrix<T> Matrix<T>::strike_out(Matrix<T> matrix, size_t row, size_t col)
{
	Matrix<T> new_matrix;
	try
	{
		if(row >= matrix.rows() || col >= matrix.cols())
			throw "Error: matrix subscript out of range";
		else
		{
			new_matrix = Matrix<T>(matrix.rows() - 1, matrix.cols() - 1);
			for(size_t i = 0, new_i = 0; i < matrix.rows(); i++)
			{
				if(i != row)
				{
					for(unsigned j = 0, new_j = 0; j < matrix.cols(); j++)
					{
						if(j != col)
						{
							new_matrix(new_i, new_j) = matrix(i, j);
							new_j++;
						}
					}
					new_i++;
				}
			}
		}
	}
	catch(char *str)
	{
		cerr << str << endl;
	}
	return new_matrix;
}
