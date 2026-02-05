#ifndef _FUNCTION_
#define _FUNCTION_

#include <iostream>
#include <stdexcept>
#include <cmath>
#include "linalg.h"

/**
 * @brief Template class for a vector.
 */
template <typename T>
class vec;

/**
 * @brief Template class for a matrix.
 */
template <typename T>
class Matrix;


// Function declarations
template <typename T>
vec<T> slice(vec<T> &v, int start, int end);

template <typename T>
Matrix<T> slice(Matrix<T> &M, int row_start, int row_end, int col_start, int col_end);

template <typename T>
Matrix<T> diagonalPreconditioner(Matrix<T> &M);

template <typename T>
Matrix<T> transpose(Matrix<T> &A);

template <typename T>
bool symmetrycheck(Matrix<T> &A);

template <typename T>
vec<T> flatten(Matrix<T> M);

template <typename T>
Matrix<T> unflatten(vec<T> v, int row, int col);

template <typename T>
Matrix<T> unflatten(vec<T> v);



// Slice a vector
template <typename T>
vec<T> slice(vec<T> &v, int start, int end)
{
    if (start < 0 || start >= v.size || end > v.size || end <= start)
    {
        throw std::invalid_argument("Invalid slice indices");
    }

    int n = end - start; 
    vec<T> v_sliced(n);

    for (int i = 0; i < n; ++i)
    {
        v_sliced(i) = v(start + i);
    }

    return v_sliced;
}

// Slice a matrix
template <typename T>
Matrix<T> slice(Matrix<T> &M, int row_start, int row_end, int col_start, int col_end)
{
    if (row_start < 0 || row_start >= M.row || row_end > M.row || row_start >= row_end ||
        col_start < 0 || col_start >= M.col || col_end > M.col || col_start >= col_end)
    {
        throw std::invalid_argument("Invalid slice indices");
    }

    int num_rows = row_end - row_start;
    int num_cols = col_end - col_start;

    Matrix<T> M_sliced(num_rows, num_cols);

    for (int i = 0; i < num_rows; ++i)
    {
        for (int j = 0; j < num_cols; ++j)
        {
            M_sliced(i, j) = M(row_start + i, col_start + j);
        }
    }

    return M_sliced;
}

// Create a diagonal preconditioner matrix
template <typename T>
Matrix<T> diagonalPreconditioner(Matrix<T> &M)
{
    Matrix<T> DPC(M.row, M.col);
    for (int i = 0; i < M.row; ++i)
    {
        DPC(i, i) = 1.0 / M(i, i);
    }
    return DPC;
}

// Transpose a matrix
template <typename T>
Matrix<T> transpose(Matrix<T> &A)
{
    Matrix<T> A_transpose(A.col, A.row);

    for (int i = 0; i < A.row; ++i)
    {
        for (int j = 0; j < A.col; ++j)
        {
            A_transpose(j, i) = A(i, j);
        }
    }

    return A_transpose;
}

// Check if a matrix is symmetric
template <typename T>
bool symmetrycheck(Matrix<T> &A)
{

    if (A.row != A.col)
    {
        std::cout << "Matrix is not symmetric (not a square matrix)." << std::endl;
        //return;
    }

    bool is_symmetric = true;
    for (int i = 0; i < A.row; ++i)
    {
        for (int j = 0; j < A.col; ++j)
        {
            if (A(i, j) != A(j, i))
            {
                is_symmetric = false;
                break;
            }
        }
        if (!is_symmetric)
            break;
    }
    
    return is_symmetric;
}

// Flatten a matrix to a vector
template <typename T>
vec<T> flatten(Matrix<T> M)
{
    int row = M.row;
    int col = M.col;
    vec<T> v(row * col);
    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < col; ++j)
        {
            v(i * col + j) = M(i, j);
        }
    }
    return v;
}

// Unflatten a vector to a matrix with specified rows and columns
template <typename T>
Matrix<T> unflatten(vec<T> v, int row, int col)
{
    if (row * col != v.size)
    {
        throw std::invalid_argument("Vector size does not match matrix dimensions");
    }

    Matrix<T> M(row, col);

    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < col; ++j)
        {
            M(i, j) = v(i * col + j);
        }
    }
    return M;
}

// Unflatten a vector to a square matrix
template <typename T>
Matrix<T> unflatten(vec<T> v)
{
    int n = static_cast<int>(std::sqrt(v.size));
    if (n * n != v.size)
    {
        throw std::invalid_argument("Vector size is not a perfect square");
    }

    Matrix<T> M(n, n);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            M(i, j) = v(i * n + j);
        }
    }
    return M;
}

#endif // _FUNCTION_
