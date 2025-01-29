#pragma once
#ifndef MATRIXOPERATIONS_H
#define MATRIXOPERATIONS_H
#include"Assets.h"

/// <summary>
/// пространство имён с функциями обращения и нахождения определителя матрицы
/// </summary>
namespace matrix_operations
{
	/// <summary>
	/// Операция нахождения определителя матрицы
	/// </summary>
	/// <typeparam name="T"></typeparam>
	/// <param name="matr">- матрица</param>
	/// <param name="n">- ранг матрицы</param>
	/// <returns></returns>
	template<typename T>
	Ftype determinant(std::vector<std::vector<T>>& matr, size_t n)
	{
		if (n == 1) {
			return matr[0][0];
		}
		else if (n == 2) {
			return matr[0][0] * matr[1][1] - matr[0][1] * matr[1][0];
		}
		else {
			Ftype det = 0;
			for (size_t j = 0; j < n; j++) {
				std::vector<std::vector<T>> submatrix(n - 1, std::vector<T>(n - 1));
				size_t row = 0;
				for (size_t i = 1; i < n; i++) {
					size_t col = 0;
					for (size_t k = 0; k < n; k++) {
						if (k != j) {
							submatrix[row][col] = matr[i][k];
							col++;
						}
					}
					row++;
				}
				det += pow(-1, j) * matr[0][j] * determinant(submatrix, n - 1);
			}
			return det;
		}
	}

	template<typename T>
	std::vector<std::vector<T>> cofactorMatrix(std::vector<std::vector<T>>& matrix, size_t n) {
		std::vector<std::vector<T>> cofactor(n, std::vector<T>(n));
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < n; j++) {
				std::vector<std::vector<T>> submatrix(n - 1, std::vector<T>(n - 1));
				size_t row = 0;
				for (size_t k = 0; k < n; k++) {
					if (k != i) {
						size_t col = 0;
						for (size_t l = 0; l < n; l++) {
							if (l != j) {
								submatrix[row][col] = matrix[k][l];
								col++;
							}
						}
						row++;
					}
				}
				cofactor[i][j] = pow(-1, i + j) * matrix_operations::determinant(submatrix, n - 1);
			}
		}
		return cofactor;
	}

	/// <summary>
	/// Операция обращения матрицы
	/// </summary>
	/// <typeparam name="T"></typeparam>
	/// <param name="matrix"> - матрица</param>
	/// <param name="n"> - ранг матрицы</param>
	/// <returns></returns>
	template<typename T>
	std::vector<std::vector<T>> inverse(std::vector<std::vector<T>>& matrix, size_t n) {
		Ftype det = determinant(matrix, n);
		if (det == 0) {
			std::cout << "Matrix is singular. Inverse does not exist." << '\n';
			return matrix;
		}
		else {
			std::vector<std::vector<T>> adjoint(n, std::vector<T>(n));
			adjoint = cofactorMatrix(matrix, n);
			std::vector<std::vector<T>> inverse(n, std::vector<T>(n));
			for (size_t i = 0; i < n; i++) {
				for (size_t j = 0; j < n; j++) {
					inverse[i][j] = adjoint[j][i] / det;
				}
			}
			return inverse;
		}
	}

	template<typename T>
	std::vector<std::vector<T>> inverse_(std::vector<std::vector<T>>& matrix, size_t n) {
		std::vector<std::vector<T>> inverse(matrix);
		for (size_t i = 0; i < n; i++)
		{
			inverse[i][i] = 1 / inverse[i][i];

		}
		return inverse;
	}

	/// <summary>
	/// Вывод матрицы
	/// </summary>
	void printMatrix(const std::vector<std::vector<Ftype>>& matrix) {
		size_t rows = matrix.size();
		size_t cols = matrix[0].size();
		std::cout << '(' << '\n';
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				std::cout << std::setw(5) << std::setprecision(3) << matrix[i][j] << " ";
			}
			std::cout << '\n';
		}
		std::cout << ')' << '\n';
	}
}
#endif // !MATRIXOPERATIONS_H