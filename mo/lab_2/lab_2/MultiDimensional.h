#pragma once
#ifndef MULTYDIMENSIONAL_H
#define MULTYDIMENSIONAL_H

#include"Assets.h"
#include"Utils.h"
#include"MatrixOperations.h"
#include"OneDimensional.h"

template<typename T>
solution gauss_zeidel_(point<T> p, Ftype(*f)(point<T>))
{
	size_t iterations = 0, i = 0;
	bool min_flag = false;
	point<Ftype> a1, a2, p_new = p;
	Ftype y1, y2, z0, z1;
	for (; iterations < constants::max_iterations && !min_flag; ++iterations)
	{
		p = p_new;																// сдвиг точки после спуска по всем координатам

		for (size_t i = 0; i < p.dimensionality(); ++i)
		{
			a1 = a2 = p_new;
			a1.updateCoord(i, (-1) * constants::eps, SHIFT);					// поиск направления в котором будет минимализироваться функция
			a2.updateCoord(i, constants::eps, SHIFT);							// 
#ifdef DEBUG
			std::cout << "a1: " << a1.to_string() << '\n';
			std::cout << "a2: " << a2.to_string() << '\n';
#endif
			y1 = f(a1);															// смотрим в какую сторону растёт функция
			y2 = f(a2);															//
			a1.updateCoord(i, (-1) * 100, SHIFT);								// расширяем область поиска
			a2.updateCoord(i, 100, SHIFT);										//
			a1 = y1 > y2 ? a2 : a1;												// выбираем направление поиска
			p_new = one_dimensional::fibonachi_kifer(a1, p_new, i, f).res;		// одномерная оптимизация по направлению i-ой координаты
			//std::cout << p_new.to_string() << '\n';
#ifdef DEBUG
			getchar();
#endif 
		}
		if (point<Ftype>::distance(p, p_new) < 2 * constants::eps)				// проверка на выход из цикла вычислений
		{
			min_flag = true;
		}
	}
	return solution(iterations, p, point<Ftype>::distance(p, p_new));
}


/// <summary>
/// Метод Хука-Дживса
/// </summary>
/// <returns></returns>
template<typename T>
solution hooke_jeevs(point<T> p, Ftype(*f)(point<T>))
{
	point<T> p_basis = p, p_new = p, p_searched = p;
	std::vector<T> h(p.dimensionality(), .3);
	bool end_flag = false;
	bool is_search_successfull = false;
	size_t iterations = 0;
	for (size_t i = 0; i < constants::max_iterations && !end_flag; ++i, ++iterations)
	{
#ifdef DEBUG
		if (iterations)
		{
			std::cout << "Iteration: " << std::setw(10) << iterations << '\n';
			std::cout << "p_basis: " << std::setw(10) << p_basis.to_string() << " f(p_b): " << std::setw(5) << f(p_basis)
				<< " | p_new " << std::setw(10) << p_new.to_string() << " f(p_n): " << std::setw(5) << f(p_new)
				<< " | steps (" << std::setprecision(4);
			for (size_t j = 0; j < p.dimensionality(); ++j)
			{
				std::cout << h[j];
				if (j < p.dimensionality() - 1)
				{
					std::cout << ", ";
				}
			}
			std::cout << ")\n";
			std::cout << std::string(50, '-') << '\n';
		}
#endif
		is_search_successfull = false;
		p_new = hooke_jeevs_utils::searchStep(p_basis, h, f, is_search_successfull);
		if (!is_search_successfull)
		{
			end_flag = hooke_jeevs_utils::updateStep(h);
		}
		else
		{
			p_basis = hooke_jeevs_utils::diagStep(p_new, p_basis);
		}
	}
	return solution(iterations, p_new, point<Ftype>::distance(p, p_new) / 2);
}


/// <summary>
/// Метод наискорейшего спуска
/// </summary>
/// <returns></returns>
template<typename T>
solution gradient_descent_(point<T> p, Ftype(*f)(point<T>))
{
	size_t iterations = 0;
	Ftype learning_rate = 2E-03;
	point<T> p_step = p;
	std::vector<Ftype> grad(p.dimensionality());

	bool end_flag = false;
	auto mult = [](Ftype lr, std::vector<Ftype> gr) {for (size_t i = 0; i < gr.size(); ++i) gr[i] *= lr; return gr; }; // умножение вектора на скаляр

	for (; iterations < constants::max_iterations && !end_flag; ++iterations)
	{
		grad = derivative_utils::gradient(p, f);											// вычисление градиента функции в точке p
		derivative_utils::norm(grad);														// нормализация градиента
		auto fun = [f, p, grad, mult](Ftype alpha) { return f(p - mult(alpha, grad)); };	// оптимизируемая функция по параметру alpha

		Ftype learning_rate = one_dimensional::dihotomia_method(fun, 0, 10);				// оптимизация параметра alpha через метод дихотомии
		p_step = p - mult(learning_rate, grad);												// шаг градиентного спуска

#ifdef DEBUG
		std::cout << "p_beg: " << p.to_string() << '\n';
		std::cout << "Learning_rate: " << learning_rate << '\n';
		std::cout << "p_step: " << p_step.to_string() << '\n';
		getchar();
#endif 

		if (point<Ftype>::distance(p, p_step) < 2 * constants::eps)							// проверка на выход из цикла вычислений
		{
			end_flag = true;
		}
		else
		{
			p = p_step;
		}
	}
	return solution(iterations, (p + p_step) / 2, point<Ftype>::distance(p, p_step) / 2);
}


template<typename T>
solution newton_method(point<T> p, Ftype(*f)(point<T>))
{
	size_t iterations = 0;
	bool counting_flag = true;
	point<T> p_i = p, p_i_1 = p;
	std::vector<std::vector<T>> hessian(p.dimensionality(), std::vector<T>(p.dimensionality(), 1));
	std::vector<std::vector<T>> inv_hessian = hessian;
	for (; iterations < constants::max_iterations && counting_flag; ++iterations)
	{
		p_i = p_i_1;
		std::vector<T> grad = derivative_utils::gradient(p_i, f);			// считаем градиент в точке p_i
		hessian = second_order_utils::hessian(f, p_i);						// считаем матрицу Гессе в точке p_i
		inv_hessian = matrix_operations::inverse_(hessian, hessian.size());	// считаем обратную матрицу Гессе
		point<T> mult = second_order_utils::multiply_(inv_hessian, grad);	// умножение градиента на матрицу Гессе
		p_i_1 = p_i - mult;													// делаем шаг

#ifdef DEBUG
		if (iterations >= 0)
		{
			std::cout << "Iteration " << iterations << '\n';
			std::cout << "p " << p_i.to_string() << " | grad (";
			for (size_t j = 0; j < p_i.dimensionality(); ++j)
			{
				std::cout << grad[j];
				if (j < p.dimensionality() - 1)
				{
					std::cout << ", ";
				}
			}
			std::cout << ") " << '\n';
			std::cout << "Hessian" << '\n';
			matrix_operations::printMatrix(hessian);
			std::cout << "Inv_Hessian" << '\n';
			matrix_operations::printMatrix(inv_hessian);
			std::cout << "Mult | " << mult.to_string() << "\n";
			std::cout << "p_i_1 | " << p_i_1.to_string() << '\n';
			std::cout << std::string(50, '-') << '\n';
			getchar();
		}
#endif

		if (point<Ftype>().distance(p_i, p_i_1) < 2 * constants::eps)		// проверка на выход из цикла вычислений
		{
			counting_flag = false;
		}
	}
	return solution(iterations, p_i, point<Ftype>::distance(p_i, p_i_1) / 2);
}

#endif // !MULTYDIMENSIONAL_H