#pragma once
#ifndef UTILS_H
#define UTILS

#include"Assets.h"
#include"Point.h"
struct solution
{
	size_t iterations;
	point<Ftype> res;
	Ftype precission;

	solution()
	{
		this->iterations = 0;
		this->precission = {};
		this->res = point<Ftype>();
	}

	solution(size_t iterations, point<Ftype> res)
	{
		this->iterations = iterations;
		this->res = res;
	}

	solution(size_t iterations, point<Ftype> res, Ftype precission)
	{
		this->iterations = iterations;
		this->res = res;
		this->precission = precission;
	}
};

template<typename T>
solution gauss_zeidel_(point<T>, Ftype(*f)(point<T>));

template<typename T>
solution gradient_descent_(point<T>, Ftype(*f)(point<T>));

template<typename T>
solution newton_method(point<T>, Ftype(*f)(point<T>));

/// <summary>
/// Пространство имён с константами и утилитной функцией вычисления чисел Фибоначчи
/// </summary>
namespace constants
{
	/// <summary>
	/// Функция вычисления n-ого числа Фибоначчи
	/// </summary>
	/// <returns></returns>
	int F(size_t n)
	{
		Ftype phi = (1 + sqrt(5)) / 2;
		return round(pow(phi, n) / sqrt(5));
	}
	Ftype eps = 0.00001;
	size_t max_iterations = 1000;
	double pi = 3.14159;
	Ftype alpha = pi / 3;
	std::vector<Ftype> shifts = { 0, 1, -4 };
}

// y_ = (y - y_0) * cos(a) + (z - z_0) * sin(a)
template<typename T>
Ftype y_(point<T> p)
{
	return (p[1] - constants::shifts[1]) * cos(constants::alpha) + (p[2] - constants::shifts[2]) * sin(constants::alpha);
}

// z_ = (z - z_0) * cos(a) - (y - y_0) * sin(a)
template<typename T>
Ftype z_(point<T> p)
{
	return (p[2] - constants::shifts[2]) * cos(constants::alpha) - (p[1] - constants::shifts[1]) * sin(constants::alpha);
}

namespace hooke_jeevs_utils
{
	template<typename T>
	point<T> searchStep(point<T> x, std::vector<T> h, Ftype(*f)(point<T>), bool& flag)
	{
		T f_x = f(x);
		point<T> x_p, x_m;
		T x_i_p, x_i_m, f_x_m, f_x_p;
		for (size_t i = 0; i < x.dimensionality(); ++i)
		{
			x_i_p = x[i] + h[i];
			x_i_m = x[i] - h[i];
			x_p.setCoordinates(x, x_i_p, i);
			x_m.setCoordinates(x, x_i_m, i);
			f_x_p = f(x_p);
			f_x_m = f(x_m);
			if (f_x_m < f_x)
			{

				x = x_m;
				f_x = f_x_m;
				flag = true;
			}
			else if (f_x_p < f_x)
			{
				x = x_p;
				f_x = f_x_p;
				flag = true;
			}
		}
		return x;
	}

	template<typename T>
	bool updateStep(std::vector<T>& h)
	{
		if (h[0] / 2 < constants::eps)
		{
			return true;
		}
		for (size_t i = 0; i < h.size(); ++i)
		{
			h[i] /= 2;
		}
		return false;
	}

	template<typename T>
	point<T> diagStep(point<T> p_prev, point<T> p_cur)
	{
		point<T> p = p_cur;
		for (size_t i = 0; i < p.dimensionality(); ++i)
		{
			p.setCoord(i, 2 * p_prev[i] - p_cur[i]);
		}
		return p;
	}
}

namespace derivative_utils
{
	/// <summary>
	/// Подсчёт частной производной по i-ой переменной, i = ind
	/// </summary>
	/// <returns></returns>
	template<typename T>
	Ftype partial_derivative(point<T> p, Ftype(*f)(point<T>), size_t ind)
	{
		point<T> delta_p(p, p[ind] + 1.0E-7, ind);
		return (f(delta_p) - f(p)) / 1.0E-7;
	}

	/// <summary>
	/// Подсчёт градиента
	/// </summary>
	/// <returns></returns>
	template<typename T>
	std::vector<Ftype> gradient(point<T> p, Ftype(*f)(point<T>))
	{
		std::vector<Ftype> grad(p.dimensionality());
		for (size_t i = 0; i < p.dimensionality(); ++i)
		{
			grad[i] = derivative_utils::partial_derivative(p, f, i);
		}
		return grad;
	}

	void norm(std::vector<Ftype>& grad)
	{
		Ftype square_len = 0;
		for (size_t i = 0; i < grad.size(); ++i)
		{
			square_len += std::pow(grad[i], 2);
		}
		Ftype len = sqrt(square_len);
		for (size_t i = 0; i < grad.size(); ++i)
		{
			grad[i] /= len;
		}
	}

	/// <summary>
	/// Подсчёт второй частной производной по i-ой переменной, i = ind
	/// </summary>
	/// <returns></returns>
	template<typename T>
	Ftype second_partial_derivative(point<T> p, Ftype(*f)(point<T>), size_t ind1, size_t ind2)
	{
		point<T> delta_p(p, p[ind2] + 1.0E-5, ind2);
		return (partial_derivative(delta_p, f, ind1) - partial_derivative(p, f, ind1)) / 1.0E-5;
	}
}

namespace second_order_utils
{
	/// <summary>
	/// Фунция вычисления матрицы Гессе
	/// </summary>
	/// <typeparam name="T"></typeparam>
	/// <param name="f"> - функция, для которой вычисляется матрица</param>
	/// <param name="p"> - точка, в которой вычисляется матрица</param>
	/// <returns></returns>
	template<typename T>
	std::vector<std::vector<T>> hessian(Ftype(*f)(point<T>), point<T> p)
	{
		std::vector<std::vector<T>> matr(p.dimensionality(), std::vector<T>(p.dimensionality()));
		for (size_t i = 0; i < p.dimensionality(); ++i)
		{
			for (size_t j = 0; j <= i; ++j)
			{
				matr[i][j] = derivative_utils::second_partial_derivative(p, f, i, j);
				matr[j][i] = matr[i][j];
			}
		}
		return matr;
	}

	/// <summary>
	/// Функция перемножения матрицы на вектор
	/// </summary>
	/// <typeparam name="T"></typeparam>
	/// <param name="matr">- матрица</param>
	/// <param name="vec">- вектор</param>
	/// <returns></returns>
	template<typename T>
	point<T> multiply(std::vector<std::vector<T>> matr, std::vector<T> vec)
	{
		std::vector<T> res(vec.size(), 0);
		for (size_t i = 0; i < vec.size(); ++i)
		{
			for (size_t j = 0; j < vec.size(); j++)
			{
				res[i] += matr[i][j] * vec[j];
			}
		}
		return res;
	}

	template<typename T>
	point<T> multiply_(std::vector<std::vector<T>> matr, std::vector<T> vec)
	{
		std::vector<T> res(vec.size(), 0);
		for (size_t i = 0; i < vec.size(); ++i)
		{
			res[i] = matr[i][i] * vec[i];
		}
		return res;
	}

};

// f(x, y_, z_) = ch(3x) - 10 / (1 + (y_)^2) + ch(z_)
template<typename T>
Ftype f(point<T> p)
{
	return std::cosh(3 * p[0]) - 10 / (1 + y_(p) * y_(p)) + std::cosh(z_(p));
}

// f(x, y_, z_) = 3x^2 - 10 / (1 + (y_)^2) + ch(z_)
template<typename T>
Ftype f_(point<T> p)
{
	return 3 * p[0] * p[0] - 10 / (1 + y_(p) * y_(p)) + std::cosh(z_(p));
}

// f(x, y_, z_) = ch(3x) - y_^2 + ch(z_)
template<typename T>
Ftype f__(point<T> p)
{
	return std::cosh(3 * p[0]) + y_(p) * y_(p) + std::cosh(z_(p));
}

// f(x, y_, z_) = 3x^2 + y_^2 + z_^2
template<typename T>
Ftype f___(point<T> p)
{
	return 3 * p[0] * p[0] + y_(p) * y_(p) + z_(p) * z_(p);
}

// f(x, y_, z_) = 3x^2 + y_^2 + z_^2
template<typename T>
Ftype _f(point<T> p)
{
	return std::cosh(1 * p[0]) + cosh(y_(p) * y_(p)) + std::cosh(z_(p));
}

void method_selector(size_t method)
{
	point<Ftype> a(2, 1.5, -3.5);
	solution result;
	switch (method)
	{
		case(1):
		{
			result = gauss_zeidel_(a, f);
			std::cout << "Минимум, полученный методом Гаусса-Зейделя: " << result.res.to_string() << " | +- " << result.precission << '\n';
			std::cout << "За " << result.iterations << " * 3 итераций" << '\n';
			break;
		}
		case(2):
		{
			result = gradient_descent_(a, f);
			std::cout << "Минимум, полученный методом градиентого спуска: " << result.res.to_string() << " | +- " << result.precission << '\n';
			std::cout << "За " << result.iterations << " итераций" << '\n';
			break;
		}
		case(3):
		{
			result = newton_method(a, f___);
			std::cout << "Минимум, полученный методом Ньютона: " << result.res.to_string() << " | +- " << result.precission << '\n';
			std::cout << "За " << result.iterations << " итераций" << '\n';
			break;
		}
	}
}


size_t get_input()
{
	std::string input;
	size_t res;
	bool flag = true;
	while (flag)
	{
		std::getline(std::cin, input);
		try
		{
			res = std::stoi(input);
			if (res < 1 || res > 3)
			{
				throw std::invalid_argument("");
			}
			flag = false;
		}
		catch (std::invalid_argument e)
		{
			std::cout << "Неверный аргумент! Повторите ввод:" << '\n';
		}
	}
	//system("cls");
	return res;
}

size_t choose_method()
{
	std::cout << "Выберите метод оптимизации:" << '\n';
	std::cout << "-> 1 - Метод Гаусса-Зейделя" << '\n';
	//std::cout << "-> 2 - Метод Хука-Дживса" << '\n';
	std::cout << "-> 2 - Градиентный спуск" << '\n';
	std::cout << "-> 3 - Метод Ньютона" << '\n';
	return get_input();
}

bool continue_flag()
{
	std::cout << "Чтобы выйти нажмите q" << '\n';
	std::cout << "Чтобы продолжить нажмите c" << '\n';
	std::cout << "Чтобы продолжить и очистить вывод введите \\cls" << '\n';
	std::string input;
	bool flag = true;
	while (flag)
	{
		std::getline(std::cin, input);
		if (input != "q" && input != "c" && input != "\\cls")
		{
			std::cout << "Некорректный ввод" << '\n';
		}
		else
		{
			flag = false;
		}
	}
	if (input == "\\cls")
	{
		system("cls");
	}
	return input != "q";
}
#endif // !UTILS_H
