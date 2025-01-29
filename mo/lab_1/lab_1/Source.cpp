/**
* Целевая функция : f(x) = 2 - cos(x)
* Отрезок : [-pi/4, pi/2]
*/

#include<iostream>
#include<vector>
#include<math.h>
#include<iomanip>
#include<string>

typedef std::vector<std::vector<double>> logger;

double f(double x) 
{
	return 2 - cos(0.1 * x);
}

struct solution
{
	int iterations;
	double res;
	logger logs;

	solution()
	{
		this->iterations = 0;
		this->logs = {};
		this->res = 0;
	}

	solution(int iterations, double res)
	{
		this->iterations = iterations;
		this->res = res;
	}

	solution(int iterations, double res, logger logs)
	{
		this->iterations = iterations;
		this->res = res;
		this->logs = logs;
	}

	void print(int method)
	{
		double precission;
		switch (method)
		{
			case(1):
			{
				std::cout << "Метод дихотомии" << '\n';
				precission = this->logs[this->logs.size() - 1][1] - this->logs[this->logs.size() - 1][0];
				break;
			}
			case(2):
			{
				std::cout << "Метод золотого сечения" << '\n';
				precission = this->logs[this->logs.size() - 1][3] - this->logs[this->logs.size() - 1][0];
				break;
			}
			case(3):
			{
				std::cout << "Метод фибонначи" << '\n';
				precission = this->logs[this->logs.size() - 1][3] - this->logs[this->logs.size() - 1][0];
				break;
			}
		}
		//double result = this->res > 0.00001 ? this->res : 0;
		std::cout << "Найденный минимум f(x): " << std::setprecision(3) << res << " +- " << precission / 2 << " | Количество итераций для нахождения минимума: " << this->iterations << '\n';
	}

	void logger_dihotomia()
	{
		std::cout << "_____________________________________________________________________" << '\n';
		std::cout << "Логи для метода дихотомии" << '\n';
		std::cout << std::setprecision(4);
		for (int i = 0; i < logs.size(); ++i) {
			std::cout << "Итерация " << std::setw(4) << i << ": a = " << std::setw(11) << logs[i][0] << " | b = " << std::setw(11) << logs[i][1] << " | x_mid = " << std::setw(11) << logs[i][2] << '\n';
		}
		std::cout << "_____________________________________________________________________" << '\n' << '\n';
	}

	void logger_golden_ratio()
	{
		std::cout << "_____________________________________________________________________" << '\n';
		std::cout << "Логи для метода золотого сечения" << '\n';
		std::cout << std::setprecision(4);
		for (int i = 0; i < logs.size(); ++i) {
			std::cout << "Итерация " << std::setw(3) << i << ": a = " << std::setw(10) << logs[i][0] << " | x_1 = " << std::setw(10)
				<< logs[i][1] << " | x_2 = " << std::setw(10) << logs[i][2] << " | b = " << std::setw(10) << logs[i][3] << " | f(x_1) = " << std::setw(10) << logs[i][4]
				<< " | f(x_2) = " << std::setw(10) << logs[i][5] << '\n';
		}
		std::cout << "_____________________________________________________________________" << '\n' << '\n';
	}

	void logger_fibonachi_kifer()
	{
		std::cout << "_____________________________________________________________________" << '\n';
		std::cout << "Логи для метода Фибоначчи-Керера" << '\n';
		std::cout << std::setprecision(3);
		for (int i = 0; i < logs.size(); ++i) {
			//" | f(x_1) = " << std::setw(10) << logs[i][4]
			//<< " | f(x_2) = " << std::setw(10) << logs[i][5] <<
			std::cout << "Итерация " << std::setw(3) << i << ": a = " << std::setw(10) << logs[i][0] << " | x_1 = " << std::setw(10)
				<< logs[i][1] << " | x_2 = " << std::setw(10) << logs[i][2] << " | b = " << std::setw(10) << logs[i][3] << " | |a - b|: " << std::setw(10) << logs[i][3] - logs[i][0] << '\n';
		}
		std::cout << "_____________________________________________________________________" << '\n' << '\n';
	}
};

namespace constants 
{
	double pi = 3.14159;
	double eps = 0.00001;
	double psi = 0.61803;
	int F(int num) 
	{
		return num <= 1? 1 : F(num - 2) + F(num - 1);
	}
}

namespace task {
	double min(double a, double mid, double b) {
		double res = 0;
		if (f(a) < f(b) && f(a) < f(mid))
		{
			res = a;
		}
		else if (f(b) < f(a) && f(b) < f(mid))
		{
			res = b;
		}
		else
		{
			res = mid;
		}
		return res;
	}

	solution dihotomia_method(double a, double b) 
	{
		logger logs;
		if (a > b) 
		{
			std::swap(a, b);
		}
		int iterations = 0;
		while (abs(b - a) > 2 * constants::eps)
		{
			double x = (b + a) / 2;
			std::vector<double> log = { a, b, (a + b) / 2 };
			if (f(x + constants::eps * 0.1) > f(x - constants::eps * 0.1)) 
			{ 
				b = x;
			}
			else 
			{
				a = x;
			}
			logs.push_back(log);
			iterations += 2;
		}
		std::vector<double> log = { a, b, (a + b) / 2};
		logs.push_back(log);
		return solution(iterations, min(a, (b + a) / 2, b), logs);
	}

	solution golden_ratio(double a, double b) 
	{
		if (a > b)
		{
			std::swap(a, b);
		}
		int iterations = 2;
		logger logs;
		double x1 = b - (b - a) * constants::psi;
		double x2 = a + (b - a) * constants::psi;
		double f1 = f(x1);
		double f2 = f(x2);
		std::vector<double> log = { a, x1, x2, b, f1, f2 };
		logs.push_back(log);
		while (abs(b - a) > 2 * constants::eps)
		{
			if (f1 > f2) 
			{
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + (b - a) * constants::psi;
				f2 = f(x2);
			}
			else 
			{
				b = x2;
				x2 = x1;
				f2 = f1;
				x1 = b - (b - a) * constants::psi;
				f1 = f(x1);
			}
			log = { a, x1, x2, b, f1, f2 };
			logs.push_back(log);
			++iterations;
		}
		return solution(iterations, min(a, (b + a) / 2, b), logs);
	}

	solution fibonachi_kifer(double a, double b)
	{
		if (a > b)
		{
			std::swap(a, b);
		}
		int iterations = 0, counter = 2;
		while (constants::F(iterations) < (b - a) / constants::eps)
		{
			++iterations;
		}
		logger logs;
		double x1 = a + (b - a) * constants::F(iterations - 2) / constants::F(iterations);
		double x2 = a + (b - a) * constants::F(iterations - 1) / constants::F(iterations);
		double f1 = f(x1);
		double f2 = f(x2);
		double eps = 0;
		std::vector<double> log = {a, x1, x2, b, f1, f2};
		logs.push_back(log);
		while (counter < iterations)
		{
			eps = std::abs(a - b) / 100;
			if (f1 > f2)
			{
				a = x1;
				x1 = x2;
				x2 = a + (b - a) * constants::F(iterations - counter - 1) / constants::F(iterations - counter);
				if (std::abs(x1 - x2) < eps)
				{
					x2 += eps;
				}
				f1 = f2;
				f2 = f(x2);
			}
			else
			{
				b = x2;
				x2 = x1;
				x1 = a + (b - a) * constants::F(iterations - counter - 2) / constants::F(iterations - counter);
				if (std::abs(x1 - x2) < eps)
				{
					x1 -= eps;
				}
				f2 = f1;
				f1 = f(x1);
			}
			log = { a, x1, x2, b, f1, f2};
			logs.push_back(log);
			++counter;
		}
		return solution(iterations, min(a, (b + a) / 2, b), logs);
	}
}

void method_selector(int method, bool is_logger)
{
	double a = (-10) * constants::pi / 4;
	double b = 5 * constants::pi / 2;
	solution result;
	switch (method)
	{
		case(1):
		{
			result = task::dihotomia_method(a, b);
			result.print(1);
			if (is_logger)
			{
				result.logger_dihotomia();
			}
			break;
		}
		case(2):
		{
			result = task::golden_ratio(a, b);
			result.print(2);
			if (is_logger)
			{
				result.logger_golden_ratio();
			}
			break;
		}
		case(3):
		{
			result = task::fibonachi_kifer(a, b);
			result.print(3);
			if (is_logger)
			{
				result.logger_fibonachi_kifer();
			}
			break;
		}
	}
}

int get_input()
{
	std::string input;
	int res;
	bool flag = true;
	while (flag)
	{
		std::getline(std::cin, input);
		try
		{
			res = std::stoi(input);
			if (res < 1 || res > 3)
			{
				throw std::invalid_argument("incorrect argument");
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

int choose_method()
{
	std::cout << "Выберите метод оптимизации:" << '\n';
	std::cout << "-> 1 - Метод дихотомии" << '\n';
	std::cout << "-> 2 - Метод золотого сечения" << '\n';
	std::cout << "-> 3 - Метод Фибоначчи-Кифера" << '\n';
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

int main() 
{
	setlocale(LC_ALL, "Russian");
	do
	{
		std::cout << "_____________________________________________________________________" << '\n' << '\n';
		method_selector(choose_method(), 1);
		std::cout << "_____________________________________________________________________" << '\n' << '\n';
	} while (continue_flag());
	return 0;
}