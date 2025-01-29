#pragma once
#ifndef ONEDIMENSIONAL_H
#define ONEDIMENSIONAL
#include"Assets.h"
#include"Utils.h"

namespace one_dimensional
{
	template<typename T>
	solution fibonachi_kifer(point<T> a, point<T> b, size_t coord_ind, Ftype(*f)(point<T>))
	{
		size_t iterations = 0, counter = 2;
		Ftype val = point<Ftype>::distance(a, b) / constants::eps;
		while (constants::F(iterations) < val)
		{
			++iterations;
		}
		point<Ftype> p1 = a + (b - a) * constants::F(iterations - 2) / constants::F(iterations);
		point<Ftype> p2 = a + (b - a) * constants::F(iterations - 1) / constants::F(iterations);
		Ftype f1 = f(p1);
		Ftype f2 = f(p2);
		Ftype eps = 0;
		while (counter < iterations)
		{
#ifdef DEBUG_1
			std::cout << '\n' << std::string(15, '-') << '\n';
			std::cout << "Fibonachi Iteration: " << std::setw(5) << counter - 1 << '\n';
			std::cout << std::setprecision(4) << "a: " << std::setw(10) << a.to_string() << " | " <<
				"p1: " << std::setw(10) << p1.to_string() << " | " <<
				"p2: " << std::setw(10) << p2.to_string() << " | " <<
				"b: " << std::setw(10) << b.to_string() << " | " <<
				"f1: " << std::setw(10) << f1 << " | " << "f2: " << std::setw(10) << f2 << '\n';
			//getchar();
#endif
			eps = point<Ftype>::distance(a, b) / 100;
			if (f1 > f2)
			{
				a.setCoord(coord_ind, p1[coord_ind]);
				p1.setCoord(coord_ind, p2[coord_ind]);
				p2.setCoord(coord_ind, (a + (b - a) * constants::F(iterations - counter - 1) / constants::F(iterations - counter))[coord_ind]);
				if (point<Ftype>::distance(a, b) < eps)
				{
					p2.updateCoord(coord_ind, eps, SHIFT);
				}
				f1 = f2;
				f2 = f(p2);
			}
			else
			{
				b.setCoord(coord_ind, p2[coord_ind]);
				p2.setCoord(coord_ind, p1[coord_ind]);
				p1.setCoord(coord_ind, (a + (b - a) * constants::F(iterations - counter - 2) / constants::F(iterations - counter))[coord_ind]);
				if (point<Ftype>::distance(a, b) < eps)
				{
					p1.updateCoord(coord_ind, (-1) * eps, SHIFT);
				}
				f2 = f1;
				f1 = f(p1);
			}
			++counter;
		}
		point<Ftype> point_res = b;
		b.setCoord(coord_ind, (a[coord_ind] + b[coord_ind]) / 2);
		solution res(iterations, point_res);
		return res;
	}

	template<typename Functional_>
	Ftype dihotomia_method(Functional_ f, Ftype a, Ftype b)
	{
		if (a > b)
		{
			std::swap(a, b);
		}
		size_t counter = 0;
		while (abs(b - a) > 2 * constants::eps)
		{
#ifdef DEBUG_1
			std::cout << '\n' << std::string(15, '-') << '\n';
			std::cout << "Iteration: " << std::setw(5) << counter << '\n';
			std::cout << "a = " << a << " | b = " << b << " | f(a) = " << f(a) << " | f(b) = " << f(b) << " | b - a = " << b - a << '\n';
#endif // DEBUG_1

			f((b + a) / 2 + constants::eps * 0.1) > f((b + a) / 2 - constants::eps * 0.1) ? b = (b + a) / 2 : a = (b + a) / 2;
			++counter;
		}
		return (b + a) / 2;
	}
}
#endif // !ONEDIMENSIONAL_H
