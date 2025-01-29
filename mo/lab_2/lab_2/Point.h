#pragma once
#ifndef POINT_H
#define POINT_H
/// <summary>
/// класс дл€ работы с вектором координат
/// </summary>
#include"Assets.h"
template<typename T>
class point
{
private:

	std::vector<T>	coords;
	Ftype			alpha;
	std::vector<T>	shifts;

public:

	point()
	{
		coords = {};
	}

	/// <summary>
	/// „астный случай дл€ данной задачи
	/// </summary>
	/// <param name="x">- координата</param>
	/// <param name="y">- координата</param>
	/// <param name="z">- координата</param>
	point(T x, T y, T z)
	{
		coords.push_back(x);
		coords.push_back(y);
		coords.push_back(z);
	}

	point(const point<T>& p)
	{
		this->coords = p.coords;
	}

	point(point<T> p, T val, size_t ind)
	{
		this->coords		= p.coords;
		this->coords[ind]	= val;
	}

	point(const std::vector<T>& coordinates) : coords(coordinates) {}

	void setTransformInfo(std::vector<T> shifts, Ftype alpha)
	{
		this->alpha			= alpha;
		this->shifts		= shifts;
	}

	/// <summary>
	/// ”станавливаем новые координаты
	/// </summary>
	/// <returns></returns>
	void setCoordinates(std::vector<T> coordinates)
	{
		this->coords(coordinates);
	}

	/// <summary>
	/// ”станавливаем новые координаты по другой точке, мен€€ в ней координату по индексу "ind"
	/// </summary>
	/// <param name="ind"></param>
	/// <returns></returns>
	void setCoordinates(point<T> p, T val, size_t ind)
	{
		this->coords		= p.coords;
		this->coords[ind]	= val;
	}

	/// <summary>
	/// ѕолучаем координаты в форме std::vector<T>
	/// </summary>
	/// <returns></returns>
	std::vector<T> getCoordinates()
	{
		return this->coords;
	}

	/// <summary>
	/// ѕолучаем трансформированные координаты в форме std::vector<T>
	/// </summary>
	/// <returns></returns>
	std::vector<T> transform()
	{
		std::vector<T> new_coord(this->coords.size());
		new_coord[0] = this->coords[0];
		new_coord[1] = (this->coords[1] - this->shifts[0]) * cos(this->alpha) + (this->coords[2] - this->shifts[1]) * sin(this->alpha);
		new_coord[2] = (this->coords[2] - this->shifts[1]) * cos(this->alpha) - (this->coords[1] - this->shifts[0]) * sin(this->alpha);
		return point<T>(new_coord);
	}

	std::vector<T> getTransformedCoordinates()
	{
		return this->transform(this->coords);
	}

	const T& operator[](size_t ind) const
	{
		if (ind < dimensionality()) {
			return coords[ind];
		}
		else {
			return coords[0];
		}
	}

	/// <summary>
	/// размерность пространства
	/// </summary>
	/// <returns></returns>
	size_t dimensionality() const
	{
		return this->coords.size();
	}

	/// <summary>
	/// ћен€ем значение координаты по индексу
	/// </summary>
	/// <param name="index"></param>
	/// <param name="shift"></param>
	/// <param name="is_shift"></param>
	void setCoord(size_t index, T new_value)
	{
		if (index > dimensionality() - 1 || index < 0)
		{
			throw std::invalid_argument("Ќекорректный индекс!");
		}
		this->coords[index] = new_value;
	}

	/// <summary>
	/// —двигаем все координаты на значение shift
	/// </summary>
	/// <returns></returns>
	void updateCoord(size_t index, Ftype shift, bool is_shift)
	{
		if (index > dimensionality() - 1 || index < 0)
		{
			throw std::invalid_argument("Ќекорректный индекс!");
		}
		this->coords[index] += shift;
	}

	/// <summary>
	/// ”множаем все координаты на скал€р
	/// </summary>
	/// <param name="p"></param>
	/// <returns></returns>
	void updateCoord(size_t index, Ftype scale)
	{
		if (index > dimensionality() - 1 || index < 0)
		{
			throw std::invalid_argument("Ќекорректный индекс!");
		}
		this->coords[index] *= scale;
	}

	/// <summary>
	/// Ќормализаци€ координат
	/// </summary>
	void normalize()
	{
		Ftype len = getLength();
		(*this) /= len;
	}

	/// <summary>
	/// ѕолучаем рассто€ние между данной точкой и точкой (0, ..., 0)
	/// </summary>
	/// <param name="p"></param>
	/// <returns></returns>
	Ftype getLength()
	{
		Ftype square_len = 0;
		for (size_t i = 0; i < dimensionality(); ++i)
		{
			square_len += std::pow(coords[i], 2);
		}
		return sqrt(square_len);
	}

	point<T> operator+(const point<T>& p) const {
		if (dimensionality() != p.dimensionality()) {
			std::cerr << "Error: Cannot add points with different dimensions!" << std::endl;
			return point<T>();
		}

		std::vector<T> result(dimensionality());
		for (size_t i = 0; i < dimensionality(); ++i) {
			result[i] = coords[i] + p.coords[i];
		}
		return point<T>(result);
	}

	point<T> operator+(Ftype shift) const
	{
		std::vector<T> new_coord = {};
		for (size_t i = 0; i < dimensionality(); ++i)
		{
			new_coord.push_back(this->coords[i] + shift);
		}
		return point<T>(new_coord);
	}

	point<T> operator-(const point<T>& p) const
	{
		if (p.dimensionality() != this->dimensionality())
		{
			throw std::range_error("Ќесоответсвие размерности векторов");
		}
		std::vector<T> new_coord = {};
		for (size_t i = 0; i < dimensionality(); ++i)
		{
			new_coord.push_back(coords[i] - p[i]);
		}
		return point<T>(new_coord);
	}

	point<T> operator-(Ftype shift) const
	{
		std::vector<T> new_coord = {};
		for (size_t i = 0; i < dimensionality(); ++i)
		{
			new_coord.push_back(this->coords[i] + shift);
		}
		return point<T>(new_coord);
	}

	point<T> operator*(Ftype scale) const
	{
		std::vector<T> new_coord = {};
		for (size_t i = 0; i < dimensionality(); ++i)
		{
			new_coord.push_back(this->coords[i] * scale);
		}
		return point<T>(new_coord);
	}

	T operator*(const point<T>& p) const
	{
		if (p.dimensionality() != this->dimensionality())
		{
			throw std::range_error("Ќесоответсвие размерности векторов");
		}
		T res{ 0 };
		for (size_t i = 0; i < dimensionality(); ++i)
		{
			res += this->coords[i] * p[i];
		}
		return res;
	}

	void operator/=(Ftype scale)
	{
		if (scale < 0.000001)
		{
			throw std::invalid_argument("ƒеление на 0");
		}
		for (size_t i = 0; i < dimensionality(); ++i)
		{
			this->coords[i] /= scale;
		}
	}

	point<T> operator/(Ftype scale)
	{
		if (scale < 0.000001)
		{
			throw std::invalid_argument("ƒеление на 0");
		}
		for (size_t i = 0; i < dimensionality(); ++i)
		{
			this->coords[i] /= scale;
		}
		return *this;
	}

	void operator+=(Ftype shift)
	{
		for (size_t i = 0; i < dimensionality(); ++i)
		{
			this->coords[i] += shift;
		}
	}

	void operator-=(Ftype shift)
	{
		for (size_t i = 0; i < dimensionality(); ++i)
		{
			this->coords[i] -= shift;
		}
	}

	void operator+=(std::vector<T> shifts)
	{
		if (shifts.size() != this->dimensionality())
		{
			throw std::range_error("Ќесоответсвие размерности векторов");
		}
		for (size_t i = 0; i < dimensionality(); ++i)
		{
			this->coords[i] += shifts[i];
		}
	}

	void operator-=(std::vector<T> shifts)
	{
		if (shifts.size() != this->dimensionality())
		{
			throw std::range_error("Ќесоответсвие размерности векторов");
		}
		for (size_t i = 0; i < dimensionality(); ++i)
		{
			this->coords[i] -= shifts[i];
		}
	}

	/// <summary>
	/// –ассто€ние между 2 точками в n-мерном пространстве
	/// </summary>
	/// <param name="p"> - точка в n-мерном пространстве</param>
	/// <returns></returns>
	static Ftype distance(point<T> a, point<T> b)
	{
		if (a.dimensionality() != b.dimensionality())
		{
			throw std::range_error("Ќесоответсвие размерности векторов");
		}
		Ftype square_dist = 0;
		for (size_t i = 0; i < a.dimensionality(); ++i)
		{
			square_dist += std::pow((a[i] - b[i]), 2);
		}
		return std::sqrt(square_dist);
	}

	/// <summary>
	/// строковое представление координат
	/// </summary>
	/// <returns></returns>
	std::string to_string() const
	{
		std::string str = "(";
		for (size_t i = 0; i < coords.size(); ++i)
		{
			str += std::to_string(coords[i]);
			if (i != coords.size() - 1)
			{
				str += ", ";
			}
		}
		str += ")";
		return str;
	}
};
#endif // !POINT_H
