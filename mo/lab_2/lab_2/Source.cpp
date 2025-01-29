#define F64
#define DEBUG
//#define DEBUG_1
#include"MultiDimensional.h"


/**
* cosh(3x) - 10 / (1 + y'^2) + cosh(z'); 
* (y_0, z_0) = (1 , -4); 
* \alpha = \pi/3; 
* X_0 = (2, 3, -1);
*/
int main()
{
	setlocale(LC_ALL, "Russian");
	do
	{
		std::cout << std::string(50, '-') << '\n';
		method_selector(choose_method());
		std::cout << std::string(50, '-') << '\n';
	} while (continue_flag());
	return 0;
}