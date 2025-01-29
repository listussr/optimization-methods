#pragma once
#ifndef ASSETS_H
#define ASSETS_H

#include<iostream>
#include<vector>
#include<iomanip>
#include<string>
#include<limits>
#include<cmath>
#define SHIFT true

#if defined(F32)
#define Ftype float
#elif defined(F64)
#define Ftype double
#else 
#define Ftype long double
#endif

#endif // !ASSETS_H

