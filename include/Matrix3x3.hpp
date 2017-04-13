#include <array>
#include <initializer_list>
#include <algorithm>
#include <cmath>

#include "ThreeVector.hpp"

#ifndef MATRIX3X3_H
#define MATRIX3X3_H

class Matrix3x3
{
private:
  std::array< std::array<float,3>,3 > mat;
  
public:
  Matrix3x3();
  Matrix3x3(std::initializer_list< std::initializer_list<float> >);
  Matrix3x3(std::array< std::array<float,3>,3 >);
  Matrix3x3(std::array<ThreeVector<float>,3>);
  
  Matrix3x3 Identity();
  
  std::array<float,3>& operator[](const unsigned long);
  ThreeVector<float> operator*(ThreeVector<float>&);
  Matrix3x3 operator+=(const Matrix3x3&);
  Matrix3x3 operator-=(const Matrix3x3&);
  friend Matrix3x3 operator+(Matrix3x3, const Matrix3x3&);
  friend Matrix3x3 operator-(Matrix3x3, const Matrix3x3&);
  
  Matrix3x3 operator*(const Matrix3x3&);
  
  Matrix3x3 operator*=(const float&);
  friend Matrix3x3 operator*(Matrix3x3, const float&);
  friend Matrix3x3 operator*(const float&, Matrix3x3);
  
  const float Determinant();
  const float Trace();
  
  Matrix3x3 Inverse();
  bool Invert();
};

// #include "Matrix3x3.cpp"

#endif