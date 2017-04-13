#include "../include/Matrix3x3.hpp"

Matrix3x3::Matrix3x3() {}

Matrix3x3::Matrix3x3(std::initializer_list< std::initializer_list<float> > input_list)
{
  unsigned dim = 0;
  for(auto column_iter : input_list)    
  {
    std::copy( column_iter.begin(), column_iter.end(), mat[dim].begin() );
    dim++;
  }
}

Matrix3x3::Matrix3x3(std::array< std::array<float,3>, 3 > input_matrix)
{
  this->mat = input_matrix;
}

Matrix3x3::Matrix3x3(std::array< ThreeVector<float>, 3> input_matrix)
{
  unsigned row_index = 0;
  for(auto input_row : input_matrix)
  {
    this->mat[row_index] = input_row.GetArray();
    row_index++;
  }
}

Matrix3x3 Matrix3x3::Identity()
{
  Matrix3x3 UnityMatrix = {{1,0,0},{0,1,0},{0,0,1}};
  return UnityMatrix;
}

std::array<float,3>& Matrix3x3::operator[](const long unsigned int index)
{
  return this->mat[index];
}

Matrix3x3 Matrix3x3::operator+=(const Matrix3x3& right_mat)
{
  for(unsigned row = 0; row < this->mat.size(); row++)
  {
    for(unsigned column = 0; column < this->mat[row].size(); column++)
    {
      this->mat[row][column] += right_mat.mat[row][column];
    }
  }
  return *this;
}

Matrix3x3 Matrix3x3::operator-=(const Matrix3x3& right_mat)
{
  for(unsigned row = 0; row < this->mat.size(); row++)
  {
    for(unsigned column = 0; column < this->mat[row].size(); column++)
    {
      this->mat[row][column] -= right_mat.mat[row][column];
    }
  }
  return *this;
}

Matrix3x3 operator+(Matrix3x3 left_mat, const Matrix3x3& right_mat)
{
  return left_mat += right_mat;
}

Matrix3x3 operator-(Matrix3x3 left_mat, const Matrix3x3& right_mat)
{
  return left_mat -= right_mat;
}

Matrix3x3 Matrix3x3::operator*(const Matrix3x3& mat_right)
{
  Matrix3x3 ResultMatrix = {{0,0,0},{0,0,0},{0,0,0}};
  for(unsigned row = 0; row < this->mat.size(); row++)
  {
    for(unsigned column = 0; column < this->mat[row].size(); column++)
    {
      for(unsigned entry = 0; entry < this->mat.size(); entry++)
      {
	ResultMatrix.mat[row][column] += this->mat[row][entry]*mat_right.mat[entry][column];
      }
    }
  }
  return ResultMatrix;
}

ThreeVector<float> Matrix3x3::operator*(ThreeVector<float>& right_vec)
{
  ThreeVector<float> result = {0.0,0.0,0.0};
  for(unsigned row = 0; row < this->mat.size(); row++)
  {
    for(unsigned column = 0; column < this->mat[row].size(); column++)
    {
      result[row] += this->mat[row][column]*right_vec[column];
    }
  }
  return result;
}

Matrix3x3 Matrix3x3::operator*=(const float& scaling_factor)
{
  for(unsigned row = 0; row < this->mat.size(); row++)
  {
    for(unsigned column = 0; column < this->mat[row].size(); column++)
    {
      this->mat[row][column]*= scaling_factor;
    }
  }
  return *this;
}

Matrix3x3 operator*(Matrix3x3 left_mat, const float& right_scal)
{
  return left_mat *= right_scal;
}

Matrix3x3 operator*(const float& left_scal, Matrix3x3 right_mat)
{
  return right_mat *= left_scal;
}

const float Matrix3x3::Determinant()
{
  float det = this->mat[0][0]*(this->mat[1][1]*this->mat[2][2] - this->mat[1][2]*this->mat[2][1])
            + this->mat[0][1]*(this->mat[1][2]*this->mat[2][0] - this->mat[1][0]*this->mat[2][2])
	    + this->mat[0][2]*(this->mat[1][0]*this->mat[2][1] - this->mat[1][1]*this->mat[2][0]);
  return det;
}

const float Matrix3x3::Trace()
{
  float tr = 0;
  for(unsigned dim = 0; dim < this->mat.size(); dim++) tr += this->mat[dim][dim];
  return tr;
}

Matrix3x3 Matrix3x3::Inverse()
{
  // Cayley-Hamilton decomposition
  if(this->Determinant() != 0.0) 
    return 1/this->Determinant()*( 0.5*(pow(this->Trace(),2) - (*this * *this).Trace()) * this->Identity() - *this*this->Trace() + *this * *this );
  else
  {
    std::cout << "Warning: Matrix is Singular!" << std::endl;
    return *this;
  }
}

bool Matrix3x3::Invert()
{ 
  if(this->Determinant() != 0.0)
  {
    *this = this->Inverse();
    return true;
  }
  else return false;
}
