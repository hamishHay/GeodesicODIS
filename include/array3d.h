#ifndef ARRAY3D_H
#define ARRAY3D_H

#include "gridConstants.h"
#include <stdexcept>


template<typename T>
class Array3D {
public:
  Array3D();
  Array3D(unsigned rows, unsigned cols, unsigned depth);

  // template<typename T>

  T& operator() (unsigned row, unsigned col, unsigned depth);        // Subscript operators often come in pairs
  T  operator() (unsigned row, unsigned col, unsigned depth) const;  // Subscript operators often come in pairs
  // ...
 ~Array3D();                              // Destructor
  Array3D(const Array3D& m);               // Copy constructor
  Array3D& operator= (const Array3D& m);   // Assignment operator
  // ...
private:
  unsigned rows_, cols_, depth_;
  T* data_;
};

template<typename T>
inline
Array3D<T>::Array3D() {}



template<typename T>
inline
Array3D<T>::Array3D(unsigned rows, unsigned cols, unsigned depth)
  : rows_ (rows)
  , cols_ (cols)
  , depth_ (depth)
  // , data_ ← initialized below after the if...throw statement
{
  // if (rows == 0 || cols == 0)
  //   throw BadIndex("Array3D constructor has 0 size");
  data_ = new T[rows * cols * depth]();
}

template<typename T>
inline
Array3D<T>::~Array3D()
{
  delete[] data_;
}

template<typename T>
inline
T& Array3D<T>::operator() (unsigned row, unsigned col, unsigned depth)
{
#ifdef DEBUG
  if (row >= rows_ || col >= cols_ || depth >= depth_)
    throw std::out_of_range("Array3D subscript out of bounds");
#endif

  return data_[(cols_*row + col) * depth_ + depth];
  // (x * SIZE2 + y) * SIZE3 + z
}

template<typename T>
inline
T Array3D<T>::operator() (unsigned row, unsigned col, unsigned depth) const
{
#ifdef DEBUG
  if (row >= rows_ || col >= cols_ || depth >= depth_)
    throw std::out_of_range("Array3D subscript out of bounds");
#endif

  return data_[(cols_*row + col) * depth_ + depth];
}

template<typename T>
inline
Array3D<T> & Array3D<T>::Array3D::operator=(const Array3D<T> &other_array)
{
  rows_ = other_array.rows_;
  cols_ = other_array.cols_;
  depth_ = other_array.depth_;
  
  data_ = new T[rows_*cols_*depth_]();

  return *this;
}




#endif
