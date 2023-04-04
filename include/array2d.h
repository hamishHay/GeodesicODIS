#ifndef ARRAY2D_H
#define ARRAY2D_H

#include "gridConstants.h"
#include <stdexcept>

template<typename T>
class Array2D {
public:
  Array2D();
  Array2D(unsigned rows, unsigned cols);

  // template<typename T>

  T& operator() (unsigned row, unsigned col);        // Subscript operators often come in pairs
  T  operator() (unsigned row, unsigned col) const;  // Subscript operators often come in pairs
  // ...
 ~Array2D();                              // Destructor
  Array2D(const Array2D& m);               // Copy constructor
  Array2D& operator= (const Array2D& m);   // Assignment operator
  // ...
private:
  unsigned rows_, cols_;
  T* data_;
};

template<typename T>
inline
Array2D<T>::Array2D()
  // , data_ ← initialized below after the if...throw statement
{
  // if (rows == 0 || cols == 0)
  //   throw BadIndex("Array2D constructor has 0 size");
//   data_ = new T[rows * cols]();
}

template<typename T>
inline
Array2D<T>::Array2D(unsigned rows, unsigned cols)
  : rows_ (rows)
  , cols_ (cols)
  // , data_ ← initialized below after the if...throw statement
{
  // if (rows == 0 || cols == 0)
  //   throw BadIndex("Array2D constructor has 0 size");
  data_ = new T[rows * cols]();
}

template<typename T>
inline
Array2D<T>::~Array2D()
{
  delete[] data_;
}

template<typename T>
inline
T& Array2D<T>::operator() (unsigned row, unsigned col)
{

#ifdef DEBUG
  if (row >= rows_ || col >= cols_)
    throw std::out_of_range("Array2D subscript out of bounds");
#endif

  return data_[cols_*row + col];
}

template<typename T>
inline
T Array2D<T>::operator() (unsigned row, unsigned col) const
{

#ifdef DEBUG
  if (row >= rows_ || col >= cols_)
    throw std::out_of_range("Array2D subscript out of bounds");
#endif

  // if (row >= rows_ || col >= cols_)
  //   throw BadIndex("const Matrix subscript out of bounds");
  return data_[cols_*row + col];
}

template<typename T>
inline
Array2D<T>::Array2D(const Array2D<T> &other_array)
    : rows_ (other_array.rows_),
      cols_ (other_array.cols_)
{
  
  // if (row >= rows_ || col >= cols_)
  //   throw BadIndex("const Matrix subscript out of bounds");
  data_ = new T[rows_*cols_]();
//   return data_ = ;
}

template<typename T>
inline
Array2D<T> & Array2D<T>::Array2D::operator=(const Array2D<T> &other_array)
{
  rows_ = other_array.rows_;
  cols_ = other_array.cols_;
  // if (row >= rows_ || col >= cols_)
  //   throw BadIndex("const Matrix subscript out of bounds");
  data_ = new T[rows_*cols_]();
//   return &this;
  return *this;
}


// Element & Element::operator=(const Element &other_element)
// {
//   this->xyz_coords[0] = other_element.xyz_coords[0];
//   this->xyz_coords[1] = other_element.xyz_coords[1];
//   this->xyz_coords[2] = other_element.xyz_coords[2];

//   // Maybe we shouldn't do this here?
//   cart2sph(xyz_coords, sph_coords);
//   // sph_coords[1] = pi*0.5 - sph_coords[1];

//   return *this;
// }


#endif
