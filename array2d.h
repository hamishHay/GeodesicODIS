#ifndef ARRAY2D_H
#define ARRAY2D_H



// class Array2D {
// public:
//   Array2D(unsigned rows, unsigned cols);
//
//   //
//
//   double& operator() (unsigned row, unsigned col);        // Subscript operators often come in pairs
//   double  operator() (unsigned row, unsigned col) const;  // Subscript operators often come in pairs
//   // ...
//  ~Array2D();                              // Destructor
//   Array2D(const Array2D& m);               // Copy constructor
//   Array2D& operator= (const Array2D& m);   // Assignment operator
//   // ...
// private:
//   unsigned rows_, cols_;
//   double* data_;
// };
//
//
//
// inline
// Array2D::Array2D(unsigned rows, unsigned cols)
//   : rows_ (rows)
//   , cols_ (cols)
//   // , data_ ← initialized below after the if...throw statement
// {
//   // if (rows == 0 || cols == 0)
//   //   throw BadIndex("Array2D constructor has 0 size");
//   data_ = new double[rows * cols];
// }
//
//
// inline
// Array2D::~Array2D()
// {
//   delete[] data_;
// }
//
//
// inline
// double& Array2D::operator() (unsigned row, unsigned col)
// {
//   // if (row >= rows_ || col >= cols_)
//   //   throw BadIndex("Array2D subscript out of bounds");
//   return data_[cols_*row + col];
// }
//
//
// inline
// double Array2D::operator() (unsigned row, unsigned col) const
// {
//   // if (row >= rows_ || col >= cols_)
//   //   throw BadIndex("const Matrix subscript out of bounds");
//   return data_[cols_*row + col];
// }

template<typename T>
class Array2D {
public:
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
Array2D<T>::Array2D(unsigned rows, unsigned cols)
  : rows_ (rows)
  , cols_ (cols)
  // , data_ ← initialized below after the if...throw statement
{
  // if (rows == 0 || cols == 0)
  //   throw BadIndex("Array2D constructor has 0 size");
  data_ = new T[rows * cols];
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
  // if (row >= rows_ || col >= cols_)
  //   throw BadIndex("Array2D subscript out of bounds");
  return data_[cols_*row + col];
}

template<typename T>
inline
T Array2D<T>::operator() (unsigned row, unsigned col) const
{
  // if (row >= rows_ || col >= cols_)
  //   throw BadIndex("const Matrix subscript out of bounds");
  return data_[cols_*row + col];
}


#endif
