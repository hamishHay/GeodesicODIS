#ifndef ARRAY3D_H
#define ARRAY3D_H

template<typename T>
class Array3D {
public:
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
Array3D<T>::Array3D(unsigned rows, unsigned cols, unsigned depth)
  : rows_ (rows)
  , cols_ (cols)
  , depth_ (depth)
  // , data_ ← initialized below after the if...throw statement
{
  // if (rows == 0 || cols == 0)
  //   throw BadIndex("Array3D constructor has 0 size");
  data_ = new T[rows * cols * depth];
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
  // if (row >= rows_ || col >= cols_)
  //   throw BadIndex("Array3D subscript out of bounds");
  // return data_[cols_*row + col];
  return data_[(cols_*row + col) * depth_ + depth];
  // (x * SIZE2 + y) * SIZE3 + z
}

template<typename T>
inline
T Array3D<T>::operator() (unsigned row, unsigned col, unsigned depth) const
{
  // if (row >= rows_ || col >= cols_)
  //   throw BadIndex("const Matrix subscript out of bounds");
  return data_[(cols_*row + col) * depth_ + depth];
}


#endif
