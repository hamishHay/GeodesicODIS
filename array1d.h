#ifndef ARRAY1D_H
#define ARRAY1D_H

template<typename T>
class Array1D {
public:
  Array1D(unsigned rows);

  // template<typename T>

  T& operator() (unsigned row);        // Subscript operators often come in pairs
  T  operator() (unsigned row) const;  // Subscript operators often come in pairs
  // ...
 ~Array1D();                              // Destructor
  Array1D(const Array1D& m);               // Copy constructor
  Array1D& operator= (const Array1D& m);   // Assignment operator
  // ...
private:
  unsigned rows_;
  T* data_;
};


template<typename T>
inline
Array1D<T>::Array1D(unsigned rows)
  : rows_ (rows)
  // , data_ ← initialized below after the if...throw statement
{
  // if (rows == 0 || cols == 0)
  //   throw BadIndex("Array1D constructor has 0 size");
  data_ = new T[rows]();
}

template<typename T>
inline
Array1D<T>::~Array1D()
{
  delete[] data_;
}

template<typename T>
inline
T& Array1D<T>::operator() (unsigned row)
{
  // if (row >= rows_ || col >= cols_)
  //   throw BadIndex("Array1D subscript out of bounds");
  return data_[row];
}

template<typename T>
inline
T Array1D<T>::operator() (unsigned row) const
{
  // if (row >= rows_ || col >= cols_)
  //   throw BadIndex("const Matrix subscript out of bounds");
  return data_[row];
}


#endif
