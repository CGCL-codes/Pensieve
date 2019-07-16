#ifndef __VECTOR_HPP__
#define __VECTOR_HPP__

#define _a Alloc()

#include <memory>
#include <iostream>
#include "myutil.h"


template< typename T, typename Alloc = std::allocator<T> >

class myVector : public std::allocator<T> {
public:
  // Constructor
  myVector() {
    _data = _a.allocate(1);
    _size = 0;
    _capacity = 1;
  }

  myVector(const std::size_t & size, const T & val, Alloc a = Alloc()) {
    if (size == 0) {
      _data = a.allocate(1);
      _size = 0;
      _capacity = 1;
    } else {
      _data = a.allocate(size);
      for (std::size_t i = 0; i < size; ++i)
        a.construct(_data + i, val);
      _size = _capacity = size;
    }
  }
  template<typename InputIterator>
  myVector(InputIterator begin, InputIterator end, Alloc a = Alloc()) {
    if ((end - begin) == 0) {
      _data = _a.allocate(1);
      _size = 0;
      _capacity = 1;
    } else {
      _size = _capacity = end - begin;
      _data = a.allocate(_size);
      std::size_t cnt = 0;
      for (InputIterator it = begin; it != end; ++it)
        a.construct(_data + (cnt++), *it);
    }
  }
  myVector(const myVector & other) {
    _size = other._size;
    _capacity = other._capacity;
    _data = _a.allocate(_capacity);
    for (std::size_t i = 0; i < _size; ++i)
      _a.construct(_data + i, other._data[i]);
  }

  // Destructor
  ~myVector() {
    for (std::size_t i = 0; i < _size; ++i)
      _a.destroy(_data + i);
    if (_capacity > 0)
      _a.deallocate(_data, _capacity);
  }

  // Copy Operator
  myVector & operator=(const myVector & other) {
    if (&other != this) {
      std::size_t i;
      //std::cout << "size" << _size << " capa" << _capacity << std::endl;
      for (i = 0; i < _size; ++i)
        _a.destroy(_data + i);
      //std::cout << "size" << _size << " capa" << _capacity << std::endl;
      if (_capacity > 0)
        _a.deallocate(_data, _capacity);
      //std::cout << "size" << _size << " capa" << _capacity << std::endl;
      _size = other._size;
      _capacity = other._capacity;
      _data = _a.allocate(_capacity);
      for (i = 0; i < _size; ++i)
        _a.construct(_data + i, other._data[i]);
    }
    return *this;
  }

  // Iterator
  typedef T *         iterator;
  typedef const T *   const_iterator;

  inline iterator begin() {
    return _data;
  }
  inline const_iterator begin() const {
    return _data;
  }
  inline iterator end() {
    return _data + _size;
  }
  inline const_iterator end() const {
    return _data + _size;
  }

  // Capacity
  inline std::size_t size() const {
    return _size;
  }
  void resize(const std::size_t & newSize) {
    std::size_t i;
    if (newSize <= _size) {
      for (i = newSize; i < _size; ++i)
        _a.destroy(_data + i);
      }
    else {
      if (newSize > _capacity) {
        std::size_t newCapacity = (std::size_t) (_capacity * get_factor());
        while (newSize > newCapacity)
          newCapacity = (std::size_t) (newCapacity*get_factor());
        reserve(newCapacity);
      }
      for (i = _size; i < newSize; ++i)
        _a.construct(_data + i, T());
    }
    _size = newSize;
  }
  void resize(const std::size_t & newSize, const T & val) {
    std::size_t i;
    if (newSize <= _size) {
      for (i = newSize; i < _size; ++i)
        _a.destroy(_data + i);
      }
    else {
      if (newSize > _capacity) {
        std::size_t newCapacity = (std::size_t) (_capacity * get_factor());
        while (newSize > newCapacity)
          newCapacity = (std::size_t) (newCapacity*get_factor());
        reserve(newCapacity);
      }
      for (i = _size; i < newSize; ++i)
        _a.construct(_data + i, val);
    }
    _size = newSize;
  }
  inline std::size_t capacity() const {
    return _capacity;
  }
  inline bool empty() const {
    return _size == 0;
  }
  void reserve(const std::size_t & newCapacity) {
    if (newCapacity > _capacity) {
      T * temp = _a.allocate(newCapacity);
      for (std::size_t i = 0; i < _size; ++i) {
        _a.construct(temp + i, _data[i]);
        _a.destroy(_data + i);
      }
      _a.deallocate(_data, _capacity);
      _capacity = newCapacity;
      _data = temp;
    }
  }

  // Element Access
  inline T & operator[](const std::size_t & index) {
    return _data[index];
  }
  inline const T & operator[](const std::size_t & index) const {
    return _data[index];
  }
  inline T & front() {
    return _data[0];
  }
  inline const T & front() const {
    return _data[0];
  }
  inline T & back() {
    return _data[_size - 1];
  }
  inline const T & back() const {
    return _data[_size - 1];
  }
  inline T * data() {
    return _data;
  }
  inline const T * data() const {
    return _data;
  }

  // Modifiers
  template<typename InputIterator>
  void assign(InputIterator begin, InputIterator end) {
    std::size_t newSize = 0;
    InputIterator it;
    for (it = begin; it != end; ++it)
      ++newSize;
    if (newSize > _capacity) {
      std::size_t newCapacity = (std::size_t) (_capacity * get_factor());
      while (newSize > newCapacity)
        newCapacity = (std::size_t) (newCapacity*get_factor());
      reserve(newCapacity);
    }
    std::size_t i;
    for (i = 0; i < _size; ++i)
      _a.destroy(_data + i);
    for (i = 0, it = begin; i < newSize; ++i, ++it)
      _a.construct(_data + i, *it);
    _size = newSize;
  }
  void assign(const std::size_t & newSize, const T & val) {
    if (newSize > _capacity) {
      std::size_t newCapacity = (std::size_t) (_capacity * get_factor());
      while (newSize > newCapacity)
        newCapacity = (std::size_t) (newCapacity*get_factor());
      reserve(newCapacity);
    }
    std::size_t i;
    for (i = 0; i < _size; ++i)
      _a.destroy(_data + i);
    for (i = 0; i < newSize; ++i)
      _a.construct(_data + i, val);
    _size = newSize;
  }
  void push_back(const T & val) {
    if (_size >= _capacity){
      reserve((std::size_t) (_capacity * get_factor()));
    }
    _a.construct(_data + (_size++), val);
  }

  void careful_push_back(const T & val) {
    if (_size >= _capacity) {
      reserve(_capacity + 1);
    }
    _a.construct(_data + (_size++), val);
  }

  void push_back(iterator begin, iterator end) {
    std::size_t newSize = 0;
    iterator it;
    
    for (it=begin; it != end; it ++) {
      ++ newSize;
    }
    // std::cout << _size << " " << newSize << "\n";
    if (newSize + _size > _capacity) {
      reserve(newSize + _size);
    }
    std::size_t i;
    for (it = begin, i=0; i < newSize; i++, it++) {
      _a.construct(_data+_size+i, *it);
    }
    _size = _size + newSize;
  }

  void pop_back() {
    // T back = _data[_size -1];
    _a.destroy(_data + (--_size));
    // return back;
  }//

  void clear() {
    myVector<T>().swap(*this);
  }
  //delete element according to index
  void index_delete(std::size_t index) {
    if (index >= _size){
      std::cout << "size erro in delete " << index << " of " << _size << std::endl;
      abort();
    }
    else {
      _data[index] = _data[_size - 1];
      _a.destroy(_data + (--_size));
    }
  }
  //add element in backtracking
  int index_addtion(const T & val, std::size_t index) {
    push_back(val);
    if (index >= _size) {
      std::cout << "size erro in add" << _size << " " << index << std::endl;
      abort();
      return -1;
    }
    // else {
      swap(index, _size - 1);
      return _size;
    // }
  }
  //swap element
  void swap(const std::size_t index1, const std::size_t index2) {
    if (index1 < _size && index2 < _size) {
      T tmp = _data[index1];
      _data[index1] = _data[index2];
      _data[index2] = tmp;
    }
    else
      std::cout << "erro" << std::endl;
  }
  //swap vector
  void swap(myVector<T> & _v) {
    std::swap(_data, _v._data);
    std::swap(_size, _v._size);
    std::swap(_capacity, _v._capacity);
  }
  //find index of element
  int find(const T & val) {
    for (std::size_t i = 0; i < _size; i++)
    {
      if (val == _data[i])
        return i;
    }
    return -1;
  }

  void printdata(int limit) {
    std::cout << "vector contains : ";
    for(auto i = 0; i < limit && i < _size; i++) {
      std::cout << _data[i] << " ";
    }
    std::cout << std::endl;
  }

  void printdata(){
    printdata(_size);
  }

  int check_increase() {
    if (_size <= 1) {
      return 0;
    }
    for (auto i = 1; i < _size; i++) {
      if (_data[i-1] > _data[i]) {
        return i;
      }
    }
    return 0;
  }

  double get_factor() {
    if (_capacity > 1000) {
      return 1.1;
    } else if (_capacity > 500) {
      return 1.2;
    } else if (_capacity > 100) {
      return 1.5;
    } else {
      return 2.0;
    }
  }

private:
  iterator _data;
  std::size_t _size, _capacity;
};

// template <typename T, typename A>

#endif
