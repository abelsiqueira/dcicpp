#ifndef dci_extra_h
#define dci_extra_h

namespace DCI {

  template < typename T > 
  void delpointer (T * p) {
    if (p != 0)
      delete p;
  }

  template < typename T > 
  void delarray (T * p) {
    if (p != 0)
      delete []p;
  }

  template < typename T >
  T Max (T a, T b) {
    if (a > b)
      return a;
    else
      return b;
  }

  template <typename T >
  T Min (T a, T b) {
    if (a > b)
      return b;
    else 
      return a;
  }

  template <typename T >
  T AbsValue (T x) {
    if (x > 0)
      return x;
    else
      return -x;
  }

}

#endif
