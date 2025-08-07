/*! \file Array3D.cpp
 *  \brief 3D array
 */
#include <cstdlib>
#include <vector>
using namespace std;

namespace cusg {
/*! \class Array3D
 *  \brief 3D array
 */  
template<class T>
class Array3D {
    size_t width, height;
    vector<T> data;
  public:
    Array3D(): width(0), height(0) { };

    Array3D(size_t x, size_t y, size_t z):
      width(x), height(y)
    {
        data.resize(x*y*z);
    }

    T& operator()(size_t x, size_t y, size_t z) {
        return data.at(x + y * width + z * width * height);
    }
    
    T& at(size_t x, size_t y, size_t z) {
        return data.at(x + y * width + z * width * height);
    }

    void clear() {
        data.clear();
    }

    void init(size_t x, size_t y, size_t z) {
        data.clear();
        width = x;
        height = y;
        data.resize(x*y*z);
    }
};
}   // end namespace cusg
