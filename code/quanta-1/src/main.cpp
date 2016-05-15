/* Quanta 1 - physics simulation of colliding balls
 * Copyright (C) 2016 Nolan Eakins
 * All rights reserved.
 */

#include <cmath>
#include <climits>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <SDL2/SDL.h>
#define GL_GLEXT_PROTOTYPES
#if MAC
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include <assert.h>
#include <stdexcept>
#include <functional>
#include <unistd.h>

#include "math-fun.hpp"

#ifdef CHEAPER_MATH
typedef float Real;

#define ator atof

#define fmodr fmodf
#define acosr acosf
#define floorr floorf
#define ceilr ceilf
#define modfr modff
#define atan2r atan2f
#define sqrtr sqrtf
#define rabs fabs

#define glVertex3r glVertex3f
#define glVertex3rv glVertex3fv
#define glNormal3rv glNormal3fv
#define glLoadMatrixr glLoadMatrixf
#define glMultMatrixr glMultMatrixf
#define glScaler glScalef
#define glTranslater glTranslatef
#define glColor3rv glColor3fv
#define glColor4rv glColor4fv
#define glLightrv glLightfv
#define glLightr glLightf
#define glMaterialrv glMaterialfv
#define glFogrv glFogfv

#else
typedef double Real;

#define ator atof

#define fmodr fmod
#define acosr acos
#define floorr floor
#define ceilr ceil
#define modfr modf
#define atan2r atan2
#define sqrtr sqrt
#define rabs fabs

#define glVertex3r glVertex3d
#define glVertex3rv glVertex3dv
#define glNormal3rv glNormal3dv
#define glLoadMatrixr glLoadMatrixd
#define glMultMatrixr glMultMatrixd
#define glScaler glScaled
#define glTranslater glTranslated
#define glColor3rv glColor3dv
#define glColor4rv glColor4dv
#define glLightr glLightf

void glLightrv(GLenum light, GLenum a, const Real v[4])
{
  float fv[4] = { (float)v[0], (float)v[1], (float)v[2], (float)v[3] };
  glLightfv(light, a, fv);
}

void glFogrv(GLenum name, const Real v[4])
{
  float fv[4] = { (float)v[0], (float)v[1], (float)v[2], (float)v[3] };
  glFogfv(name, fv);
}

void glMaterialrv(GLenum face, GLenum a, const Real v[4])
{
  float fv[4] = { (float)v[0], (float)v[1], (float)v[2], (float)v[3] };
  glMaterialfv(face, a, fv);
}

#endif /* CHEAPER_MATH */

#define DELEGATE(type, method, member) \
  type method() const { return member.method(); }

#define ASSERT_GL_ERROR { \
    GLenum error = glGetError(); \
    if(error != GL_NO_ERROR) { \
      std::cerr << "GL error:" << __LINE__ << ": " << gluErrorString(error) << std::endl; \
    } \
  }

#define warn(condition) { \
  if(!(condition)) { std::cout << "Warning:" << __LINE__ << ": " << #condition << std::endl; } \
}

#define APP_TITLE "Quanta 1"
#define APP_LOGO "quanta.bmp"
#define APP_WIDTH 800
#define APP_HEIGHT 600
const std::string APP_TEXTURE_PATH("assets/textures/");

class IndexOutOfBoundsError
{
  const void *_array;
  int _index, _max;

public:
  IndexOutOfBoundsError(const void *array, int index, int max)
    : _array(array), _index(index), _max(max)
  {
  }

  const void *array() const { return _array; }
  int index() const {  return _index; }
  int max() const { return _max; }
};

#ifdef MAC
using std::isnan;
using std::isinf;
#endif

class Vec3;
Vec3 operator* (Real a, const Vec3 &v);

class Vec3
{
#ifdef SLOW_VECTOR
  typedef Real RealArray[4];
#else
  typedef Real RealArray __attribute__((vector_size(sizeof(Real) * 4)));
  typedef int BoolArray __attribute__((vector_size(sizeof(int) * 4)));
#endif

  RealArray _v;
  
public:
  static const Vec3 X, Y, Z, W, Origin;
  static const Vec3 nX, nY, nZ, nW;

  Vec3(const RealArray other)
  {
#ifdef SLOW_VECTOR
    for(int i = 0; i < 4; i++) {
      _v[i] = other[i];
    }
#else
    _v = other;
#endif
  }
  
  Vec3(Real x = 0.0f, Real y = 0.0f, Real z = 0.0f, Real w = 1.0f)
  {
    _v[0] = x;
    _v[1] = y;
    _v[2] = z;
    _v[3] = w;
  }
  Vec3(const Vec3 &other)
  {
    set(other);
  }

  Real x() const { return _v[0]; }
  void setX(Real n) { _v[0] = n; }
  Real y() const { return _v[1]; }
  void setY(Real n) { _v[1] = n; }
  Real z() const { return _v[2]; }
  void setZ(Real n) { _v[2] = n; }
  Real w() const { return _v[3]; }
  void setW(Real n) { _v[3] = n; }

#ifdef SLOW_VECTOR
#define SCALAR_OP(OP) \
  template<typename T> \
  Vec3 operator OP (const T &n) const \
  { \
    return Vec3(_v[0] OP n, _v[1] OP n, _v[2] OP n, _v[3] OP n); \
  }
#else
#define SCALAR_OP(OP) \
  template<typename T> \
  Vec3 operator OP (const T &n) const \
  { \
    return Vec3(_v OP (Real)n);                 \
  }
#endif

  SCALAR_OP(+);
  SCALAR_OP(-);
  SCALAR_OP(*);
  SCALAR_OP(/);

  Vec3 &operator += (const Vec3 &other)
  {
#ifdef SLOW_VECTOR
    for(int i = 0; i < 4; i++) _v[i] += other._v[i];
#else
    _v = _v + other._v;
#endif
    return *this;
  }

  Vec3 &operator -= (const Vec3 &other)
  {
#ifdef SLOW_VECTOR
    for(int i = 0; i < 4; i++) _v[i] -= other._v[i];
#else
    _v = _v - other._v;
#endif
    return *this;
  }

  Vec3 operator - () const
  {
    return *this * -1.0f;
  }

#undef SCALAR_OP

  Vec3 operator % (Real n) const
  {
    return Vec3(fmodr(x(), n), fmodr(y(), n), fmodr(z(), n), fmodr(w(), n));
  }

#ifdef SLOW_VECTOR
#define VECTOR_OP(OP) \
  Vec3 operator OP (const Vec3 &other) const \
  { \
    return Vec3(_v[0] OP other._v[0], _v[1] OP other._v[1], _v[2] OP other._v[2], _v[3] OP other._v[3]); \
  }
#else
#define VECTOR_OP(OP) \
  Vec3 operator OP (const Vec3 &other) const \
  { \
    return Vec3(_v OP other._v); \
  }
#endif

  VECTOR_OP(+);
  VECTOR_OP(-);
  VECTOR_OP(*);

#undef VECTOR_OP

  bool operator == (const Vec3 &other) const
  {
    //return (*this - other).magnitude() < 0.001;
    // return fabs(x() - other.x()) < 0.001 &&
    //   fabs(y() - other.y()) < 0.001 &&
    //   fabs(z() - other.z()) < 0.001 &&
    //   fabs(w() - other.w()) < 0.001;
#ifdef SLOW_VECTOR
    return x()==other.x() && y() == other.y() && z() == other.z(); // && w() == other.w();
#else
    BoolArray a = _v == other._v;
    return a[0] && a[1] && a[2];
#endif
  }

  bool operator != (const Vec3 &other) const
  {
#ifdef SLOW_VECTOR
    return x()!=other.x() || y() != other.y() || z() != other.z(); // || w() != other.w();
#else
    BoolArray a = _v != other._v;
    return a[0] || a[1] || a[2];
#endif
  }

  Vec3 cross(const Vec3 &other) const
  {
#ifdef SLOW_VECTOR
    return Vec3(_v[1] * other._v[2] - _v[2] * other._v[1],
                _v[2] * other._v[0] - _v[0] * other._v[2],
                _v[0] * other._v[1] - _v[1] * other._v[0]);
#else
    RealArray a = { _v[1], _v[2], _v[0] };
    RealArray b = { other._v[2], other._v[0], other._v[1] };
    RealArray c = { _v[2], _v[0], _v[1] };
    RealArray d = { other._v[1], other._v[2], other._v[0] };
    return Vec3(a * b - c * d);
#endif
  }

  Real dot(const Vec3 &other) const
  {
#ifdef SLOW_VECTOR
    return x()*other.x() + y()*other.y() + z()*other.z();
#else
    RealArray a = _v * other._v;
    return a[0] + a[1] + a[2];
#endif
  }

  Vec3 normalize() const
  {
    return *this / magnitude() * Vec3(1.0, 1.0, 1.0, 0.0);;
  }

  Real magnitude() const
  {
    return sqrt(magnitude_squared());
  }

  Real magnitude_squared() const
  {
#ifdef SLOW_VECTOR
    return (_v[0] * _v[0] + _v[1] * _v[1] + _v[2] * _v[2]);
#else
    RealArray a = _v * _v;
    return a[0] + a[1] + a[2];
#endif
  }

  Real scalarProjection(const Vec3 &v) const
  {
    return dot(v.normalize());
  }

  Vec3 projectOnto(const Vec3 &v) const
  {
    return v.normalize() * scalarProjection(v);
  }

  Vec3 rejection(const Vec3 &v) const
  {
    return *this - projectOnto(v);
  }

  Vec3 reflectBy(const Vec3 &normal) const
  {
    return -2.0f * projectOnto(normal) + *this;
  }

  Vec3 negate() const
  {
    return -(*this);
  }

  bool isNaN() const {
    for(int i = 0; i < 4; i++) {
      if(isnan(_v[i])) return true;
    }
    return false;
  }

  const Real *get() const { return (Real *)&_v; }

  void set(const Vec3 &other)
  {
#ifdef SLOW_VECTOR
    for(int i = 0; i < 4; i++) {
      _v[i] = other._v[i];
    }
#else
    _v = other._v;
#endif
  }

  operator const Real * () const
  {
    return get();
  }

  Vec3 clamp(Real min = 0.0f, Real max = 1.0f) const
  {
    return Vec3(::clamp(x(), max, min),
                ::clamp(y(), max, min),
                ::clamp(z(), max, min),
                ::clamp(w(), max, min));
  }

  Vec3 clamp(const Vec3 &min, const Vec3 &max) const
  {
    return Vec3(::clamp(x(), max.x(), min.x()),
                ::clamp(y(), max.y(), min.y()),
                ::clamp(z(), max.z(), min.z()),
                ::clamp(w(), max.w(), min.w()));
  }

  Vec3 abs() const
  {
    return Vec3(fabs(x()), fabs(y()), fabs(z()), fabs(w()));
  }

  static Real randf()
  {
    return (std::rand() & 0xFFFFFF) / (Real)0xFFFFFF;
  }

  static Vec3 random(const Vec3 &min, const Vec3 &max)
  {
    Vec3 delta = max - min;
    Vec3 offset = delta / 2.0f;
    return Vec3(delta.x() != 0.0f ? randf() * delta.x() : 0.0f,
                delta.y() != 0.0f ? randf() * delta.y() : 0.0f,
                delta.z() != 0.0f ? randf() * delta.z() : 0.0f,
                delta.w() != 0.0f ? randf() * delta.w() : 0.0f) - offset;
  }

  static Vec3 random()
  {
    return random(Vec3(-1024.0, -1024.0, -1024.0),
                  Vec3(1024.0, 1024.0, 1024.0)).normalize();
  }

  Vec3 floor() const
  {
    return Vec3(::floorr(x()), ::floorr(y()), ::floorr(z()), ::floorr(w()));
  }
};

const Vec3 Vec3::X(1, 0, 0), Vec3::Y(0, 1, 0), Vec3::Z(0, 0, 1), Vec3::W(0, 0, 0, 1), Vec3::Origin(0, 0, 0);
const Vec3 Vec3::nX(-1, 0, 0), Vec3::nY(0, -1, 0), Vec3::nZ(0, 0, -1), Vec3::nW(0, 0, 0, -1);

Vec3 operator* (int a, const Vec3 &v)
{
  return v * a;
}

Vec3 operator* (float a, const Vec3 &v)
{
  return v * (Real)a;
}

Vec3 operator* (double a, const Vec3 &v)
{
  return v * (Real)a;
}

class Matrix
{
#ifdef SLOW_VECTOR
  typedef Real RealArray[16];
#else
  typedef Real RealArray __attribute__((vector_size(sizeof(Real) * 16)));
#endif
  RealArray _values;

public:
  Matrix()
  {
#if 1
    for(int i = 0; i < 4; i++) {
      for(int j = 0; j < 4; j++) {
        _values[i*4+j] = (Real)0.0;
      }
    }
#else
    _values = { (Real)0.0, (Real)0.0, (Real)0.0, (Real)0.0,
                (Real)0.0, (Real)0.0, (Real)0.0, (Real)0.0,
                (Real)0.0, (Real)0.0, (Real)0.0, (Real)0.0,
                (Real)0.0, (Real)0.0, (Real)0.0, (Real)0.0
    };
#endif
  }

  Matrix(const Matrix &m)
  {
#ifdef SLOW_VECTOR    
    for(int i = 0; i < 4; i++) {
      for(int j = 0; j < 4; j++) {
        _values[i*4+j] = m._values[i*4+j];
      }
    }
#else
    _values = m._values;
#endif
  }

  Matrix(const RealArray other)
  {
#ifdef SLOW_VECTOR
    for(int i = 0; i < 4; i++) {
      for(int j = 0; j < 4; j++) {
        _values[i*4+j] = other[i*4+j];
      }
    }
#else
    _values = other;
#endif
  }

  Matrix(Real a, Real b, Real c, Real d,
         Real e, Real f, Real g, Real h,
         Real i, Real j, Real k, Real l,
         Real m, Real n, Real o, Real p)
  {
    set(0, 0, a);
    set(1, 0, b);
    set(2, 0, c);
    set(3, 0, d);

    set(0, 1, e);
    set(1, 1, f);
    set(2, 1, g);
    set(3, 1, h);

    set(0, 2, i);
    set(1, 2, j);
    set(2, 2, k);
    set(3, 2, l);

    set(0, 3, m);
    set(1, 3, n);
    set(2, 3, o);
    set(3, 3, p);
  }

#ifdef SLOW_VECTOR
#define BINARY_OP(OP) \
  Matrix operator OP (const Matrix &m) const \
  { \
    Matrix ret; \
    \
    for(int i = 0; i < 4; i++) { \
      for(int j = 0; j < 4; j++) { \
        ret.set(i, j, get(i, j) OP m.get(i, j)); \
      } \
    } \
    \
    return ret;\
  }
#else
#define BINARY_OP(OP) \
  Matrix operator OP (const Matrix &m) const \
  { \
    return Matrix(_values OP m._values); \
  }

#endif
  BINARY_OP(+);
  BINARY_OP(-);

#undef BINARY_OP

  Matrix operator*(const Matrix &m) const
  {
    Matrix ret;

    for(int i = 0; i < 4; i++) {
      for(int j = 0; j < 4; j++) {
        for(int k = 0; k < 4; k++) {
          ret.set(i, j, ret.get(i, j) + get(i, k) * m.get(k, j));
        }
      }
    }

    return ret;
  }

  Vec3 operator*(const Vec3 &v) const
  {
    return Vec3(get(0, 0) * v.x() + get(1, 0) * v.y() + get(2, 0) * v.z() + get(3, 0) * v.w(),
                get(0, 1) * v.x() + get(1, 1) * v.y() + get(2, 1) * v.z() + get(3, 1) * v.w(),
                get(0, 2) * v.x() + get(1, 2) * v.y() + get(2, 2) * v.z() + get(3, 2) * v.w(),
                get(0, 3) * v.x() + get(1, 3) * v.y() + get(2, 3) * v.z() + get(3, 3) * v.w());
  }

  Matrix operator*(Real n) const
  {
#ifdef SLOW_VECTOR
    Matrix ret;
    for(int i = 0; i < 4; i++) {
      for(int j = 0; j < 4; j++) {
        ret.set(i, j, get(i, j) * n);
      }
    }
    return ret;
#else
    return Matrix(_values * n);
#endif
  }

  Matrix operator/(Real n) const
  {
#ifdef SLOW_VECTOR
    Matrix ret;
    for(int i = 0; i < 4; i++) {
      for(int j = 0; j < 4; j++) {
        ret.set(i, j, get(i, j) / n);
      }
    }
    return ret;
#else
    return Matrix(_values / n);
#endif
  }

  Matrix transpose() const
  {
    Matrix ret;
    for(int i = 0; i < 4; i++) {
      for(int j = 0; j < 4; j++) {
        ret.set(i, j, get(j, i));
      } 
    }
    return ret;
  }

  Real sign(int x, int y) const
  {
    return (x*4+(y+x%2)) % 2 ? (Real)-1.0 : (Real)1.0;
  }

  Real sign3x3(int x, int y) const
  {
    return (x*3+y) % 2 ? (Real)-1.0 : (Real)1.0;
  }

  Real cofactor(int x, int y) const
  {
    return sign(x, y) * minorM(x, y).determinant3x3();
  }

  Real cofactor2x2(int x, int y) const
  {
    return sign3x3(x, y) * minorM(x, y).determinant2x2();
  }

  Matrix cofactorMatrix() const
  {
    Matrix m;
    for(int i = 0; i < 4; i++) {
      for(int j = 0; j < 4; j++) {
        m.set(i, j, cofactor(i, j));
      }
    }
    return m;
  }

  Real determinant() const
  {
    Real d = (Real)0.0;

    for(int i = 0; i < 4; i++) {
      d += get(i, 0) * cofactor(i, 0);
    }

    return d;
  }

  Real determinant3x3() const
  {
    Real d = (Real)0.0;
    for(int i = 0; i < 3; i++) {
      d += get(i, 0) * cofactor2x2(i, 0);
    }
    return d;
  }

  Real determinant2x2() const
  {
    return get(0,0) * get(1,1) - get(1,0) * get(0,1);
  }

  Matrix minorM(int col, int row) const 
  {
    Matrix r;

    for(int i = 0, x = 0; i < 3; i++, x++) {
      if(x == col) x++;
      for(int j = 0, y = 0; j < 3; j++, y++) {
        if(y == row) y++;
        r.set(i, j, get(x, y));
      }
    }
    return r;
  }

  Matrix adjoint() const
  {
    return cofactorMatrix().transpose();
  }

  Matrix invert() const
  {
    return adjoint() * ((Real)1.0 / determinant());
  }

  bool hasNaN() const
  {
    for(int i = 0; i < 4; i++) {
      for(int j = 0; j < 4; j++) {
        if(isnan(get(i, j))) {
          return true;
        }
      }
    }

    return false;
  }

  static Matrix rotateOnto(const Vec3 &a, const Vec3 &b)
  {
    return rotate(acosr(a.dot(b)), a.cross(b));
  }

  // https://en.wikipedia.org/wiki/Cross_product#Conversion_to_matrix_multiplication
  static Matrix cross(const Vec3 &v)
  {
    return Matrix(0, -v.z(), v.y(), 0,
                  v.z(), 0, -v.x(), 0,
                  -v.y(), v.x(), 0, 0,
                  0, 0, 0, 0);
  }

  static Matrix rotate(Real angle, const Vec3 &axis)
  {
    Real c = cos(angle), s = sin(angle), t = (Real)1.0 - c;
    Real x = axis.x(), y = axis.y(), z = axis.z();

    return Matrix(t*x*x + c,   t*y*x - s*z, t*z*x + s*y, 0,
                  t*x*y + s*z, t*y*y + c  , t*z*y - s*x, 0,
                  t*x*z - s*y, t*y*z + s*x, t*z*z + c,   0,
                  0, 0, 0, 1);
  }

  static Matrix rotateX(Real angle)
  {
    Real s = sin(angle), c = cos(angle);

    return Matrix(1, 0, 0, 0,
                  0, c, -s, 0,
                  0, s, c, 0,
                  0, 0, 0, 1);
  }

  static Matrix rotateY(Real angle)
  {
    Real s = sin(angle), c = cos(angle);

    return Matrix(c, 0, s, 0,
                  0, 1, 0, 0,
                  -s, 0, c, 0,
                  0, 0, 0, 1);
  }

  static Matrix rotateZ(Real angle)
  {
    Real s = sin(angle), c = cos(angle);

    return Matrix(c, -s, 0, 0,
                  s, c, 0, 0,
                  0, 0, 1, 0,
                  0, 0, 0, 1);
  }

  static Matrix rotate(const Vec3 &angles)
  {
    //return rotateX(angles.x()) * rotateY(angles.y()) * rotateZ(angles.z());
    return rotateZ(angles.z()) * rotateX(angles.x()) * rotateY(angles.y());
  }

  static Matrix translation(const Vec3 &v)
  {
    return Matrix(1, 0, 0, v.x(),
                  0, 1, 0, v.y(),
                  0, 0, 1, v.z(),
                  0, 0, 0, 1); 
  }

  static Matrix translation(Real x, Real y, Real z)
  {
    return translation(Vec3(x, y, z));
  }

  static Matrix orient(const Vec3 &r, const Vec3 &u, const Vec3 &f, const Vec3 &o = Vec3::Origin)
  {
    return Matrix(r.x(), u.x(), -f.x(), o.x(),
                  r.y(), u.y(), -f.y(), o.y(),
                  r.z(), u.z(), -f.z(), o.z(),
                  0, 0, 0, 1);
  }

  static Matrix lookAt(const Vec3 &forward, const Vec3 &up)
  {
      Vec3 axis = forward.normalize();
      Vec3 right = up.cross(axis).normalize();
      Vec3 real_up = axis.cross(right).normalize();
      return orient(right, real_up, axis);
  }

  static Matrix lookAt(const Vec3 &target, const Vec3 &from, const Vec3 &up)
  {
      return lookAt(target - from, up) * translation(from);
  }

  static Matrix scale(const Vec3 &v)
  {
    return Matrix(v.x(), 0, 0, 0,
                  0, v.y(), 0, 0,
                  0, 0, v.z(), 0,
                  0, 0, 0, v.w());
  }

  static Matrix identity()
  {
    return Matrix(1, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0, 1, 0,
                  0, 0, 0, 1);
  }

  Real get(int x, int y) const
  {
    if(x < 0 || x >= 4)
      throw IndexOutOfBoundsError(this, x, 4);
    if(y < 0 || y >= 4)
      throw IndexOutOfBoundsError(this, y, 4);

    return _values[x*4+y];
  }

  void set(int x, int y, Real v)
  {
    if(x < 0 || x >= 4)
      throw IndexOutOfBoundsError(this, x, 4);
    if(y < 0 || y >= 4)
      throw IndexOutOfBoundsError(this, y, 4);

    _values[x*4+y] = v;
  }

  const Real *get() const {
    return (Real *)&_values;
  }

  operator const Real * () const
  {
    return get();
  }
};

std::ostream &operator<<(std::ostream &s, const Vec3 &v)
{
  s << "<Vec3: " << v.x() << ", " << v.y() << ", " << v.z() << ", " << v.w() << ">";
  return s;
}

std::ostream &operator<<(std::ostream &s, const Matrix &m)
{
  s << "<Matrix: ";

  for(int i = 0; i < 4; i++) {
    s << "[ " << m.get(0, i) << ", " << m.get(1, i) << ", " << m.get(2, i) << ", " << m.get(3, i) << " ], ";
  }

  s << ">";

  return s;
}

namespace Colors
{
  const Vec3 Clear(1, 1, 1, 0);
  const Vec3 Opaque(1, 1, 1, 1);
  const Vec3 White(1, 1, 1);
  const Vec3 Black(0, 0, 0);
  const Vec3 Red(1, 0, 0);
  const Vec3 Green(0, 1, 0);
  const Vec3 Blue(0, 0, 1);
  const Vec3 Yellow(1, 1, 0);
  const Vec3 Cyan(0, 1, 1);
  const Vec3 Magenta(1, 0, 1);

  const int NumColors(6);

  const Vec3 &color(int n)
  {
    static const Vec3 *palette[NumColors] = { &Red, &Yellow, &Green, &Cyan, &Blue, &Magenta };
    return *palette[abs(n) % NumColors];
  }

  const Vec3 &color()
  {
    return color(rand());
  }

  Vec3 color(Real n)
  {
    if(isinf(n)) {
      return White;
    } else if(isnan(n)) {
      return Black;
    } else {
      Vec3 a = color((int)floorr(n)), b = color((int)ceilr(n));
      Real i, t = modfr(n, &i);
      return a * ((Real)1.0 - t) + b * t;
    }
  }
};

class Material
{
  Vec3 _diffuse, _specular, _emission, _ambient;
  int _shininess;
  
public:
  Material(const Vec3 &diffuse, const Vec3 &specular, int shininess, const Vec3 &ambient)
    : _diffuse(diffuse), _specular(specular), _emission(), _shininess(shininess), _ambient(ambient)
  {
  }

  const Vec3 &diffuse() const { return _diffuse; }
  const Vec3 &specular() const { return _specular; }
  const Vec3 &emission() const { return _emission; }
  const Vec3 &ambient() const { return _ambient; }
  int shininess() const { return _shininess; }
};

class Texture
{
  std::string _path;
  bool _repeat_x, _repeat_y;
  SDL_Surface *_surface;

public:
  Texture(const std::string &path, bool repeat_x = true, bool repeat_y = true)
    : _path(path), _repeat_x(repeat_x), _repeat_y(repeat_y), _surface(NULL)
  {
  }

  virtual ~Texture()
  {
    unload();
  }

  const std::string &path() const { return _path; }

  bool loaded() const { return _surface != NULL; }

  void unload()
  {
    if(_surface) {
      SDL_FreeSurface(_surface);
      _surface = NULL;
    }
  }

  bool load()
  {
    SDL_Surface *bmp = SDL_LoadBMP(_path.c_str());
    if(bmp != NULL) {
      unload();
      _surface = bmp;
      return true;
    } else {
      std::cerr << "Error loading " << _path << ": " << SDL_GetError() << std::endl;
      return false;
    }
  }

  SDL_Surface *surface() const { return _surface; }

  bool repeatsX() const { return _repeat_x; }
  bool repeatsY() const { return _repeat_y; }

  int width() const { return _surface->w; }
  int height() const { return _surface->h; }
};

class TexturedMaterial
{
  Material _material;
  Texture _texture;

public:
  TexturedMaterial(const std::string &path, const Vec3 &diffuse = Colors::White, const Vec3 &specular = Colors::Clear, int shininess = 0, const Vec3 &ambient = Colors::Black)
    : _material(diffuse, specular, shininess, ambient), _texture(path, true, true)
  {
  }

  const Material &material() const { return _material; }
  Material &material() { return _material; }

  const Texture &texture() const { return _texture; }
  Texture &texture() { return _texture; }

  DELEGATE(const Vec3 &, diffuse, _material);
};

class ArgumentError
{
public:
  ArgumentError() { }
};

class NotFoundError
{
public:
  const void *root;

  NotFoundError(const void *r, const void*) : root(r) { }
};

class ReallyNotFoundError
{
public:
  ReallyNotFoundError(const void *, const void*) { }
};

class Triangle
{
  Vec3 _verts[3], _normal;

public:
  Triangle(const Vec3 &a, const Vec3 &b, const Vec3 &c)
  {
    _verts[0] = a;
    _verts[1] = b;
    _verts[2] = c;
    updateNormal();
  }

  Vec3 project(const Vec3 &point) const
  {
    return point.projectOnto(normal());
  }

  Vec3 baryCoordinate(const Vec3 &point) const
  {
    Vec3 d = project(point - _verts[0]);
    Vec3 v = point - d;
    Real d_dot = d.dot(normal());
    Real d_sign = d_dot < 0.0f ? -1.0f : 1.0f;
    Vec3 v0 = _verts[1] - _verts[0], v1 = _verts[2] - _verts[0], v2 = v - _verts[0];

    Real d00 = v0.dot(v0),
      d01 = v0.dot(v1),
      d11 = v1.dot(v1),
      d20 = v2.dot(v0),
      d21 = v2.dot(v1);

    Real invDenom = 1.0f / (d00 * d11 - d01 * d01);
    Real s = (d11 * d20 - d01 * d21) * invDenom;
    Real t = (d00 * d21 - d01 * d20) * invDenom;

    return Vec3(s, d_sign * d.magnitude(), t);
  }

  const Vec3 &vertex(int n) const
  {
    if(n < 0 || n >= 3) throw IndexOutOfBoundsError(this, n, 3);
    return _verts[n];
  }

  void setVertex(int n, const Vec3 &v)
  {
    _verts[n] = v;
  }

  void updateNormal()
  {
    _normal = (_verts[1] - _verts[0]).cross(_verts[2] - _verts[0]).normalize();
  }

  const Vec3 &normal() const
  {
    return _normal;
  }

  Vec3 cartCoordinate(const Vec3 &bary) const
  {
    Real t = 1.0f - bary.x() - bary.z();

    return Vec3(t * _verts[0].x() + bary.x() * _verts[1].x() + bary.z()*_verts[2].x(),
                t * _verts[0].y() + bary.x() * _verts[1].y() + bary.z()*_verts[2].y(),
                t * _verts[0].z() + bary.x() * _verts[1].z() + bary.z()*_verts[2].z()) +
      Vec3(normal() * bary.y());
  }

  static bool sameSideAs(const Vec3 &p1, const Vec3 &p2, const Vec3 &a, const Vec3 &b)
  {
    Vec3 cp1 = (b - a).cross(p1 - a);
    Vec3 cp2 = (b - a).cross(p2 - a);
    return cp1.dot(cp2) >= 0.0f;
  }

  bool contains(const Vec3 &v) const
  {
    //    0
    //    /\
    //   / v\
    // 2+----+1
    Vec3 b = baryCoordinate(v);
    return
      (b.x() >= (Real)0.0 && b.x() <= (Real)1.0) &&
      (b.z() >= (Real)0.0 && b.z() <= (Real)1.0) &&
      b.x() + b.z() <= (Real)1.01;
  }

  Vec3 closestPoint(const Vec3 &v) const
  {
    Vec3 bary = baryCoordinate(v);
    Real u = clamp(1.0 - bary.x() - bary.z(), 1.0);
    bary.setX(clamp(bary.x(), (Real)1.0));
    bary.setY(clamp(bary.y(), (Real)0.0));
    bary.setZ(clamp(bary.z(), (Real)(1.0 - bary.x() - u)));
    return cartCoordinate(bary);
  }

  Real distanceTo(const Vec3 &v) const
  {
    return (v - closestPoint(v)).magnitude();
  }

  Real distanceToSquared(const Vec3 &v) const
  {
    return (v - closestPoint(v)).magnitude_squared();
  }
};

class Camera
{
  Matrix _matrix;

public:
  Camera(const Vec3 &position, const Vec3 &forward = Vec3::nZ, const Vec3 &up = Vec3::Y)
  {
    setPosition(position);
    lookAt(forward, up);
  }

  Vec3 position() const //{ return _position; }
  {
    return Vec3(_matrix.get(3,0), _matrix.get(3,1), _matrix.get(3,2));
  }

  void setPosition(const Vec3 &v)
  {
    for(int i = 0; i < 3; i++) {
      _matrix.set(3, i, v[i]);
    }
    _matrix.set(3, 3, 1);
  }

  Real x() const { return position().x(); }
  Real y() const { return position().y(); }
  Real z() const { return position().z(); }

  void setMatrix(const Matrix &m)
  {
    _matrix = m;
  }

  virtual void moveBy(const Vec3 &delta)
  {
    //_position += _matrix * delta;
    _matrix = Matrix::translation(delta) * _matrix;
  }

  void increasePosition(const Vec3 &delta)
  {
    //_position += delta;
    for(int i = 0; i < 3; i++) {
      _matrix.set(3, i, _matrix.get(3, i) + delta[i]);
    }
  }

  void transformTo(const Matrix &m)
  {
    _matrix = _matrix * m;
  }

  void transformBy(const Matrix &m)
  {
    _matrix = m * _matrix;
  }

  virtual void rotateBy(const Vec3 &angle) {
    transformBy(Matrix::rotate(angle));
  }

  //const Vec3 &angle() const { return _angle; }

  Real yaw() const { return M_PI + acosr(Vec3::X.dot(right())); }

  /*
  Real pitch() const { return _angle.x(); }
  Real setPitch(Real n) { _angle.setX(n); updateDirections(); }
  Real setYaw(Real n) { _angle.setY(n); updateDirections(); }
  Real roll() const { return _angle.z(); }
  Real setRoll(Real n) { _angle.setZ(n); updateDirections(); }
  */

  virtual void lookAt(const Vec3 &forward, const Vec3 &up)
  {
    _matrix = Matrix::lookAt(forward, up) * translationMatrix();
  }

  void setUp(const Vec3 &v)
  {
    lookAt(forward(), v);
  }

  Vec3 up() const {
    return Vec3(_matrix.get(1, 0), _matrix.get(1, 1), _matrix.get(1, 2));
  }

  void setForward(const Vec3 &v)
  {
    lookAt(v, up());
  }

  Vec3 forward() const {
    // todo figure why this needs to be negated
    return -Vec3(_matrix.get(2, 0), _matrix.get(2, 1), _matrix.get(2, 2));
  }

  Vec3 right() const
  {
    return Vec3(_matrix.get(0, 0), _matrix.get(0, 1), _matrix.get(0, 2));
  }

  Matrix matrix() const
  {
    return _matrix;
    //return rotationMatrix() * translationMatrix();
  }

  Matrix rotationMatrix() const
  {
    Matrix m = _matrix * translationMatrix().invert();
    return m;
  }

  void setRotationMatrix(const Matrix &m)
  {
    _matrix = m * translationMatrix();
  }

  Matrix translationMatrix() const
  {
    //Real n = negate ? -1.0 : 1.0;
    return Matrix::translation(position());
  }

  Vec3 focalPoint(Real render_distance_ratio) const
  {
    return position() + forward() * render_distance_ratio * 0.5f;
  }

  bool canSee(const Vec3 &point, const Vec3 &normal) const
  {
    return isFacing(normal) && !wouldClip(point);
  }

  bool isFacing(const Vec3 &normal) const
  {
    return forward().dot(normal) <= 0.2f;
  }

  bool wouldClip(const Vec3 &point) const
  {
    Vec3 v = matrix().invert() * point;
    return v.z() > 0.0;
  }
};

class FPSCamera: public Camera
{
  Vec3 _angles, _plane;
  bool _flying;

public:
  FPSCamera(const Vec3 &position, const Vec3 &angles)
    : Camera(position), _angles(angles), _plane(Vec3::Y), _flying(false)
  {
    updateMatrix();
  }

  virtual void moveBy(const Vec3 &delta)
  {
    increasePosition(Matrix::rotate(_angles.y(), _plane) * delta);
  }

  void setAngles(const Vec3 &angles)
  {
    _angles = angles.clamp(Vec3(-M_PI/2.0, -INFINITY, -INFINITY), Vec3(M_PI/2.0, INFINITY, INFINITY));
    updateMatrix();
  }

  virtual void rotateBy(const Vec3 &angles)
  {
    setAngles(_angles + angles);
  }

  virtual void lookAt(const Vec3 &forward, const Vec3 &up)
  {
    Matrix m(Matrix::lookAt(forward, up));
    _plane = up;
    setAngles(Vec3(atan2r(m.get(2,1), m.get(2,2)),
                   atan2r(-m.get(2,0), sqrtr(m.get(2,1)*m.get(2,1)+m.get(2,2)*m.get(2,2))),
                   M_PI + atan2r(m.get(1,0), m.get(0, 0))));
  }

  void updateMatrix()
  {
    setRotationMatrix(Matrix::rotate(_angles) * Matrix::orient(_plane.cross(Vec3::Z), _plane, -Vec3::Z));
  }

  const Vec3 &plane() const { return _plane; }

  void setPlane(const Vec3 &n)
  {
    _plane = n.normalize();
    updateMatrix();
  }
};

class Renderer
{
public:
  typedef std::map<const Texture *, GLuint> TextureIdMap;
private:
  TextureIdMap _texture_handles;
public:
  Renderer() { }

  Renderer(const Renderer &r)
  {
    throw "fuck";
  }

  virtual ~Renderer()
  {
    for(TextureIdMap::iterator it = _texture_handles.begin(); it != _texture_handles.end(); it++) {
      glDeleteTextures(1, &it->second);
    }
  }

  void bindMaterial(GLenum face, const Material &m)
  {
    glMateriali(face, GL_SHININESS, m.shininess());
    glMaterialrv(face, GL_SPECULAR, m.specular());
    glMaterialrv(face, GL_EMISSION, m.emission());
    glMaterialrv(face, GL_AMBIENT, m.ambient());
    glMaterialrv(face, GL_DIFFUSE, m.diffuse());
  }

  void unbindMaterial(GLenum face)
  {
    glMateriali(face, GL_SHININESS, 127);
    glMaterialrv(face, GL_SPECULAR, Colors::White);
    glMaterialrv(face, GL_EMISSION, Colors::Black);
    glMaterialrv(face, GL_AMBIENT, Colors::Black);
    glMaterialrv(face, GL_DIFFUSE, Colors::White);
  }

  unsigned int findTextureHandle(const Texture &tex)
  {
    return _texture_handles[&tex];
  }

  void bindTexture(const Texture &tex)
  {
    glBindTexture(GL_TEXTURE_2D, findTextureHandle(tex));
  }

  void unbindTexture()
  {
    glBindTexture(GL_TEXTURE_2D, 0);
  }

  bool loadTexture(Texture &tex)
  {
    GLuint tex_id;
    GLenum mode = GL_RGBA;

    if(!tex.loaded()) {
      if(!tex.load()) {
        return false;
      }
    }

    int bpp;
    Uint32 rmask, gmask, bmask, amask;
    SDL_PixelFormatEnumToMasks(SDL_PIXELFORMAT_ABGR8888, &bpp, &rmask, &gmask, &bmask, &amask);
    SDL_Surface *surf = SDL_CreateRGBSurface(0, tex.width(), tex.height(), bpp, rmask, gmask, bmask, amask);
    if(surf == NULL) {
      return false;
    }

    SDL_SetSurfaceAlphaMod(surf, 0xFF);
    SDL_SetSurfaceBlendMode(surf, SDL_BLENDMODE_NONE);
    SDL_BlitSurface(tex.surface(), NULL, surf, NULL);

    switch(surf->format->BytesPerPixel) {
      case 3:
        mode = GL_RGB;
        break;
      case 4:
        mode = GL_RGBA;
        break;
      default:
        std::cerr << "Unknown mode in " << tex.path() << "\t" << (int)surf->format->BytesPerPixel << std::endl;
        return false;
    }

    //glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_ROW_LENGTH, surf->w);
    ASSERT_GL_ERROR;
    glGenTextures(1, &tex_id);
    ASSERT_GL_ERROR;
    glBindTexture(GL_TEXTURE_2D, tex_id);
    ASSERT_GL_ERROR;
    //glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
    //ASSERT_GL_ERROR;
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    ASSERT_GL_ERROR;
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
    ASSERT_GL_ERROR;
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, tex.repeatsX() ? GL_REPEAT : GL_CLAMP_TO_EDGE);
    ASSERT_GL_ERROR;
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, tex.repeatsY() ? GL_REPEAT : GL_CLAMP_TO_EDGE);
    ASSERT_GL_ERROR;
    glTexImage2D(GL_TEXTURE_2D, 0, mode, surf->w, surf->h, 0, mode, GL_UNSIGNED_BYTE, surf->pixels);
    ASSERT_GL_ERROR;
    glGenerateMipmap(GL_TEXTURE_2D);
    ASSERT_GL_ERROR;

    glBindTexture(GL_TEXTURE_2D, 0);

    _texture_handles[&tex] = tex_id;

    SDL_FreeSurface(surf);

    return true;
  }

  void renderCameraAxis(const Camera &camera)
  {
    glBegin(GL_LINES);
    glColor3rv(Colors::Red);
    glVertex3rv(camera.position() + camera.right() * 1024);
    glVertex3rv(camera.position());
    glColor3rv(Colors::Green);
    glVertex3rv(camera.position() + camera.up() * 1024);
    glVertex3rv(camera.position());
    glColor3rv(Colors::Blue);
    glVertex3rv(camera.position() + camera.forward() * 1024);
    glVertex3rv(camera.position());


    glColor3rv(Colors::Yellow);
    glVertex3rv(camera.position() + Vec3(1, 0, 0) * 1024);
    glVertex3rv(camera.position());
    glColor3rv(Colors::Cyan);
    glVertex3rv(camera.position() + Vec3(0, 1, 0) * 1024);
    glVertex3rv(camera.position());
    glColor3rv(Colors::Magenta);
    glVertex3rv(camera.position() + Vec3(0, 0, 1) * 1024);
    glVertex3rv(camera.position());

    glEnd();
  }

  void renderOriginAxis()
  {
    glBegin(GL_LINES);
    glColor3rv(Colors::Red);
    glVertex3r(1024, 0, 0);
    glVertex3r(0, 0, 0);
    glColor3rv(Colors::Green);
    glVertex3r(0, 1024, 0);
    glVertex3r(0, 0, 0);
    glColor3rv(Colors::Blue);
    glVertex3r(0, 0, 1024);
    glVertex3r(0, 0, 0);
    glEnd();
  }

  void renderAxis(const Vec3 &origin, const Vec3 &right = Vec3::X, const Vec3 &up = Vec3::Y)
  {
    glPushMatrix(); {
      glMultMatrixr(Matrix::translation(origin));
      glMultMatrixr(Matrix::lookAt(up.cross(right), up));
      // todo rotate to fit right

      renderOriginAxis();

    } glPopMatrix();
  }

  void renderVector(const Vec3 &origin, const Vec3 &v, const Vec3 &color = Colors::White)
  {
    glLineWidth(16.0);
    glColor3rv(color);
    glBegin(GL_LINES);
    glVertex3rv(origin);
    glVertex3rv(origin + v);
    glEnd();
  }

  void renderNormal(const Vec3 &origin, const Vec3 &normal, const Vec3 &color = Colors::Green)
  {
    renderVector(origin, normal * 512.0, color);
  }

  void renderPole(const Camera &cam, const Vec3 &position, Real value, const Vec3 &color)
  {
    renderVector(position, Vec3(0, 1, 0) * value, color);
  }

  void renderPoint(const Vec3 &position, const Vec3 &color = Colors::White)
  {
    glBegin(GL_POINTS);
    glColor4rv(color);
    glVertex3rv(position);
    glEnd();
  }

  void renderLine(const Vec3 &a, const Vec3 &b, const Vec3 &color = Colors::White)
  {
    glLineWidth(4);
    glBegin(GL_LINES);
    glColor3rv(color);
    glVertex3rv(a);
    glVertex3rv(b);
    glEnd();
  }
};

class Icosahedron
{
public:
  class Face
  {
    Vec3 *_verts[3];

  public:
    Face()
    {
      for(int i = 0; i < 3; i++) {
        _verts[i] = NULL;
      }
    }

    Face(Vec3 *a, Vec3 *b, Vec3 *c)
    {
      _verts[0] = a;
      _verts[1] = b;
      _verts[2] = c;
    }

    const Vec3 &vertex(int i) const
    {
      if(i < 0 || i >= 3) throw IndexOutOfBoundsError(this, i, 3);
      return *_verts[i];
    }

    Vec3 center() const
    {
      return (vertex(0) + vertex(1) + vertex(2)) / 3.0;
    }
  };

  static const int NumPoints = 12;
  static const int NumFaces = 20;

private:

  Vec3 _points[NumPoints];
  Face _faces[NumFaces];
  Real _radius;

public:
  Icosahedron(Real radius = 1.0)
    : _radius(radius)
  {
    generate();
  }

  void generate()
  {
    generatePoints();
    generateFaces();
  }

  void generatePoints()
  {
    Real one = 1.0, psi = (1.0 + sqrt(5.0)) / 2.0;

    _points[0] = Vec3(0, one, psi);
    _points[1] = Vec3(0, -one, psi);
    _points[2] = Vec3(0, one, -psi);
    _points[3] = Vec3(0, -one, -psi);

    _points[4] = Vec3(psi, 0, one);
    _points[5] = Vec3(psi, 0, -one);
    _points[6] = Vec3(-psi, 0, one);
    _points[7] = Vec3(-psi, 0, -one);

    _points[8] = Vec3(one, psi, 0);
    _points[9] = Vec3(-one, psi, 0);
    _points[10] = Vec3(one, -psi, 0);
    _points[11] = Vec3(-one, -psi, 0);

    for(int i = 0; i < NumPoints; i++) {
      _points[i] = _points[i].normalize() * _radius;
    }
  }

  static const int NeighborIndex[NumFaces][3];
  static const int PointNeighborIndex[NumFaces][3][3];
  static const int VertexIndex[NumFaces][3];

  void generateFaces()
  {
    for(int face = 0; face < NumFaces; face++) {
      const int *index = VertexIndex[face];
      _faces[face] = Face(&_points[index[0]], &_points[index[1]], &_points[index[2]]);
    }
  }

  const Vec3 &vertex(int v) const
  {
    if(v < 0 || v >= NumPoints) throw IndexOutOfBoundsError(this, v, NumPoints);
    return _points[v];
  }

  Face &face(int f)
  {
    if(f < 0 || f >= NumFaces) throw IndexOutOfBoundsError(this, f, NumFaces);
    return _faces[f];
  }

  const Face &face(int f) const
  {
    if(f < 0 || f >= NumFaces) throw IndexOutOfBoundsError(this, f, NumFaces);
    return _faces[f];
  }

  const Face &neighbor(int f, int neighbor) const
  {
    if(f < 0 || f >= NumFaces) throw IndexOutOfBoundsError(this, f, NumFaces);
    if(neighbor < 0 || neighbor >= 3) throw IndexOutOfBoundsError(this, neighbor, 3);
    return face(NeighborIndex[f][neighbor]);
  }

  const Face &pointNeighbor(int f, int point, int neighbor) const
  {
    if(f < 0 || f >= NumFaces) throw IndexOutOfBoundsError(this, f, NumFaces);
    if(point < 0 || point >= 3) throw IndexOutOfBoundsError(this, point, 3);
    if(neighbor < 0 || neighbor >= 3) throw IndexOutOfBoundsError(this, neighbor, 3);
    return face(PointNeighborIndex[f][point][neighbor]);
  }

  Real radius() const { return _radius; }
};

const int Icosahedron::NeighborIndex[NumFaces][3] = {
  { 12, 1, 13 },
  { 15, 0, 14 },
  { 17, 3, 16 },
  { 18, 2, 19 },
  { 12, 5, 16 },
  { 17, 4, 13 },
  { 18, 7, 14 },
  { 15, 6, 19 },
  { 12, 9, 14 },
  { 18, 8, 16 },
  { 15, 11, 13 },
  { 17, 10, 19 },
  { 4, 8, 0 },
  { 0, 10, 5 },
  { 1, 8, 6 },
  { 7, 10, 1 },
  { 2, 9, 4 },
  { 5, 11, 2 },
  { 6, 9, 3 },
  { 3, 11, 7 },
};

// todo
const int Icosahedron::PointNeighborIndex[NumFaces][3][3] = {
  { // face 0
    { 0, 4, 5 }, // point 0
    { 0, 8, 14 }, // point 8
    { 0, 10, 15 } // point 9
  },
  { // face 1
    { 1, 6, 7 }, // point 2
    { 1, 10, 13 }, // point 9
    { 1, 8, 12 } // point 8
  },
  { // face 2
    { 2, 4, 5 }, // point 1
    { 2, 11, 19 }, // point 11
    { 2, 9, 18 }, // point 10
  },
  { // face 3
    { 3, 6, 7 }, // point 3
    { 3, 9, 16 }, // point 10
    { 3, 11, 17 }, // point 11
  },
  { // face 4
    { 4, 8, 9 }, // point 4
    { 0, 4, 13 }, // point 0
    { 2, 4, 17 }, // point 1
  },
  { // face 5
    { 5, 10, 11 }, // point 6
    { 2, 5, 16 }, // point 1
    { 0, 5, 12 }, // point 0
  },
  { // face 6
    { 6, 8, 9 }, // point 5
    { 3, 6, 19 }, // point 3
    { 1, 6, 15 }, // point 2
  },
  { // face 7
    { 7, 10, 11 }, // point 7
    { 1, 7, 14 }, // point 2
    { 3, 7, 18 }, // point 3
  },
  { // face 8
    { 0, 1, 8 }, // point 8
    { 4, 8, 16 }, // point 4
    { 6, 8, 18 }, // point 5
  },
  { // face 9
    { 2, 3, 9 }, // point 10
    { 6, 9, 14 }, // point 5
    { 4, 9, 12 }, // point 4
  },
  { // face 10
    { 0, 1, 10 }, // point 9
    { 7, 10, 19 }, // point 7
    { 5, 10, 17 }, // point 6
  },
  { // face 11
    { 2, 3, 11 }, // point 11
    { 5, 11, 13 }, // point 6
    { 7, 11, 15 }, // point 7
  },
  { // face 12
    { 5, 12, 13 }, // point 0
    { 9, 12, 16 }, // point 4
    { 1, 12, 14 }, // point 8
  },
  { // face 13
    { 4, 12, 13 }, // point 0
    { 1, 13, 15 }, // point 9
    { 11, 13, 17 }, // point 6
  },
  { // face 14
    { 7, 14, 15 }, // point 2
    { 0, 12, 14 }, // point 8
    { 9, 14, 18 }, // point 5
  },
  { // face 15
    { 6, 14, 15 }, // point 2
    { 11, 15, 19 }, // point 7
    { 0, 13, 15 }, // point 9
  },
  { // face 16
    { 5, 16, 17 }, // point 1
    { 3, 16, 18 }, // point 10
    { 8, 12, 16 }, // point 4
  },
  { // face 17
    { 4, 16, 17 }, // point 1
    { 10, 13, 17 }, // point 6
    { 3, 17, 19 }, // point 11
  },
  { // face 18
    { 7, 18, 19 }, // point 3
    { 8, 14, 18 }, // point 5
    { 2, 16, 18 }, // point 10
  },
  { // face 19
    { 6, 18, 19 }, // point 3
    { 2, 17, 19 }, // point 11
    { 10, 15, 19 }, // point 7
  }
};

const int Icosahedron::VertexIndex[NumFaces][3] = {
  { 0, 8, 9 },
  { 2, 9, 8 },
  { 1, 11, 10 },
  { 3, 10, 11 },

  { 4, 0, 1 },
  { 6, 1, 0 },
  { 5, 3, 2 },
  { 7, 2, 3 },

  { 8, 4, 5 },
  { 10, 5, 4 },
  { 9, 7, 6 },
  { 11, 6, 7 },

  { 0, 4, 8 },
  { 0, 9, 6 },
  { 2, 8, 5 },
  { 2, 7, 9 },

  { 1, 10, 4 },
  { 1, 6, 11 },
  { 3, 5, 10 },
  { 3, 11, 7 },
};

class IcosahedronRenderer
{
  Renderer &_renderer;

public:
  IcosahedronRenderer(Renderer &r)
  : _renderer(r)
  {
  }

  void render(const Camera &cam, const Icosahedron &iso)
  {
    glBegin(GL_TRIANGLES);

    for(int i = 0; i < Icosahedron::NumFaces; i++) {
      const Icosahedron::Face &face = iso.face(i);
      glNormal3rv(face.vertex(0).normalize());
      glVertex3rv(face.vertex(0));
      glNormal3rv(face.vertex(1).normalize());
      glVertex3rv(face.vertex(1));
      glNormal3rv(face.vertex(2).normalize());
      glVertex3rv(face.vertex(2));
    }

    glEnd();
  }

  void renderNeighbors(const Camera &cam, const Icosahedron &iso)
  {
    glPushMatrix(); {
      glScalef(1.01, 1.01, 1.01);

      glBegin(GL_LINES);
      for(int face = 0; face < Icosahedron::NumFaces; face++) {
        for(int n = 0; n < 3; n++) {
          glColor3rv(Colors::Black);
          glVertex3rv(iso.face(face).center());
          glColor3rv(Colors::color(n));
          glVertex3rv(iso.neighbor(face, n).center());
        }
      }
      glEnd();
    } glPopMatrix();
  }

  void renderPointNeighbors(const Camera &cam, const Icosahedron &iso, int point)
  {
    glPushMatrix(); {
      glScalef(1.01, 1.01, 1.01);

      glBegin(GL_LINES);
      for(int face = 0; face < Icosahedron::NumFaces; face++) {
        for(int n = 0; n < 3; n++) {
          glColor3rv(Colors::Black);
          glVertex3rv(iso.face(face).vertex(point));
          glColor3rv(Colors::color(n));
          glVertex3rv(iso.pointNeighbor(face, point, n).center());
        }
      }
      glEnd();
    } glPopMatrix();
  }
};

class PlayerCommandState
{
  bool _forward, _backward, _left, _right, _up, _down;

public:
  PlayerCommandState()
    : _forward(false), _backward(false), _left(false), _right(false), _up(false), _down(false)
  {
  }

  void movingForward(bool yes) { _forward = yes; }
  bool isMovingForward() const { return _forward; }

  void movingBackward(bool yes) { _backward = yes; }
  bool isMovingBackward() const { return _backward; }

  void movingLeft(bool yes) { _left = yes; }
  bool isMovingLeft() const { return _left; }

  void movingRight(bool yes) { _right = yes; }
  bool isMovingRight() const { return _right; }

  void movingUp(bool yes) { _up = yes; }
  bool isMovingUp() const { return _up; }

  void movingDown(bool yes) { _down = yes; }
  bool isMovingDown() const { return _down; }

  Vec3 motionVector()
  {
    Vec3 dir(0, 0, 0, 0);

    if(isMovingForward())
      dir += Vec3(0.0, 0.0, -1.0);
    if(isMovingBackward())
      dir += Vec3(0.0, 0.0, 1.0);
    if(isMovingLeft())
      dir += Vec3(-1.0, 0.0, 0.0);
    if(isMovingRight())
      dir += Vec3(1.0, 0.0, 0.0);
    if(isMovingUp())
      dir += Vec3(0.0, 1.0, 0.0);
    if(isMovingDown())
      dir += Vec3(0.0, -1.0, 0.0);

    return dir.normalize();
  }
};

namespace XPad
{
  enum Axis {
    AXIS_LEFT_X,
    AXIS_LEFT_Y,
    AXIS_LEFT_TRIGGER,
    AXIS_RIGHT_X,
    AXIS_RIGHT_Y,
    AXIS_RIGHT_TRIGGER
  };

  const int AXIS_MIN = -32768, AXIS_MAX = 32768, AXIS_RANGE = AXIS_MAX - AXIS_MIN;
  const int AXIS_LEFT_DEAD_ZONE = 10000, AXIS_RIGHT_DEAD_ZONE = 4000;

  enum Button {
    BUTTON_A,
    BUTTON_B,
    BUTTON_BLACK,
    BUTTON_X,
    BUTTON_Y,
    BUTTON_WHITE,
    BUTTON_BACK,
    BUTTON_START
  };
};

class Quantum
{
  Vec3 _position, _velocity, _color;
  Real _scale, _mass;

public:
  Quantum(const Vec3 &position, const Vec3 &velocity, Real scale = 100, const Vec3 &color = Colors::White, Real mass = 1.0)
    : _position(position), _velocity(velocity), _color(color), _scale(scale), _mass(mass)
  {
  }

  const Vec3 &position() const { return _position; }
  Quantum &setPosition(const Vec3 &v) { _position = v; return *this; }
  const Vec3 &velocity() const { return _velocity; }
  Quantum &setVelocity(const Vec3 &v) { _velocity = v; return *this; }
  const Vec3 &color() const { return _color; }
  Quantum &setColor(const Vec3 &c) { _color = c; return *this; }
  Real mass() const { return _mass; }
  Quantum &setMass(Real m) { _mass = m; return *this; }

  Real universeSize() const { return _scale; }

  virtual void update(Real dt)
  {
    Vec3 r = (_position + _velocity * dt);

#if defined(WRAPAROUND_CUBE)
    if(r.x() > _scale) {
      r.setX(-_scale);
    } else if(r.x() < -_scale) {
      r.setX(_scale);
    }

    if(r.y() > _scale) {
      r.setY(-_scale);
    } else if(r.y() < -_scale) {
      r.setY(_scale);
    }

    if(r.z() > _scale) {
      r.setZ(-_scale);
    } else if(r.z() < -_scale) {
      r.setZ(_scale);
    }
#elif defined(WRAPAROUND_SPHERE)
    if(r.magnitude() >= _scale) {
      r = r.normalize() * -_scale;
    }
#endif

    _position = r;
  }
};

class Quanta
{
  typedef std::vector<Quantum *> vector;
  vector _quanta;
  Real _quantum_radius;
  bool _swap_color;

public:
  Quanta(Real quantum_radius)
    : _quantum_radius(quantum_radius), _swap_color(false)
  {
  }

  virtual ~Quanta()
  {
    for(int i = 0; i < _quanta.size(); i++) {
      delete _quanta[i];
    }
  }

  void push_back(Quantum *q)
  {
    _quanta.push_back(q);
  }

  Quantum &operator[] (int n)
  {
    return *_quanta[n];
  }

  const Quantum &operator[] (int n) const
  {
    return *_quanta[n];
  }

  int size() const { return _quanta.size(); }
  Real quantum_radius() const { return _quantum_radius; }

  bool swapColor() const {
    return _swap_color;
  }

  bool swapColor(bool yes) {
    return _swap_color = yes;
  }

  virtual void update(Real dt)
  {
    typedef std::pair<Real, int> collision_t;
    struct collision_compare: public std::less<collision_t>
    {
       bool operator() (const collision_t &a, const collision_t &b) {
         return a.first < b.first;
       }
    };
    typedef std::set<collision_t, collision_compare> seen_list_t;
    seen_list_t colliders;

    bool *seen = new bool[size()];
    for(int i = 0; i < size(); i++) {
      seen[i] = false;
    }

    // test for collisions
    for(int i = 0; i < size(); i++) {
      // find any colliders
      colliders.clear();

      for(int j = 0; j < size(); j++) {
        if(i != j && seen[i] != true && seen[j] != true) {
          Real t = collidesAt((*this)[i], (*this)[j], dt);
          if(!isnan(t)) {
            colliders.insert(collision_t(t, j));
          }
        }
      }

      // update position and velocity
      Real t_remaining = dt;

      if(colliders.size() > 0) {
        for(auto c: colliders) { // todo verify positioning and timing of 2nd and subsequent collision
          Real t = collidesAt((*this)[i], (*this)[c.second], t_remaining);
          if(!isnan(t)) {
            collide((*this)[i], (*this)[c.second], t, t_remaining);
            t_remaining -= t;
          }
          //seen[c.second] = true;
        }
      }

      (*this)[i].update(t_remaining);
      seen[i] = true;
    }

    delete[] seen;
  }

  Real collidesAt(const Quantum &a, const Quantum &b, Real dt)
  {
    Vec3 dp = (a.position() - b.position());
    Real dist = dp.magnitude();
    Vec3 dv = (a.velocity() - b.velocity());
    Real rel_speed = dv.magnitude();

    // quick escapes
    if((dp / dist).dot(dv / rel_speed) > 0.0 ||
       (dist - 2.0*_quantum_radius) > rel_speed) {
      //std::cerr << "No collision between " << a.position() << ":" << a.velocity() << " and " << b.position() << ":" << b.velocity() << "? " << (dp / dist).dot(dv / rel_speed) << "\t" << (a.position() - b.position()).magnitude() << std::endl;
      return NAN;
    }

    // Thanks Qalculate!
    // av^2 t^2 + 2*av*ax*t - 2*av*bv*t^2 - 2 * av*bx*t + ax^2 - 2*ax*bv*t - 2*ax*bx + bv^2 * t^2 + 2*bv *bx*t + bx^2
    // av^2 t^2 - 2*av*bv*t^2 + bv^2 * t^2 - 2 * av*bx*t - 2*ax*bv*t + 2*bv *bx*t + 2*av*ax*t + ax^2 + bx^2  - 2*ax*bx
    // t^2(av^2 - 2*av*bv + bv^2) + t(-2 * av*bx* - 2*ax*bv* + 2*bv *bx* + 2*av*ax) + (ax^2 + bx^2  - 2*ax*bx)

    Vec3 iv = a.position() * a.position() - 2.0 * a.position() * b.position() + b.position() * b.position();
    Vec3 jv = 2.0*a.position() * a.velocity() - 2.0 * a.position() * b.velocity() - 2.0 * b.position() * a.velocity() + 2.0 * b.position() * b.velocity();
    Vec3 kv = a.velocity() * a.velocity() - 2.0 * a.velocity() * b.velocity() + b.velocity() * b.velocity();

    // All that was for a distance squared formula, so add the components up.
    Real i = -4.0 * _quantum_radius * _quantum_radius + iv.x() + iv.y() + iv.z();
    Real j = jv.x() + jv.y() + jv.z();
    Real k = kv.x() + kv.y() + kv.z();

    // quadratic equation
    Real bac = j*j - 4.0*i*k;
    if(bac < 0.0) {
      return NAN;
    }

    Real tp = (-j + sqrt(bac)) / (2.0*k);
    Real tm = (-j - sqrt(bac)) / (2.0*k);

    // pick the t
    if((isnan(tp) || tp < 0.0 || tp > dt) && (isnan(tm) || tm < 0.0 || tm > dt)) {
      //std::cerr << "No collision between " << a.position() << ":" << a.velocity() << " and " << b.position() << ":" << b.velocity() << "? " << tp << " " << tm << minimum(tp, tm) << "\t" << (a.position() - b.position()).magnitude() << std::endl;
      return NAN;
    }

    //std::cerr << a.position() << ":" << a.velocity() << " and " << b.position() << ":" << b.velocity() << " " << "\t" << minimum(tp, tm) << "\t" << (a.position() - b.position()).magnitude() << std::endl;

    return minimum(tp, tm);
  }

  bool collide(Quantum &a, Quantum &b, Real t, Real dt)
  {
    Vec3 ac = a.position() + a.velocity() * t;
    Vec3 bc = b.position() + b.velocity() * t;

#ifdef BADMATH
    Vec3 n = ac - bc;
    Vec3 na = a.velocity().projectOnto(n);
    Vec3 nb = b.velocity().projectOnto(n);
    Vec3 va = (a.velocity() - na) + nb;
    Vec3 vb = (b.velocity() - nb) + na;
#else
    // courtesy of http://www.gamasutra.com/view/feature/131424/pool_hall_lessons_fast_accurate_.php
    Vec3 n = (ac - bc).normalize();
    Real aa = a.velocity().dot(n);
    Real ab = b.velocity().dot(n);
    Real p = (2.0 * (aa - ab)) / (a.mass() + b.mass());
    Vec3 va = a.velocity() - p * b.mass() * n;
    Vec3 vb = b.velocity() + p * a.mass() * n;
#endif

    if(va.isNaN() || vb.isNaN()) return false;

    warn(((a.mass()*a.velocity() + b.mass()*b.velocity()) - (a.mass()*va + b.mass()*vb)).magnitude() <= 0.01);

    a.setVelocity(va);
    a.setPosition(ac);
    b.setVelocity(vb);
    b.setPosition(bc);

    /*
    Vec3 ca = b.color().reflectBy(n).normalize();
    Vec3 cb = a.color().reflectBy(n).normalize();
    a.setColor(ca);
    b.setColor(cb);
    */

    if(_swap_color) {
      Vec3 c = a.color();
      a.setColor(b.color());
      b.setColor(c);
    }

    return true;
  }
};

class QuantaRenderer
{
public:
  Icosahedron _icosa;
  IcosahedronRenderer _icosa_ren;
  Renderer &_renderer;

  QuantaRenderer(Renderer &renderer)
    : _icosa(1.0), _icosa_ren(renderer), _renderer(renderer)
  {
  }

  void render(const Camera &camera, const Quantum &quantum, Real radius, Real alpha = 1.0)
  {
     glPushMatrix();
     glTranslater(quantum.position().x(), quantum.position().y(), quantum.position().z());
     glScaler(radius, radius, radius);
     glColor3rv(quantum.color().abs() * Vec3(1.0, 1.0, 1.0, alpha));
     _icosa_ren.render(camera, _icosa);
     glPopMatrix();
  }

  void render(const Camera &camera, const Quanta &quanta)
  {
    for(int i = 0; i < quanta.size(); i++) {
      Quantum q = quanta[i];
      glPushMatrix(); {
        glTranslater(q.position().x(), q.position().y(), q.position().z());
        glScaler(quanta.quantum_radius(), quanta.quantum_radius(), quanta.quantum_radius());
        glColor3rv(q.color().abs());
        _icosa_ren.render(camera, _icosa);
      } glPopMatrix();
    }
  }

  void renderVelocities(const Camera &camera, const Quanta &quanta)
  {
    for(int i = 0; i < quanta.size(); i++) {
      Quantum q = quanta[i];
      _renderer.renderVector(q.position(), q.velocity() * _icosa.radius(), q.color());
    }
  }
};

class ScreenShooter
{
  int _window_width, _window_height, _bpp, _shots;
  unsigned char *_buffer;
  pid_t _pid;

  enum Mode
  {
    OFF,
    SINGLE,
    MOVIE,
    UNKNOWN
  };

  Mode _mode;

public:
  class Error
  {
    public:
      const ScreenShooter &shooter;

      Error(const ScreenShooter &s)
        : shooter(s)
      {
      }
  };

  ScreenShooter(int w, int h)
    : _buffer(NULL), _bpp(3), _shots(0), _mode(OFF)
  {
    _pid = getpid();
    setWindowSize(w, h);
  }

  ~ScreenShooter()
  {
    if(_buffer != NULL) delete[] _buffer;
  }

  void setWindowSize(int w, int h)
  {
    _window_width = w;
    _window_height = h;
    if(_buffer != NULL) {
      delete[] _buffer;
    }
    _buffer = new unsigned char[_bpp * w * h];
  }

  void prepForShot()
  {
    _mode = SINGLE;
  }

  void prepForMovie()
  {
    _mode = MOVIE;
  }

  void toggleMovie()
  {
    if(_mode == MOVIE) {
      _mode = OFF;
    } else {
      _mode = MOVIE;
    }
  }

  void takeScreenShot()
  {
    switch(_mode) {
      case SINGLE:
        _mode = OFF;
        // move to next case
      case MOVIE:
        doScreenShot();
      default:
        break;
    }
  }

private:
  void doScreenShot()
  {
    // read the pixels to a buffer
    glReadPixels(0, 0, _window_width, _window_height, GL_RGB, GL_UNSIGNED_BYTE, _buffer);
    ASSERT_GL_ERROR;
    mirrorBuffer();

    // convert the buffer to an image format
    SDL_Surface *surf = SDL_CreateRGBSurfaceFrom(_buffer, _window_width, _window_height, _bpp*8, _window_width * _bpp, 0xFF, 0xFF00, 0xFF0000, 0);
    if(surf == NULL) {
      throw Error(*this);
    }

    auto path = nextShotPath();
    if(SDL_SaveBMP(surf, path.c_str()) == 0) {
      std::cerr << "Screenshot " << path << " saved." << std::endl;
    } else {
      std::cerr << "Error saving screenshot " << path << std::endl;
      throw Error(*this);
    }

    SDL_FreeSurface(surf);

    _shots++;
  }

  void mirrorBuffer()
  {
    const int pitch = _window_width * _bpp;
    unsigned char _line[pitch];

    for(int y = 0; y < _window_height / 2; y++ ) {
      unsigned char *top = &_buffer[y * pitch], *bottom = &_buffer[(_window_height - y) * pitch];
      memcpy(_line, top, pitch);
      memcpy(top, bottom, pitch);
      memcpy(bottom, _line, pitch);
    }
  }

  std::string nextShotPath()
  {
    // compute a path to save a single shot $PWD/quanta-$PID-$FRAME.(png|bmp)
    return std::string("quanta-" + std::to_string((unsigned int)_pid) + "-" + std::to_string(_shots) + ".bmp");
  }
};

int window_width = APP_WIDTH, window_height = APP_HEIGHT;

using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  unsigned int frames = 0, ms_per_frame = 0, frame_start = 0;
  int joystick_index = 0;
  Real near_plane = 4.0, far_plane = 1024.0 * 8.0;
  int num_quanta = 125, universe_size = 2048, generation = 0;
  Real quantum_speed = 0.0, quantum_forward_speed = 0.0, quantum_radius = 0.5f, quantum_mass = 1.0;
  Real mp = quantum_radius * 2.5f;
  Vec3 min_position(-mp, -mp, -mp), max_position(mp, mp, mp);

  if(argc > 1) {
    num_quanta = atoi(argv[1]);
  }
  if(argc > 2) {
    Real n = ator(argv[2]);
    min_position = Vec3(-n, -n, -n);
    max_position = Vec3(n, n, n);
  }
  if(argc > 3) {
    universe_size = ator(argv[3]);;
  }
  if(argc > 4) {
    Real n = ator(argv[4]);
    quantum_speed = n;
  }
  if(argc > 5) {
    quantum_radius = ator(argv[5]);
  }
  if(argc > 6) {
    quantum_mass = ator(argv[6]);
  }

  if(SDL_Init(SDL_INIT_VIDEO|SDL_INIT_JOYSTICK) != 0) {
    cout << "SDL Init Error: " << SDL_GetError() << endl;
    return -1;
  }

  //SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

  int num_displays = SDL_GetNumVideoDisplays();
  int display = 1;
  cout << "Using display " << display << "/" << num_displays << endl;

  SDL_Window *window = SDL_CreateWindow(APP_TITLE, SDL_WINDOWPOS_CENTERED_DISPLAY(display), SDL_WINDOWPOS_CENTERED_DISPLAY(display), APP_WIDTH, APP_HEIGHT, SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
  if(window == NULL) {
    cout << "Error creating window: " << SDL_GetError() << endl;
    return -1;
  }

  SDL_GLContext gl_context = SDL_GL_CreateContext(window);
  if(!gl_context) {
    cout << "Error creating gl context: " << SDL_GetError() << endl;
    SDL_DestroyWindow(window);
    return -1;
  }

  Renderer renderer;

  /*
  Texture logo(APP_LOGO);
  if(!renderer.loadTexture(logo)) {
    std::cout << "Error loading " << APP_LOGO << ": " << SDL_GetError() << endl;
    SDL_DestroyWindow(window);
    SDL_Quit();
    return -1;
  }
  */

  SDL_Joystick *joystick = SDL_JoystickOpen(joystick_index);
  if(joystick == NULL) {
    std::cout << "Error opening joystick " << joystick_index << ": " << SDL_GetError() << endl;
    std::cout << "Try one of the following: " << std::endl;
    for(int i = 0; i < SDL_NumJoysticks(); i++) {
      std::cout << "\t" << SDL_JoystickNameForIndex(i) << std::endl;
    }
  } else {
    std::cout << "Joystick " << SDL_JoystickName(joystick) << " opened:" << std::endl;
    std::cout << "\tAxes: " << SDL_JoystickNumAxes(joystick) << std::endl;
    std::cout << "\tBalls: " << SDL_JoystickNumBalls(joystick) << std::endl;
    std::cout << "\tHats: " << SDL_JoystickNumHats(joystick) << std::endl;
    std::cout << "\tButtons: " << SDL_JoystickNumButtons(joystick) << std::endl;
  }

  glMatrixMode(GL_TEXTURE);
  glLoadIdentity();

  /*
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, 0.0, 1.0, 1.0, 0.0, 100.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glEnable(GL_TEXTURE_2D);
  Real w = 1.0, h = 1.0;
  renderer.bindTexture(logo);

  glBegin(GL_QUADS);
  glTexCoord2f(0.0, 0.0);
  glVertex2f(-1.0, -1.0);
  glTexCoord2f(0.0, h);
  glVertex2f(-1.0, 1.0);
  glTexCoord2f(w, h);
  glVertex2f(1.0, 1.0);
  glTexCoord2f(w, 0.0);
  glVertex2f(1.0, -1.0);
  glEnd();

  renderer.unbindTexture();
  SDL_GL_SwapWindow(window);
  */

  Real spawn_height = 0.0; //-3016.0;
  FPSCamera fps_camera(Vec3(0, spawn_height, 0), Vec3());
  Camera *camera = &fps_camera;
  Camera *view_camera = camera;

  Quanta quanta(quantum_radius);
  QuantaRenderer quanta_renderer(renderer);

  PlayerCommandState player_command_state;

  SDL_Event event;
  bool done = false;

  glMatrixMode(GL_TEXTURE);
  glLoadIdentity();

  glMatrixMode(GL_COLOR);
  glLoadIdentity();

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);

  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  /*
  glEnable(GL_FOG);
  glFogrv(GL_FOG_COLOR, Colors::Black);
  glFogi(GL_FOG_MODE, GL_LINEAR);
  glFogi(GL_FOG_START, near_plane);
  glFogi(GL_FOG_END, far_plane / 2.0);
  */
  
  glEnable(GL_LIGHTING);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_COLOR_MATERIAL);

  //glLightModelrv(GL_LIGHT_MODEL_AMBIENT, s.sky());
  //glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

  Vec3 light_color(1.0, 1.0, 1.0);
  glLightr(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0.1); //brightness);
  glLightrv(GL_LIGHT0, GL_DIFFUSE, light_color);
  glLightrv(GL_LIGHT0, GL_SPECULAR, Colors::White);
  glLightrv(GL_LIGHT0, GL_AMBIENT, light_color * 0.01);

  glEnable(GL_LIGHT0);

  Real speed = 16.0;
  Real droll = 0.0, dpitch = 0.0, dyaw = 0.0;

  //Vec3 joy_sensitivity(50, 50, 150);
  Vec3 joy_sensitivity(0.8, 0.8, 0.25);
  Real joy_pitch = NAN, joy_yaw = NAN, joy_roll = NAN;

  ScreenShooter screen_shooter(window_width, window_height);

  while(!done) {
    frame_start = SDL_GetTicks();
    droll = dpitch = dyaw = 0.0;

    while(SDL_PollEvent(&event)) {
      switch(event.type) {
      case SDL_QUIT:
        done = true;
        break;
      case SDL_KEYDOWN:
      case SDL_KEYUP:
        switch(event.key.keysym.scancode) {
        case SDL_SCANCODE_ESCAPE:
          done = true;
          break;
        case SDL_SCANCODE_F7:
          if(event.type == SDL_KEYDOWN) {
            if(SDL_GetModState() & KMOD_ALT) {
              screen_shooter.toggleMovie();
            } else {
              screen_shooter.prepForShot();
            }
          }
          break;
        case SDL_SCANCODE_UP:
        case SDL_SCANCODE_W:
          player_command_state.movingForward(event.type == SDL_KEYDOWN);
          break;
        case SDL_SCANCODE_DOWN:
        case SDL_SCANCODE_S:
          player_command_state.movingBackward(event.type == SDL_KEYDOWN);
          break;
        case SDL_SCANCODE_LEFT:
        case SDL_SCANCODE_A:
          player_command_state.movingLeft(event.type == SDL_KEYDOWN);
          break;
        case SDL_SCANCODE_RIGHT:
        case SDL_SCANCODE_D:
          player_command_state.movingRight(event.type == SDL_KEYDOWN);
          break;
        case SDL_SCANCODE_SPACE:
          player_command_state.movingUp(event.type == SDL_KEYDOWN);
          break;
        case SDL_SCANCODE_LSHIFT:
          player_command_state.movingDown(event.type == SDL_KEYDOWN);
          break;
        case SDL_SCANCODE_Q:
        case SDL_SCANCODE_PAGEDOWN:
          droll = 1.0 * 2.0 * M_PI / 360.0;
          break;
        case SDL_SCANCODE_E:
        case SDL_SCANCODE_PAGEUP:
          droll = -(1.0 * 2.0 * M_PI / 360.0);
          break;
        case SDL_SCANCODE_R:
          if(event.type == SDL_KEYDOWN) {
            speed *= 2;
            std::cout << "Speed: " << speed << std::endl;
          }
          break;
        case SDL_SCANCODE_F:
          if(event.type == SDL_KEYDOWN) {
            speed *= 0.5;
            std::cout << "Speed: " << speed << std::endl;
          }
          break;
        case SDL_SCANCODE_F1:
          std::cout << "FPS: " << frames << " " << SDL_GetTicks() / 1000.0 << "\t" << frames / ((Real)SDL_GetTicks() / 1000.0) << "\t" << 1.0 / (ms_per_frame / 1000.0) << std::endl;
          break;
        case SDL_SCANCODE_F3:
          break;
        case SDL_SCANCODE_F4:
          break;
        case SDL_SCANCODE_F5:
          if(event.type == SDL_KEYDOWN) {
            if(view_camera == camera) {
              view_camera = new Camera(*camera);
              std::cout << "Camera decoupled" << std::endl;
            } else {
              delete view_camera;
              view_camera = camera;
              std::cout << "Camera coupled" << std::endl;
            }
          }
          break;
        case SDL_SCANCODE_F6:
          if(event.type == SDL_KEYDOWN) {
            if(quanta.swapColor(!quanta.swapColor())) {
              std::cout << "Swapping colors on collision." << std::endl;
            } else {
              std::cout << "Not swapping colors on collision." << std::endl;
            }
          }
          break;
        case SDL_SCANCODE_B:
          camera->moveBy(Vec3(0, 0, -128));
          break;
        case SDL_SCANCODE_0:
          if(event.type == SDL_KEYDOWN) {
            std::cout << "Reoriented the camera" << std::endl;
            camera->lookAt(-Vec3::Z, Vec3::Y);
          }
          break;
        case SDL_SCANCODE_G:
          if(camera == &fps_camera && event.type == SDL_KEYDOWN) {
            std::cout << "Flipped camera upside down." << std::endl;
            fps_camera.setPlane(Matrix::rotateZ(M_PI) * fps_camera.plane());
          } break;
        case SDL_SCANCODE_V:
          if(event.type == SDL_KEYDOWN) {
            if(camera == &fps_camera) {
              Camera *c = new Camera(*camera);
              camera = c;
            } else {
              fps_camera.setPosition(camera->position());
              fps_camera.lookAt(-Vec3::Z, Vec3::Y);
              delete camera;
              camera = &fps_camera;
            }
            view_camera = camera;
            std::cout << "Toggled flight" << std::endl;
          } break;
        case SDL_SCANCODE_F2:
          if(event.type == SDL_KEYDOWN) {
            std::cout << "T:\t" << SDL_GetTicks() << std::endl;
            std::cout << "Camera:" << std::endl;
            std::cout << "  Matrix: " << camera->matrix() << std::endl;
            std::cout << "  Position: " << camera->position() << std::endl;
            std::cout << "  Yaw: " << camera->yaw() << std::endl;
            std::cout << "  Right: " << camera->right() << "\tUp: " << camera->up() << "\tFront: " << camera->forward() << std::endl;
            std::cout << "  Rotation: " << camera->rotationMatrix() << camera->rotationMatrix() * camera->rotationMatrix().invert() << std::endl;
            std::cout << "  Translation: " << camera->translationMatrix() << std::endl;
          }
          break;
        case SDL_SCANCODE_RETURN:
          if(event.type == SDL_KEYDOWN) {
            quanta.push_back(new Quantum(camera->forward() * near_plane * 1.5f + camera->position(), camera->forward() * quantum_speed, universe_size, Colors::color(), quantum_mass));
          }
          break;
        case SDL_SCANCODE_EQUALS:
          if(event.type == SDL_KEYDOWN) {
            quantum_speed += 0.1;
            std::cout << "Quantum speed: " << quantum_speed << std::endl;
          }
          break;
        case SDL_SCANCODE_MINUS:
          if(event.type == SDL_KEYDOWN) {
            quantum_speed -= 0.1;
            std::cout << "Quantum speed: " << quantum_speed << std::endl;
          }
          break;
        case SDL_SCANCODE_U:
          if(event.type == SDL_KEYDOWN) {
            quantum_mass -= 1.0;
            std::cout << "Quantum mass: " << quantum_mass << std::endl;
          }
          break;
        case SDL_SCANCODE_I:
          if(event.type == SDL_KEYDOWN) {
            quantum_mass += 1.0;
            std::cout << "Quantum mass: " << quantum_mass << std::endl;
          }
          break;
        case SDL_SCANCODE_LEFTBRACKET:
          if(event.type == SDL_KEYDOWN) {
            num_quanta -= 1;
            std::cout << "Number of quanta decreased to " << num_quanta << std::endl;
          }
          break;
        case SDL_SCANCODE_RIGHTBRACKET:
          if(event.type == SDL_KEYDOWN) {
            num_quanta += 1;
            std::cout << "Number of quanta increased to " << num_quanta << std::endl;
          }
          break;
        case SDL_SCANCODE_P:
          if(event.type == SDL_KEYDOWN) {
            quantum_forward_speed += 0.1;
            std::cout << "Quantum forward speed: " << quantum_forward_speed << std::endl;
          }
          break;
        case SDL_SCANCODE_O:
          if(event.type == SDL_KEYDOWN) {
            quantum_forward_speed -= 0.1;
            std::cout << "Quantum forward speed: " << quantum_forward_speed << std::endl;
          }
          break;
        case SDL_SCANCODE_APOSTROPHE:
          if(event.type == SDL_KEYDOWN) {
            Vec3 color = Colors::color(generation++);
            for(int i = 0; i < num_quanta; i++) {
              quanta.push_back(new Quantum(camera->forward() * near_plane * 1.5f + camera->position() + Vec3::random(min_position, max_position),
                                       quantum_forward_speed * camera->forward() + quantum_speed * Vec3::random(),
                                       universe_size,
                                       color,
                                       quantum_mass));
            }
          }
          break;
        case SDL_SCANCODE_L:
          if(event.type == SDL_KEYDOWN) {
            Vec3 color = Colors::color(generation++);
            Real d = (max_position - min_position).magnitude();
            Real maxq = powf(num_quanta, 1.0/3.0);
            for(int x = 0; x < maxq; x++) {
              for(int y = 0; y < maxq; y++) {
                for(int z = 0; z < maxq; z++) {
                  quanta.push_back(new Quantum(camera->forward() * near_plane * 1.5f + camera->position() + (Vec3(x, y, z) / maxq * d - d / 2),
                                               quantum_forward_speed * camera->forward() + quantum_speed * Vec3::random(),
                                               universe_size,
                                               color,
                                               quantum_mass));
                }
              }
            }
          }
          break;
        }
        break;
      case SDL_WINDOWEVENT:
        switch(event.window.event) {
        case SDL_WINDOWEVENT_RESIZED:
          window_width = event.window.data1;
          window_height = event.window.data2;
          glViewport(0, 0, window_width, window_height);
          screen_shooter.setWindowSize(window_width, window_height);
          break;
        }
        break;
      case SDL_MOUSEMOTION:
#ifdef MAC
        if(event.motion.state & SDL_BUTTON_LMASK) {
#else
        if(event.motion.state & SDL_BUTTON_RMASK) {
#endif /* MAC */
          dyaw -= event.motion.xrel * 2.0 * M_PI / 360.0;
          dpitch -= event.motion.yrel * 2.0 * M_PI / 360.0;
          //std::cout << "MouseMotion: " << event.motion.xrel << "\t" << event.motion.yrel << endl;
        }
        break;
      case SDL_MOUSEBUTTONDOWN:
        //SDL_SetWindowGrab();
        switch(event.button.button) {
          case SDL_BUTTON_LEFT:
            if(event.button.clicks == 1) {
              quanta.push_back(new Quantum(camera->forward() * near_plane * 1.5f + camera->position(), camera->forward() * quantum_speed, universe_size, Colors::color(), quantum_mass));
            }
            break;
          case SDL_BUTTON_RIGHT:
            SDL_SetRelativeMouseMode((SDL_bool)true);
            break;
          default:
            std::cerr << "Unknown button: " << event.button.button << std::endl;
            break;
         }
        break;
      case SDL_MOUSEBUTTONUP:
        switch(event.button.button) {
          case SDL_BUTTON_RIGHT:
            SDL_SetRelativeMouseMode((SDL_bool)false);
            break;
        }
        break;
      case SDL_MOUSEWHEEL:
        if(event.wheel.y > 0) {
          quantum_speed += 0.1;
          std::cout << "Quantum speed: " << quantum_speed << std::endl;
        } else if(event.wheel.y < 0) {
          quantum_speed -= 0.1;
          std::cout << "Quantum speed: " << quantum_speed << std::endl;
        }
        if(event.wheel.x > 0) {
          quantum_mass += 1.0;
          std::cout << "Quantum mass: " << quantum_mass << std::endl;
        } else if(event.wheel.x < 0) {
          quantum_mass -= 1.0;
          std::cout << "Quantum mass: " << quantum_mass << std::endl;
        }
        break;
      case SDL_JOYAXISMOTION:
        //std::cout << "JoyAxisMotion: " << event.jaxis.which << "\t" << (int)event.jaxis.axis << "\t" << event.jaxis.value << std::endl;
        switch(event.jaxis.axis) {
        case XPad::AXIS_LEFT_X:
          if(event.jaxis.value > XPad::AXIS_LEFT_DEAD_ZONE) {
            player_command_state.movingLeft(false);
            player_command_state.movingRight(true);
          } else if(event.jaxis.value < -XPad::AXIS_LEFT_DEAD_ZONE) {
            player_command_state.movingLeft(true);
            player_command_state.movingRight(false);
          } else {
            player_command_state.movingLeft(false);
            player_command_state.movingRight(false);
          }
          break;
        case XPad::AXIS_LEFT_Y:
          if(event.jaxis.value > XPad::AXIS_LEFT_DEAD_ZONE) {
            player_command_state.movingForward(false);
            player_command_state.movingBackward(true);
          } else if(event.jaxis.value < -XPad::AXIS_LEFT_DEAD_ZONE) {
            player_command_state.movingForward(true);
            player_command_state.movingBackward(false);
          } else {
            player_command_state.movingForward(false);
            player_command_state.movingBackward(false);
          }
          break;
        case XPad::AXIS_RIGHT_Y:
          if(rabs(event.jaxis.value) > XPad::AXIS_RIGHT_DEAD_ZONE) {
            joy_pitch = event.jaxis.value * 2.0 * M_PI / (Real)XPad::AXIS_MAX * joy_sensitivity.x();
          } else {
            joy_pitch = NAN;
          }
          break;
        case XPad::AXIS_RIGHT_X:
          if(rabs(event.jaxis.value) > XPad::AXIS_RIGHT_DEAD_ZONE) {
            joy_yaw = event.jaxis.value * 2.0 * M_PI / (Real)XPad::AXIS_MAX * joy_sensitivity.y();
          } else {
            joy_yaw = NAN;
          }
          break;
        case XPad::AXIS_LEFT_TRIGGER:
          if(event.jaxis.value > XPad::AXIS_MIN) {
            joy_roll = -((event.jaxis.value - XPad::AXIS_MIN) * 2.0 * M_PI / (Real)XPad::AXIS_MAX * joy_sensitivity.z());
          } else {
            joy_roll = NAN;
          }
          break;
        case XPad::AXIS_RIGHT_TRIGGER:
          if(event.jaxis.value > XPad::AXIS_MIN) {
            joy_roll = (event.jaxis.value - XPad::AXIS_MIN) * 2.0 * M_PI / (Real)XPad::AXIS_MAX * joy_sensitivity.z();
          } else {
            joy_roll = NAN;
          }
          break;
        }
        break;
      case SDL_JOYBUTTONUP:
      case SDL_JOYBUTTONDOWN:
        switch(event.jbutton.button) {
        case XPad::BUTTON_A:
          player_command_state.movingUp(event.jbutton.state == SDL_PRESSED);
          break;
        case XPad::BUTTON_X:
          player_command_state.movingDown(event.jbutton.state == SDL_PRESSED);
          break;
        case XPad::BUTTON_B:
          if(event.jbutton.state == SDL_PRESSED) {
            speed *= 2.0;
            std::cout << "Speed: " << speed << std::endl;
          }
          break;
        case XPad::BUTTON_Y:
          if(event.jbutton.state == SDL_PRESSED) {
            speed *= 0.5;
            std::cout << "Speed: " << speed << std::endl;
          }
          break;
        case XPad::BUTTON_WHITE:
          if(event.type == SDL_JOYBUTTONDOWN) {
            if(camera == &fps_camera) {
              Camera *c = new Camera(*camera);
              camera = c;
            } else {
              fps_camera.setPosition(camera->position());
              fps_camera.lookAt(-Vec3::Z, Vec3::Y);
              delete camera;
              camera = &fps_camera;
            }
            view_camera = camera;
            std::cout << "Toggled flight" << std::endl;
          } break;
        case XPad::BUTTON_BACK:
          if(event.type == SDL_JOYBUTTONDOWN) {
            done = true;
          }
        }
      case SDL_JOYHATMOTION:
        //std::cout << "JoyHatMotion: " << event.jhat.which << "\t" << (int)event.jhat.hat << "\t" << (int)event.jhat.value << std::endl;
        break;
      default:
        //std::cout << "Unknown event: " << event.type << std::endl;
        break;
      }
    }

    Real sec_per_frame = ms_per_frame / 1000.0;

    if(!isnan(joy_roll)) {
      droll -= joy_roll * sec_per_frame;
    }
    if(!isnan(joy_pitch)) {
      dpitch -= joy_pitch * sec_per_frame;
    }
    if(!isnan(joy_yaw)) {
      dyaw -= joy_yaw * sec_per_frame;
    }

    Vec3 dir = player_command_state.motionVector();
    if(!dir.isNaN()) {
      camera->moveBy(dir * speed * (ms_per_frame / 1000.0));
    }
    camera->rotateBy(Vec3(dpitch, dyaw, droll));

    quanta.update(1.0);

    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(80.0, (Real)window_width / (Real)window_height, near_plane, far_plane);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glLoadMatrixr(camera->matrix().invert());
    //glMultMatrixf(camera->rotationMatrix().invert());

    glLightrv(GL_LIGHT0, GL_POSITION, Vec3(0,  far_plane * 2, 0, 1));

    Camera c(*view_camera);
    c.setPosition(Vec3());

    glPushMatrix(); {
      quanta_renderer.render(*view_camera, quanta);
      //quanta_renderer.renderVelocities(*view_camera, quanta);
    } glPopMatrix();

    SDL_GL_SwapWindow(window);
    ASSERT_GL_ERROR;

    screen_shooter.takeScreenShot();

    frames++;
    ms_per_frame = SDL_GetTicks() - frame_start;
  }

  SDL_GL_DeleteContext(gl_context);
  SDL_DestroyWindow(window);
  SDL_Quit();

  return 0;
}
