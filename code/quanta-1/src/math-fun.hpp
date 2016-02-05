#ifndef _MATH_FUN_HPP_
#define _MATH_FUN_HPP_

template<typename T>
T clamp(const T &value, const T &max, const T &min = 0)
{
    if(value < min) return min;
    else if(value > max) return max;
    else return value;
}

template<typename A>
A interpolate(const A &a, const A &b, float t)
{
  return a * (1.0f - t) + b * t;
}

template<typename A>
A cubic_interpolate(const A &a, const A &b, float t)
{
  float t2 = t * t, t3 = t2 * t;
  return (2*t3 - 3*t2 + 1) * a + (-2 * t3 + 3 * t2) * b;
}

float signed_ceil(float v)
{
  if(v < 0.0f) {
    return floorf(v);
  } else {
    return ceilf(v);
  }
}

#endif /* _MATH_FUN_HPP */
