#pragma once

// Matrixklasse mit LR-Zerlegung, entnommen aus
// http://www.igpm.rwth-aachen.de/Download/ss13/na2/na2-base.h Numerische
// Analysis II, SS 2013, Prof. Dr. Wolfgang Dahmen, Dr. Markus Bachmayr

#ifdef T8_ENABLE_MRA

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace t8_mra
{

class vec {
  std::vector<double> data;

 public:
  vec () = default;
  vec (size_t n): data (n, 0.0)
  {
  }
  vec (std::initializer_list<double> l): data (l)
  {
  }

  vec (const vec&) = default;
  vec&
  operator= (const vec&)
    = default;
  vec (vec&&) = default;
  vec&
  operator= (vec&&)
    = default;

  vec&
  operator= (const std::initializer_list<double>& l)
  {
    if (l.size () != data.size ())
      throw std::out_of_range ("number elements in t8_mra::util::vec does not fit to number columns "
                               "and number rows");

    std::copy_n (l.begin (), l.size (), data.begin ());

    return *this;
  }

  double&
  operator() (size_t i)
  {
    if (i >= data.size ())
      throw std::out_of_range ("index in t8_mra::util::vec::operator() is out of range");

    return data[i];
  }

  double
  operator() (size_t i) const
  {
    if (i >= data.size ())
      throw std::out_of_range ("index in t8_mra::util::vec::operator() is out of range");

    return data[i];
  }

  size_t
  size () const noexcept
  {
    return data.size ();
  }

  vec&
  operator+= (const vec& y)
  {
    if (y.size () < data.size ())
      throw std::logic_error ("lengths in t8_mra::util::vec::operator+= does not fit");

    for (auto i = 0u; i < data.size (); i++)
      data[i] += y (i);

    return *this;
  }

  vec&
  operator-= (const vec& y)
  {
    if (y.size () < data.size ())
      throw std::logic_error ("lengths in t8_mra::util::vec::operator+= does not fit");
    for (auto i = 0u; i < data.size (); i++)
      data[i] -= y (i);

    return *this;
  }

  vec&
  operator*= (double v)
  {
    for (auto i = 0u; i < data.size (); i++)
      data[i] *= v;

    return *this;
  }

  void
  resize (size_t n)
  {
    data.clear ();
    data.resize (n);
  }
};

inline double
inner (const vec& v1, const vec& v2)
{
  if (v1.size () != v2.size ())
    throw std::logic_error ("lengths in t8_mra::util::inner does not fit");

  auto res = 0.0;
  for (auto i = 0u; i < v1.size (); i++)
    res += v1 (i) * v2 (i);

  return res;
}

inline double
l1norm (const vec& v)
{
  auto res = 0.0;
  for (auto i = 0u; i < v.size (); i++)
    res += std::abs (v (i));

  return res;
}

inline double
l2norm (const vec& v)
{
  return std::sqrt (inner (v, v));
}

inline double
linftynorm (const vec& v)
{
  auto res = 0.0;
  for (auto i = 0u; i < v.size (); i++) {
    const auto tmp = std::abs (v (i));

    if (tmp > res)
      res = tmp;
  }

  return res;
}

}  // namespace t8_mra

#endif
