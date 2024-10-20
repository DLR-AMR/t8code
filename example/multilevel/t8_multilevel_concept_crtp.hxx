#pragma once

template <class Derived_t>
class crtp {
  Derived_t&
  underlying ()
  {
    return static_cast<Derived_t&> (*this);
  }
  Derived_t const&
  underlying () const
  {
    return static_cast<Derived_t const&> (*this);
  }
};
