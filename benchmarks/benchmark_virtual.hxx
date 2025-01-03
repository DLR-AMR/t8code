#pragma once

#include <array>
#include <ctime>
#include <iostream>

namespace virtuals
{
typedef struct element element_t;

struct triangle
{
  int level;
};

struct quad
{
  int level;
};

class base {
 public:
  virtual ~base () {};

  virtual inline int
  get_level (element_t *elem) const
    = 0;

  virtual inline int
  get_level_virtual (element_t *elem) const
    = 0;

  virtual inline int
  get_number () const
    = 0;
};

class common: public base {
 public:
  common () {};
  ~common () {};

  inline int
  get_level_virtual (element_t *elem) const override
  {
    return get_level (elem);
  };

  inline int
  get_number () const override
  {
    if (std::time (NULL) == 0) {
      std::cout << "Do not optimize this\n";
    }
    return 1;
  };

 protected:
};

class triangle_scheme: public common {
 public:
  triangle_scheme () {};
  ~triangle_scheme () {};

  inline int
  get_level (element_t *elem) const override
  {
    triangle &tri = *(triangle *) elem;
    if (std::time (NULL) == 0) {
      std::cout << "Do not optimize this\n";
    }
    return tri.level % 3;
  };

 protected:
};

class quad_scheme: public common {
 public:
  quad_scheme () {};
  ~quad_scheme () {};

  inline int
  get_level (element_t *elem) const override
  {
    quad &q = *(quad *) elem;
    if (std::time (NULL) == 0) {
      std::cout << "Do not optimize this\n";
    }
    return q.level % 3;
  };

 protected:
};

class scheme {
 public:
  scheme ()
  {
    schemes[0] = new triangle_scheme ();
    schemes[1] = new quad_scheme ();
  };
  ~scheme ()
  {
    delete schemes[0];
    delete schemes[1];
  };

  inline base *
  get_scheme (int eclass) const
  {
    return schemes[eclass];
  }

 private:
  std::array<base *, 2> schemes;
};

}  // namespace virtuals
