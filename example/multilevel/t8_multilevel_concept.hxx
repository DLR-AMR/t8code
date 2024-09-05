struct triangle
{
  int level;
  int orientation;
};

struct quad
{
  int level;
};

// Multilevel element, which adds the actual hierarchical level of the element
// The level inside the element will determine its size
template <class eclass>
struct multilevel_element
{
  eclass elem;
  int multilevel_level;
};

// Opaque pointer to cast the element type into
typedef struct element element_t;

// Using CRTP to avoid virtual function calls
template <class derived_scheme_t>
class scheme_base {
 public:
  ~scheme_base ()
  {
  }

  int
  get_level (element_t *elem)
  {
    // cast derived class into base class to avoid virtual functions
    return static_cast<derived_scheme_t const &> (*this).get_level (elem);
  };

  int
  get_num_children (element_t *elem)
  {
    // cast derived class into base class to avoid virtual functions
    return static_cast<derived_scheme_t const &> (*this).get_num_children (elem);
  };

  int
  get_num_vertices ()
  {
    // cast derived class into base class to avoid virtual functions
    return static_cast<derived_scheme_t const &> (*this).get_num_vertices ();
  };

 private:
  // This way the derived class and only the derived class can use this constructor.
  // This way you get an error when doing: `class triangle_scheme: public scheme_base <quad_scheme>`
  scheme_base () {};
  friend derived_scheme_t;
};

// inherits from base which is a template for this class
class triangle_scheme: public scheme_base<triangle_scheme> {
 public:
  triangle_scheme () {};

  int
  get_level (element_t *elem)
  {
    elem_type *tri = (elem_type *) elem;
    return tri->level;
  };

  int
  get_num_children (element_t *elem)
  {
    return 4;
  };

  int
  get_num_vertices ()
  {
    return 3;
  };

 protected:
  // When the multilevel class inherits from this it needs to know the element type of this scheme
  using elem_type = triangle;

 private:
};

// inherits from base which is a template for this class
class quad_scheme: public scheme_base<quad_scheme> {
 public:
  quad_scheme () {};

  int
  get_level (element_t *elem)
  {
    elem_type *q = (elem_type *) elem;
    return q->level;
  };

  int
  get_num_children (element_t *elem)
  {
    return 4;
  };

  int
  get_num_vertices ()
  {
    return 4;
  };

 protected:
  // When the multilevel class inherits from this it needs to know the element type of this scheme
  using elem_type = quad;

 private:
};

template <class scheme_impl_t>
class multilevel_scheme: public scheme_base<multilevel_scheme<scheme_impl_t>>, private scheme_impl_t {
 public:
  using scheme_impl_t::scheme_impl_t;
  using elem_type = typename scheme_impl_t::elem_type;
  ~multilevel_scheme () {};

  int
  get_level (element_t *elem)
  {
    multilevel_element<elem_type> *m_elem = (multilevel_element<elem_type> *) elem;
    return m_elem->multilevel_level;
  };

  int
  get_num_children (element_t *elem)
  {
    multilevel_element<elem_type> *m_elem = (multilevel_element<elem_type> *) elem;
    const int elem_level = scheme_impl_t::get_level ((element_t *) &(m_elem->elem));
    if (elem_level == get_level (elem)) {
      return scheme_impl_t::get_num_children ((element_t *) &(m_elem->elem)) + 1;
    }
    else {
      return 0;
    }
  };

  int
  get_num_vertices ()
  {
    return scheme_impl_t::get_num_vertices ();
  };

 private:
};
