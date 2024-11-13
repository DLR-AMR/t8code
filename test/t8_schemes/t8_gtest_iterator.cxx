#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <vector>

/** 
 * \class scheme_iterators
 * Class to iterate over all eclasses of all schemes an return the t8_eclass_scheme_c*.
 */
class scheme_iterators {
 public:
  /**
     * Initialize the iterator with a list of schemes.
     * \param [in] schemes The list of schemes to iterate over.
    */
  scheme_iterators (const std::vector<const t8_scheme_cxx*>& schemes): schemes (schemes)
  {
  }

  /**
     * \struct Iterator
     * Iterator to iterate over all eclasses of all schemes.
    */
  struct Iterator
  {
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = t8_eclass_scheme_c*;
    using pointer = t8_eclass_scheme_c**;
    using reference = t8_eclass_scheme_c*&;

    /**
         * Constructor for the iterator.
         * \param [in] schemes The list of schemes to iterate over.
         * \param [in] is_end Flag to indicate if the iterator is at the end.
        */
    Iterator (const std::vector<const t8_scheme_cxx*>& schemes, bool is_end = false)
      : schemes (schemes), scheme_index (is_end ? schemes.size () : 0), eclass_index (0)
    {
      if (!is_end && !schemes.empty ()) {
        eclass_count = schemes[scheme_index]->num_eclasses;
      }
    }

    /**
         * Dereference operator.
         * \return The current eclass scheme or tree scheme ts.
        */
    t8_eclass_scheme_c*
    operator* () const
    {
      const t8_scheme_cxx* current_scheme = schemes[scheme_index];
      return current_scheme->eclass_schemes[eclass_index];
    }

    /**
        *  Prefix increment operator to move the iterator to the next element. 
        * \return A reference to the updated iterator. 
        */
    Iterator&
    operator++ ()
    {
      if (++eclass_index >= eclass_count) {
        eclass_index = 0;
        if (++scheme_index < schemes.size ()) {
          eclass_count = schemes[scheme_index]->num_eclasses;
        }
      }
      return *this;
    }

    /**
        * Inequality operator to compare two iterators.
        * \param [in] other Another iterator to compare with.
        * \return True if the iterators are not equal, false otherwise. 
        */
    bool
    operator!= (const Iterator& other) const
    {
      return scheme_index != other.scheme_index || eclass_index != other.eclass_index;
    }

   private:
    const std::vector<const t8_scheme_cxx*>& schemes;
    size_t scheme_index;
    size_t eclass_index;
    size_t eclass_count;
  };

  Iterator
  begin () const
  {
    return Iterator (schemes);
  }
  Iterator
  end () const
  {
    return Iterator (schemes, true);
  }

 private:
  const std::vector<const t8_scheme_cxx*>& schemes;
};
