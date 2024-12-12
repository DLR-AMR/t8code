#ifndef T8_STANDALONE_ELEMENTS_HXX
#define T8_STANDALONE_ELEMENTS_HXX

constexpr uint8_t T8_ELEMENT_DIM[T8_ECLASS_COUNT] = { 0, 1, 2, 2, 3, 3, 3, 3 };
constexpr uint8_t T8_ELEMENT_MAXLEVEL[T8_ECLASS_COUNT] = { 255, 30, 29, 29, 21, 21, 21, 18 };

constexpr uint8_t T8_ELEMENT_NUM_CHILDREN[T8_ECLASS_COUNT] = { 1, 2, 4, 4, 8, 8, 8, 10 };

constexpr uint8_t T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_COUNT] = { 0, 0, 0, 1, 0, 3, 1, 2 };

typedef uint32_t t8_element_coord_t;
typedef uint8_t t8_element_level_t;
typedef int8_t t8_cube_id_t;

template <t8_eclass_t eclass_T>
using t8_element_type_t = std::bitset<T8_ELEMENT_NUM_EQUATIONS[eclass_T]>;

template <t8_eclass_t eclass_T>
using t8_element_coords_t = std::array<t8_element_coord_t, T8_ELEMENT_DIM[eclass_T]>;

#endif /* T8_STANDALONE_ELEMENTS_HXX */
