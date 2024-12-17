#ifndef T8_STANDALONE_HXX
#define T8_STANDALONE_HXX

#include <t8_schemes/t8_scheme.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_implementation.hxx>

T8_EXTERN_C_BEGIN ();

/** Return the standalone element implementation of t8code. */
t8_scheme *
t8_scheme_new_standalone (void);

/** Check whether a given eclass_scheme is one of the standalone schemes.
 * \param [in] scheme   A (pointer to a) scheme
 * \param [in] eclass   The eclass to check
 * \return              True (non-zero) if \a scheme is one of the default schemes,
 *                      false (zero) otherwise.
 */
int
t8_eclass_scheme_is_standalone (const t8_scheme *scheme, const t8_eclass_t eclass);

T8_EXTERN_C_END ();

#endif /* !T8_STANDALONE_HXX */
