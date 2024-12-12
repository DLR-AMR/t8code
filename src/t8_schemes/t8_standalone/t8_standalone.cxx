#include <t8_schemes/t8_standalone/t8_standalone.hxx>
#include <t8_schemes/t8_scheme_builder.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

t8_scheme *
t8_scheme_new_standalone (void)
{
  t8_scheme_builder builder;

  builder.add_eclass_scheme<t8_standalone_scheme_c<T8_ECLASS_VERTEX>> ();
  builder.add_eclass_scheme<t8_standalone_scheme_c<T8_ECLASS_LINE>> ();
  builder.add_eclass_scheme<t8_standalone_scheme_c<T8_ECLASS_QUAD>> ();
  builder.add_eclass_scheme<t8_standalone_scheme_c<T8_ECLASS_HEX>> ();
  builder.add_eclass_scheme<t8_standalone_scheme_c<T8_ECLASS_TRIANGLE>> ();
  builder.add_eclass_scheme<t8_standalone_scheme_c<T8_ECLASS_TET>> ();
  builder.add_eclass_scheme<t8_standalone_scheme_c<T8_ECLASS_PRISM>> ();
  builder.add_eclass_scheme<t8_standalone_scheme_c<T8_ECLASS_PYRAMID>> ();
  return builder.build_scheme ();
}

int
t8_eclass_scheme_is_standalone (const t8_scheme *scheme, const t8_eclass_t eclass)
{
  switch (eclass) {
  case T8_ECLASS_VERTEX:
    return scheme->check_eclass_scheme_type<t8_standalone_scheme_c> (T8_ECLASS_VERTEX);
  case T8_ECLASS_LINE:
    return scheme->check_eclass_scheme_type<t8_standalone_scheme_c> (T8_ECLASS_LINE);
  case T8_ECLASS_QUAD:
    return scheme->check_eclass_scheme_type<t8_standalone_scheme_c> (T8_ECLASS_QUAD);
  case T8_ECLASS_TRIANGLE:
    return scheme->check_eclass_scheme_type<t8_standalone_scheme_c> (T8_ECLASS_TRIANGLE);
  case T8_ECLASS_HEX:
    return scheme->check_eclass_scheme_type<t8_standalone_scheme_c> (T8_ECLASS_HEX);
  case T8_ECLASS_TET:
    return scheme->check_eclass_scheme_type<t8_standalone_scheme_c> (T8_ECLASS_TET);
  case T8_ECLASS_PRISM:
    return scheme->check_eclass_scheme_type<t8_standalone_scheme_c> (T8_ECLASS_PRISM);
  case T8_ECLASS_PYRAMID:
    return scheme->check_eclass_scheme_type<t8_standalone_scheme_c> (T8_ECLASS_PYRAMID);
  default:
    SC_ABORT_NOT_REACHED ();
  }
  return 0; /* Default return value false */
}

T8_EXTERN_C_END ();
