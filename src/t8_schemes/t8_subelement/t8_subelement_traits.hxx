
#pragma once

#include <t8_schemes/t8_standalone/t8_standalone_elements.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_implementation.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_schemes/t8_subelement/t8_subelement_type.hxx>

struct t8_subelementquad_scheme;

struct t8_subelementtri_scheme;

template <typename TScheme>
struct t8_subelement_traits;

template <>
struct t8_subelement_traits<t8_subelementquad_scheme>
{
  using UnderlyingScheme = t8_standalone_scheme<T8_ECLASS_QUAD>;
  using SubelementType = t8_subelement_element<t8_standalone_element<T8_ECLASS_QUAD>>;
};

template <>
struct t8_subelement_traits<t8_subelementtri_scheme>
{
  using UnderlyingScheme = t8_default_scheme_tri;
  using SubelementType = t8_subelement_element<t8_dtri>;
};
