
#pragma once

#include <t8_schemes/t8_standalone/t8_standalone_elements.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_implementation.hxx>
#include <t8_eclass/t8_eclass.h>
#include <t8_schemes/t8_subelement/t8_subelement_type.hxx>

// Forward declaration reicht hier
struct t8_subelementquad_scheme;

// Primäres Template (undefiniert — erzeugt klaren Fehler bei fehlendem Trait)
template <typename TScheme>
struct t8_subelement_traits;

// Spezialisierung für quad — BEVOR t8_subelementquad_scheme definiert wird
template <>
struct t8_subelement_traits<t8_subelementquad_scheme>
{
  using UnderlyingScheme = t8_standalone_scheme<T8_ECLASS_QUAD>;
  using SubelementType = t8_subelement_element<t8_standalone_element<T8_ECLASS_QUAD>>;
};
