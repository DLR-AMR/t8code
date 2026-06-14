#ifdef T8_ENABLE_MRA

#include <cmath>
#include "t8_mra/num/basis_functions.hxx"

/* Note from Veli Ünlü about rights: These are the basis functions from
 * the bachelor thesis in mathematics of Florian Sieglar, "Konstruktion von Multiwavelets auf Dreiecksgittern"
 * from November 2013.
 */

namespace t8_mra
{
double
scaling_function (int i, double tau1, double tau2)
{
  switch (i) {
  case 0:
    return std::sqrt (2.);
  case 1:
    return 6. * tau2 - 2.;
  case 2:
    return 4. * (0.5 - tau1 - 0.5 * tau2) * std::sqrt (3.);
  case 3:
    return (10. * (tau2 * tau2 + 1. / 10. - (4. / 5.) * tau2)) * std::sqrt (6.);
  case 4:
    return (30.
            * (tau2 * (1. - tau1 - tau2) - 1. / 10. - (2. / 5.) * tau2 + (1. / 5.) * tau1 + (1. / 2.) * tau2 * tau2))
           * std::sqrt (2.);
  case 5:
    return (6.
            * ((1. - tau1 - tau2) * (1. - tau1 - tau2) - 5. / 6. + (2. / 3.) * tau2 + tau1 + (1. / 6.) * tau2 * tau2
               + tau2 * (1. - tau1 - tau2)))
           * std::sqrt (30.);
  case 6:
    return (70. * (tau2 * tau2 * tau2 - 1. / 35. + (3. / 7.) * tau2 - (9. / 7.) * tau2 * tau2)) * std::sqrt (2.);
  case 7:
    return (84.
            * (tau2 * tau2 * (1. - tau1 - tau2) + 1. / 42. + (11. / 42.) * tau2 - (1. / 21.) * tau1
               - (11. / 14.) * tau2 * tau2 - (4. / 7.) * tau2 * (1. - tau1 - tau2) + (1. / 2.) * tau2 * tau2 * tau2))
           * std::sqrt (6.);
  case 8:
    return (84.
            * (tau2 * (1. - tau1 - tau2) * (1. - tau1 - tau2) + 5. / 42. + (1. / 14.) * tau2 - (1. / 7.) * tau1
               - (5. / 14.) * tau2 * tau2 - (1. / 7.) * (1. - tau1 - tau2) * (1. - tau1 - tau2)
               + (1. / 6.) * tau2 * tau2 * tau2 - (8. / 7.) * tau2 * (1. - tau1 - tau2)
               + tau2 * tau2 * (1. - tau1 - tau2)))
           * std::sqrt (10.);
  case 9:
    return (40.
            * ((1. - tau1 - tau2) * (1. - tau1 - tau2) * (1. - tau1 - tau2) + 11. / 20. - (9. / 20.) * tau2
               - (3. / 5.) * tau1 - (3. / 20.) * tau2 * tau2 - (3. / 2.) * (1. - tau1 - tau2) * (1. - tau1 - tau2)
               + (1. / 20.) * tau2 * tau2 * tau2 - (6. / 5.) * tau2 * (1. - tau1 - tau2)
               + (3. / 5.) * tau2 * tau2 * (1. - tau1 - tau2)
               + (3. / 2.) * tau2 * (1. - tau1 - tau2) * (1. - tau1 - tau2)))
           * std::sqrt (14.);
  default:;
  };
}


}  // namespace t8_mra

#endif
