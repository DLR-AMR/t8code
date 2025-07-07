#include <cmath>
#include "basis_functions.hxx"
//T8_EXTERN_C_BEGIN ();

/* Note from Veli Ünlü about rights: These are the basis functions from
 * the bachelor thesis in mathematics of Florian Sieglar, "Konstruktion von Multiwavelets auf Dreiecksgittern"
 * from November 2013.
 */

double
skalierungsfunktion (int i, double tau1, double tau2)
{
  switch (i) {
  case 0:
    return sqrt (2.);
  case 1:
    return 6. * tau2 - 2.;
  case 2:
    return 4. * (0.5 - tau1 - 0.5 * tau2) * sqrt (3.);
  case 3:
    return (10. * (tau2 * tau2 + 1. / 10. - (4. / 5.) * tau2)) * sqrt (6.);
  case 4:
    return (30.
            * (tau2 * (1. - tau1 - tau2) - 1. / 10. - (2. / 5.) * tau2 + (1. / 5.) * tau1 + (1. / 2.) * tau2 * tau2))
           * sqrt (2.);
  case 5:
    return (6.
            * ((1. - tau1 - tau2) * (1. - tau1 - tau2) - 5. / 6. + (2. / 3.) * tau2 + tau1 + (1. / 6.) * tau2 * tau2
               + tau2 * (1. - tau1 - tau2)))
           * sqrt (30.);
  case 6:
    return (70. * (tau2 * tau2 * tau2 - 1. / 35. + (3. / 7.) * tau2 - (9. / 7.) * tau2 * tau2)) * sqrt (2.);
  case 7:
    return (84.
            * (tau2 * tau2 * (1. - tau1 - tau2) + 1. / 42. + (11. / 42.) * tau2 - (1. / 21.) * tau1
               - (11. / 14.) * tau2 * tau2 - (4. / 7.) * tau2 * (1. - tau1 - tau2) + (1. / 2.) * tau2 * tau2 * tau2))
           * sqrt (6.);
  case 8:
    return (84.
            * (tau2 * (1. - tau1 - tau2) * (1. - tau1 - tau2) + 5. / 42. + (1. / 14.) * tau2 - (1. / 7.) * tau1
               - (5. / 14.) * tau2 * tau2 - (1. / 7.) * (1. - tau1 - tau2) * (1. - tau1 - tau2)
               + (1. / 6.) * tau2 * tau2 * tau2 - (8. / 7.) * tau2 * (1. - tau1 - tau2)
               + tau2 * tau2 * (1. - tau1 - tau2)))
           * sqrt (10.);
  case 9:
    return (40.
            * ((1. - tau1 - tau2) * (1. - tau1 - tau2) * (1. - tau1 - tau2) + 11. / 20. - (9. / 20.) * tau2
               - (3. / 5.) * tau1 - (3. / 20.) * tau2 * tau2 - (3. / 2.) * (1. - tau1 - tau2) * (1. - tau1 - tau2)
               + (1. / 20.) * tau2 * tau2 * tau2 - (6. / 5.) * tau2 * (1. - tau1 - tau2)
               + (3. / 5.) * tau2 * tau2 * (1. - tau1 - tau2)
               + (3. / 2.) * tau2 * (1. - tau1 - tau2) * (1. - tau1 - tau2)))
           * sqrt (14.);
  default:;
  };
}

double
muttermultiwavelets (int p, int i, double tau1, double tau2, int e)
{
  switch (p) {
  case 1:
    switch (i) {
    case 0:
      if (e == 0)
        return 0.;
      else if (e == 1)
        return (2. / 3.) * sqrt (3.);
      else if (e == 2)
        return -(4. / 3.) * sqrt (3.);
      else
        return (2. / 3.) * sqrt (3.);
    case 1:
      if (e == 0)
        return 0.;
      else if (e == 1)
        return 2.;
      else if (e == 2)
        return 0.;
      else
        return -2.;
    case 2:
      if (e == 0)
        return -sqrt (2.) * sqrt (3.);
      else if (e == 1)
        return (1. / 3.) * sqrt (2.) * sqrt (3.);
      else if (e == 2)
        return (1. / 3.) * sqrt (2.) * sqrt (3.);
      else
        return (1. / 3.) * sqrt (2.) * sqrt (3.);
    default:;
    };
  case 2:
    switch (i) {
    case 0:
      if (e == 0)
        return -(1. / 5.) * sqrt (2.) * (-1. + 8. * tau2) * sqrt (5.);
      else if (e == 1)
        return -(1. / 5.) * sqrt (2.) * (-3. + 16. * tau2) * sqrt (5.);
      else if (e == 2)
        return (3. / 5.) * sqrt (2.) * (-5. + 8. * tau2) * sqrt (5.);
      else
        return -(1. / 5.) * sqrt (2.) * (-3. + 16. * tau2) * sqrt (5.);
    case 1:
      if (e == 0)
        return (4. / 5.) * sqrt (2.) * (-1. + tau2 + 2. * tau1) * sqrt (15.);
      else if (e == 1)
        return -(1. / 5.) * sqrt (2.) * (-1. - 8. * tau2 + 4. * tau1) * sqrt (15.);
      else if (e == 2)
        return (8. / 5.) * sqrt (2.) * (-1. + tau2 + 2. * tau1) * sqrt (15.);
      else
        return -(1. / 5.) * sqrt (2.) * (-3. + 12. * tau2 + 4. * tau1) * sqrt (15.);
    case 2:
      if (e == 0)
        return -2. * sqrt (2.) + 4. * sqrt (2.) * tau2;
      else if (e == 1)
        return -9. * sqrt (2.) + 8. * sqrt (2.) * tau2 + 12. * sqrt (2.) * tau1;
      else if (e == 2)
        return 0.;
      else
        return 3. * sqrt (2.) - 4. * sqrt (2.) * tau2 - 12. * sqrt (2.) * tau1;
    case 3:
      if (e == 0)
        return (2. / 35. * (16. * tau2 - 7.)) * sqrt (35.);
      else if (e == 1)
        return -(8. / 35. * (7. * tau2 - 1.)) * sqrt (35.);
      else if (e == 2)
        return -(2. / 35. * (48. * tau2 - 35.)) * sqrt (35.);
      else
        return -(8. / 35. * (7. * tau2 - 1.)) * sqrt (35.);
    case 4:
      if (e == 0)
        return -(32. / 965. * (-1. + tau2 + 2. * tau1)) * sqrt (2895.);
      else if (e == 1)
        return (2. / 965. * (16. * tau1 - 92. * tau2 + 1.)) * sqrt (2895.);
      else if (e == 2)
        return (176. / 965. * (-1. + tau2 + 2. * tau1)) * sqrt (2895.);
      else
        return (2. / 965. * (16. * tau1 + 108. * tau2 - 17.)) * sqrt (2895.);
    case 5:
      if (e == 0)
        return -(2. / 203. * (7. + 8. * tau2)) * sqrt (406.);
      else if (e == 1)
        return -(1. / 203. * (84. * tau1 - 59. - 28. * tau2)) * sqrt (406.);
      else if (e == 2)
        return -(4. / 203. * (-7. + 9. * tau2)) * sqrt (406.);
      else
        return (1. / 203. * (84. * tau1 - 25. + 112. * tau2)) * sqrt (406.);
    case 6:
      if (e == 0)
        return (144. / 5597. * (-1. + tau2 + 2. * tau1)) * sqrt (11194.);
      else if (e == 1)
        return -(1. / 5597. * (-763. + 916. * tau1 + 716. * tau2)) * sqrt (11194.);
      else if (e == 2)
        return -(20. / 5597. * (-1. + tau2 + 2. * tau1)) * sqrt (11194.);
      else
        return -(1. / 5597. * (-153. + 916. * tau1 + 200. * tau2)) * sqrt (11194.);
    case 7:
      if (e == 0)
        return (22. / 29.) * sqrt (3.) * (3. * tau2 - 1.) * sqrt (29.);
      else if (e == 1)
        return -(2. / 29.) * sqrt (3.) * (7. * tau2 + 8. * tau1 - 7.) * sqrt (29.);
      else if (e == 2)
        return (2. / 29.) * sqrt (3.) * (-7. + 9. * tau2) * sqrt (29.);
      else
        return (2. / 29.) * sqrt (3.) * (tau2 + 8. * tau1 - 1.) * sqrt (29.);
    case 8:
      if (e == 0)
        return -(66. / 29. * (-1. + tau2 + 2. * tau1)) * sqrt (29.);
      else if (e == 1)
        return -(2. / 29. * (-7. + 10. * tau1 - 7. * tau2)) * sqrt (29.);
      else if (e == 2)
        return (14. / 29. * (-1. + tau2 + 2. * tau1)) * sqrt (29.);
      else
        return -(2. / 29. * (-3. + 10. * tau1 + 17. * tau2)) * sqrt (29.);
    default:;
    };
  case 3:
    switch (i) {
    case 0:
      if (e == 0)
        return -(8. / 21.) * sqrt (2.) * (1. - 15. * tau2 + 30. * tau2 * tau2) * sqrt (7.);
      else if (e == 1)
        return -(2. / 21.) * sqrt (2.) * (7. + 180. * tau2 * tau2 - 90. * tau2) * sqrt (7.);
      else if (e == 2)
        return (8. / 21.) * sqrt (2.) * (22. - 75. * tau2 + 60. * tau2 * tau2) * sqrt (7.);
      else
        return -(2. / 21.) * sqrt (2.) * (7. + 180. * tau2 * tau2 - 90. * tau2) * sqrt (7.);
    case 1:
      if (e == 0)
        return -(8. / 7. * (-3. * tau2 + 1. + 4. * tau2 * tau1 - 2. * tau1 + 2. * tau2 * tau2)) * sqrt (2.)
               * sqrt (21.);
      else if (e == 1)
        return (2. / 21.) * sqrt (2.) * (5. - 36. * tau2 * tau2 - 30. * tau2 - 12. * tau1 + 96. * tau2 * tau1)
               * sqrt (21.);
      else if (e == 2)
        return -(40. / 7. * (-3. * tau2 + 1. + 4. * tau2 * tau1 - 2. * tau1 + 2. * tau2 * tau2)) * sqrt (2.)
               * sqrt (21.);
      else
        return (2. / 21.) * sqrt (2.) * (7. + 132. * tau2 * tau2 - 78. * tau2 - 12. * tau1 + 96. * tau2 * tau1)
               * sqrt (21.);
    case 2:
      if (e == 0)
        return (8. / 21.) * sqrt (2.)
               * (7. - 21. * tau2 - 24. * tau1 + 18. * tau2 * tau2 + 24. * tau2 * tau1 + 24. * tau1 * tau1)
               * sqrt (35.);
      else if (e == 1)
        return -(2. / 21.) * sqrt (2.)
               * (-5. + 126. * tau2 - 12. * tau1 - 108. * tau2 * tau2 - 144. * tau2 * tau1 + 24. * tau1 * tau1)
               * sqrt (35.);
      else if (e == 2)
        return (16. / 7.) * sqrt (2.) * (6. * tau2 * tau1 - 2. * tau2 + tau2 * tau2 + 6. * tau1 * tau1 - 6. * tau1 + 1.)
               * sqrt (35.);
      else
        return -(2. / 21.) * sqrt (2.)
               * (7. - 54. * tau2 - 36. * tau1 + 60. * tau2 * tau2 + 192. * tau2 * tau1 + 24. * tau1 * tau1)
               * sqrt (35.);
    case 3:
      if (e == 0)
        return (8. * (-3. * tau2 + 1. + 4. * tau2 * tau1 - 2. * tau1 + 2. * tau2 * tau2)) * sqrt (2.);
      else if (e == 1)
        return (2. / 3.
                * (-90. * tau2 + 55. + 144. * tau2 * tau1 - 168. * tau1 + 120. * tau1 * tau1 + 36. * tau2 * tau2))
               * sqrt (2.);
      else if (e == 2)
        return 0.0;
      else
        return -(2. / 3. * (-18. * tau2 + 7. + 96. * tau2 * tau1 - 72. * tau1 + 120. * tau1 * tau1 + 12. * tau2 * tau2))
               * sqrt (2.);
    case 4:
      if (e == 0)
        return (5. / 148533. * (-169. + 456. * tau2 + 600. * tau2 * tau2)) * sqrt (99022.);
      else if (e == 1)
        return (1. / 148533. * (931. - 15480. * tau2 + 38520. * tau2 * tau2)) * sqrt (99022.);
      else if (e == 2)
        return (1. / 148533. * (29363. - 87000. * tau2 + 62040. * tau2 * tau2)) * sqrt (99022.);
      else
        return (1. / 148533. * (931. - 15480. * tau2 + 38520. * tau2 * tau2)) * sqrt (99022.);
    case 5:
      if (e == 0)
        return -(20. / 119. * (-4. * tau1 + 2. + 14. * tau2 * tau1 + 7. * tau2 * tau2 - 9. * tau2)) * sqrt (714.);
      else if (e == 1)
        return (1. / 357. * (-60. * tau1 + 17. + 600. * tau2 * tau1 - 900. * tau2 * tau2)) * sqrt (714.);
      else if (e == 2)
        return (100. / 119. * (-6. * tau1 + 3. + 10. * tau2 * tau1 + 5. * tau2 * tau2 - 8. * tau2)) * sqrt (714.);
      else
        return (1. / 357. * (-60. * tau1 + 43. + 600. * tau2 * tau1 + 1500. * tau2 * tau2 - 660. * tau2)) * sqrt (714.);
    case 6:
      if (e == 0)
        return (2. / 13630367358069.
                * (11771011. - 45917916. * tau1 + 45917916. * tau1 * tau1 - 12669510. * tau2 - 19247958. * tau2 * tau2
                   + 45917916. * tau2 * tau1))
               * sqrt (45434557860230.);
      else if (e == 1)
        return (1. / 13630367358069.
                * (-10972727. + 2121900. * tau1 + 17144952. * tau1 * tau1 + 195960240. * tau2 - 93890124. * tau2 * tau2
                   - 263285352. * tau2 * tau1))
               * sqrt (45434557860230.);
      else if (e == 2)
        return (4. / 4543455786023.
                * (4510579. - 31531434. * tau1 + 31531434. * tau1 * tau1 - 8359658. * tau2 + 3758359. * tau2 * tau2
                   + 31531434. * tau2 * tau1))
               * sqrt (45434557860230.);
      else
        return (1. / 13630367358069.
                * (8294125. - 36411804. * tau1 + 17144952. * tau1 * tau1 - 103736916. * tau2 + 186540180. * tau2 * tau2
                   + 297575256. * tau2 * tau1))
               * sqrt (45434557860230.);
    case 7:
      if (e == 0)
        return (8. / 710311.) * sqrt (2.)
               * (979. * tau2 * tau2 - 1572. * tau2 + 593. - 1186. * tau1 + 1958. * tau2 * tau1) * sqrt (3551555.);
      else if (e == 1)
        return -(1. / 2130933.) * sqrt (2.)
               * (-104904. * tau2 * tau2 + 9911. + 104652. * tau2 + 39984. * tau1 * tau1 - 45876. * tau1
                  - 104688. * tau2 * tau1)
               * sqrt (3551555.);
      else if (e == 2)
        return -(12. / 710311.) * sqrt (2.) * (6. * tau2 * tau2 + 7. - 13. * tau2 - 14. * tau1 + 12. * tau2 * tau1)
               * sqrt (3551555.);
      else
        return (1. / 2130933.) * sqrt (2.)
               * (39768. * tau2 * tau2 + 4019. - 34128. * tau2 + 39984. * tau1 * tau1 - 34092. * tau1
                  + 184656. * tau2 * tau1)
               * sqrt (3551555.);
    case 8:
      if (e == 0)
        return (6. / 8002149779950753.) * sqrt (2.)
               * (5585291. - 74690060. * tau2 * tau2 + 34995320. * tau2 + 99240520. * tau2 * tau1 - 99240520. * tau1
                  + 99240520. * tau1 * tau1)
               * sqrt (40010748899753765.);
      else if (e == 1)
        return -(3. / 8002149779950753.) * sqrt (2.)
               * (891343943. + 629296536. * tau2 * tau2 - 1507965276. * tau2 + 2094287056. * tau2 * tau1
                  - 2411211524. * tau1 + 1572337984. * tau1 * tau1)
               * sqrt (40010748899753765.);
      else if (e == 2)
        return -(12. / 8002149779950753.) * sqrt (2.)
               * (730663. + 248486. * tau2 * tau2 - 1037301. * tau2 + 7247964. * tau2 * tau1 - 7247964. * tau1
                  + 7247964. * tau1 * tau1)
               * sqrt (40010748899753765.);
      else
        return -(3. / 8002149779950753.) * sqrt (2.)
               * (52470403. + 107347464. * tau2 * tau2 - 147142664. * tau2 + 1050388912. * tau2 * tau1
                  - 733464444. * tau1 + 1572337984. * tau1 * tau1)
               * sqrt (40010748899753765.);
    case 9:
      if (e == 0)
        return (6. / 3267475795710091.
                * (-227001675. * tau2 - 31397600. * tau1 + 338286776. * tau2 * tau2 + 34827437.
                   + 31397600. * tau2 * tau1 + 31397600. * tau1 * tau1))
               * sqrt (16337378978550455.);
      else if (e == 1)
        return (6. / 3267475795710091.
                * (165910845. * tau2 - 7535320. * tau1 - 320509544. * tau2 * tau2 - 6185613. - 84778720. * tau2 * tau1
                   + 11385920. * tau1 * tau1))
               * sqrt (16337378978550455.);
      else if (e == 2)
        return (6. / 3267475795710091.
                * (267081993. - 715539635. * tau2 - 132269280. * tau1 + 470342264. * tau2 * tau2
                   + 132269280. * tau2 * tau1 + 132269280. * tau1 * tau1))
               * sqrt (16337378978550455.);
      else
        return (6. / 3267475795710091.
                * (65895605. * tau2 - 224344904. * tau2 * tau2 - 2335013. + 107550560. * tau2 * tau1
                   + 11385920. * tau1 * tau1 - 15236520. * tau1))
               * sqrt (16337378978550455.);
    case 10:
      if (e == 0)
        return -(1. / 67813809.) * sqrt (3.)
               * (80383. - 160766. * tau1 - 232909. * tau2 + 152526. * tau2 * tau2 + 305052. * tau2 * tau1)
               * sqrt (22604603.);
      else if (e == 1)
        return (1. / 67813809.) * sqrt (3.)
               * (-2545. + 7014. * tau1 - 4195. * tau2 + 474982. * tau2 * tau2 - 232116. * tau2 * tau1
                  + 9800. * tau1 * tau1)
               * sqrt (22604603.);
      else if (e == 2)
        return (5. / 67813809.) * sqrt (3.)
               * (106401. - 212802. * tau1 - 265819. * tau2 + 159418. * tau2 * tau2 + 318836. * tau2 * tau1)
               * sqrt (22604603.);
      else
        return -(1. / 67813809.) * sqrt (3.)
               * (14269. - 26614. * tau1 - 262925. * tau2 + 716898. * tau2 * tau2 + 251716. * tau2 * tau1
                  + 9800. * tau1 * tau1)
               * sqrt (22604603.);
    case 11:
      if (e == 0)
        return (4. / 2090724953684663.) * sqrt (2.)
               * (302997438. * tau1 * tau1 - 302997438. * tau1 + 98047139. + 235688625. * tau2 * tau2
                  - 293824914. * tau2 + 302997438. * tau2 * tau1)
               * sqrt (10453624768423315.);
      else if (e == 1)
        return (4. / 2090724953684663.) * sqrt (2.)
               * (24456488. * tau1 * tau1 - 9860609. - 5548132. * tau1 - 40311757. * tau2 * tau2 + 236858050. * tau2
                  - 332175154. * tau2 * tau1)
               * sqrt (10453624768423315.);
      else if (e == 2)
        return -(4. / 2090724953684663.) * sqrt (2.)
               * (475844038. * tau1 * tau1 - 475844038. * tau1 + 67323083. + 59252053. * tau2 * tau2 - 127308866. * tau2
                  + 475844038. * tau2 * tau1)
               * sqrt (10453624768423315.);
      else
        return (4. / 2090724953684663.) * sqrt (2.)
               * (24456488. * tau1 * tau1 + 9047747. - 43364844. * tau1 + 316319885. * tau2 * tau2 - 138681948. * tau2
                  + 381088130. * tau2 * tau1)
               * sqrt (10453624768423315.);
    case 12:
      if (e == 0)
        return (1. / 11229548361.
                * (-227207. * tau2 + 402996. * tau2 * tau1 + 201498. * tau2 * tau2 + 25709. - 51418. * tau1))
               * sqrt (4031407861599.);
      else if (e == 1)
        return -(1. / 11229548361.
                 * (374455. * tau2 - 466428. * tau2 * tau1 - 174814. * tau2 * tau2 + 19645. - 96078. * tau1
                    + 88120. * tau1 * tau1))
               * sqrt (4031407861599.);
      else if (e == 2)
        return (5. / 11229548361.
                * (-43663. * tau2 + 51140. * tau2 * tau1 + 25570. * tau2 * tau2 + 18093. - 36186. * tau1))
               * sqrt (4031407861599.);
      else
        return (1. / 11229548361.
                * (-172135. * tau2 + 642668. * tau2 * tau1 + 379734. * tau2 * tau2 + 11687. - 80162. * tau1
                   + 88120. * tau1 * tau1))
               * sqrt (4031407861599.);
    case 13:
      if (e == 0)
        return (1. / 1537694113448282.
                * (-2998716. * tau2 - 6553980. * tau1 + 2549628. * tau2 * tau2 + 6553980. * tau1 * tau1 + 1333543.
                   + 6553980. * tau2 * tau1))
               * sqrt (404363576778211096835.);
      else if (e == 1)
        return (1. / 1537694113448282.
                * (1335668. * tau2 - 7463684. * tau1 - 4154804. * tau2 * tau2 + 4931668. * tau1 * tau1 + 2702051.
                   - 89012. * tau2 * tau1))
               * sqrt (404363576778211096835.);
      else if (e == 2)
        return -(1. / 1537694113448282.
                 * (937117. - 2441404. * tau2 - 64436. * tau1 + 1561004. * tau2 * tau2 + 64436. * tau1 * tau1
                    + 64436. * tau2 * tau1))
               * sqrt (404363576778211096835.);
      else
        return (1. / 1537694113448282.
                * (-1152996. * tau2 - 2399652. * tau1 + 865876. * tau2 * tau2 + 4931668. * tau1 * tau1 + 170035.
                   + 9952348. * tau2 * tau1))
               * sqrt (404363576778211096835.);
    case 14:
      if (e == 0)
        return (6. / 739755559. * (1711. - 5671. * tau2 - 3422. * tau1 + 7920. * tau2 * tau1 + 3960. * tau2 * tau2))
               * sqrt (1083741893935.);
      else if (e == 1)
        return -(2. / 739755559.
                 * (30961. - 56675. * tau2 - 78022. * tau1 + 48072. * tau1 * tau1 + 72192. * tau2 * tau1
                    + 25752. * tau2 * tau2))
               * sqrt (1083741893935.);
      else if (e == 2)
        return (2. / 739755559. * (517. - 1237. * tau2 - 1034. * tau1 + 1440. * tau2 * tau1 + 720. * tau2 * tau2))
               * sqrt (1083741893935.);
      else
        return (2. / 739755559.
                * (1011. - 2605. * tau2 - 18122. * tau1 + 48072. * tau1 * tau1 + 23952. * tau2 * tau1
                   + 1632. * tau2 * tau2))
               * sqrt (1083741893935.);
    case 15:
      if (e == 0)
        return -(1. / 59474.) * sqrt (5.)
               * (-985. - 3420. * tau1 + 12132. * tau2 + 3420. * tau1 * tau1 + 3420. * tau2 * tau1
                  - 22932. * tau2 * tau2)
               * sqrt (29737.);
      else if (e == 1)
        return -(1. / 59474.) * sqrt (5.)
               * (2947. - 7204. * tau1 - 4844. * tau2 + 4340. * tau1 * tau1 + 5804. * tau2 * tau1 + 2492. * tau2 * tau2)
               * sqrt (29737.);
      else if (e == 2)
        return (1. / 59474.) * sqrt (5.)
               * (-3139. - 3412. * tau1 + 8548. * tau2 + 3412. * tau1 * tau1 + 3412. * tau2 * tau1
                  - 5636. * tau2 * tau2)
               * sqrt (29737.);
      else
        return -(1. / 59474.) * sqrt (5.)
               * (83. - 1476. * tau1 - 516. * tau2 + 4340. * tau1 * tau1 + 2876. * tau2 * tau1 + 1028. * tau2 * tau2)
               * sqrt (29737.);
    case 16:
      if (e == 0)
        return (9. / 131.) * sqrt (2.) * sqrt (5.)
               * (122. * tau2 * tau1 + 61. * tau2 * tau2 + 25. - 86. * tau2 - 50. * tau1) * sqrt (131.);
      else if (e == 1)
        return (1. / 131.) * sqrt (2.) * sqrt (5.)
               * (-54. * tau2 * tau1 - 123. * tau2 * tau2 + 96. * tau1 * tau1 + 53. + 74. * tau2 - 146. * tau1)
               * sqrt (131.);
      else if (e == 2)
        return (1. / 131.) * sqrt (2.) * sqrt (5.)
               * (138. * tau2 * tau1 + 49. - 118. * tau2 + 69. * tau2 * tau2 - 98. * tau1) * sqrt (131.);
      else
        return -(1. / 131.) * sqrt (2.) * sqrt (5.)
               * (246. * tau2 * tau1 + 27. * tau2 * tau2 + 96. * tau1 * tau1 + 3. - 26. * tau2 - 46. * tau1)
               * sqrt (131.);
    case 17:
      if (e == 0)
        return -(1. / 227.) * sqrt (2.) * sqrt (3.) * sqrt (5.)
               * (582. * tau2 * tau1 - 258. * tau2 + 33. * tau2 * tau2 + 582. * tau1 * tau1 - 582. * tau1 + 145.)
               * sqrt (227.);
      else if (e == 1)
        return (1. / 227.) * sqrt (2.) * sqrt (3.) * sqrt (5.)
               * (-102. * tau2 * tau1 + 19. + 82. * tau2 - 33. * tau2 * tau2 + 42. * tau1 * tau1 - 58. * tau1)
               * sqrt (227.);
      else if (e == 2)
        return (15. / 227.) * sqrt (2.) * sqrt (3.) * sqrt (5.)
               * (6. * tau2 * tau1 - 2. * tau2 + tau2 * tau2 + 6. * tau1 * tau1 - 6. * tau1 + 1.) * sqrt (227.);
      else
        return (1. / 227.) * sqrt (2.) * sqrt (3.) * sqrt (5.)
               * (186. * tau2 * tau1 + 3. - 46. * tau2 + 111. * tau2 * tau2 + 42. * tau1 * tau1 - 26. * tau1)
               * sqrt (227.);
    default:;
    };
  case 4:
    switch (i) {
    case 0:
      if (e == 0)
        return (1. / 51.) * sqrt (2.) * (-11. - 1176. * tau2 * tau2 + 264. * tau2 + 1344. * tau2 * tau2 * tau2)
               * sqrt (51.);
      else if (e == 1)
        return (1. / 51.) * sqrt (2.) * (-15. - 1512. * tau2 * tau2 + 336. * tau2 + 1792. * tau2 * tau2 * tau2)
               * sqrt (51.);
      else if (e == 2)
        return -(5. / 51.) * sqrt (2.) * (-97. - 840. * tau2 * tau2 + 504. * tau2 + 448. * tau2 * tau2 * tau2)
               * sqrt (51.);
      else
        return (1. / 51.) * sqrt (2.) * (-15. - 1512. * tau2 * tau2 + 336. * tau2 + 1792. * tau2 * tau2 * tau2)
               * sqrt (51.);
    case 1:
      if (e == 0)
        return (12. / 17.) * sqrt (2.) * (2. - 9. * tau2 + 7. * tau2 * tau2 - 4. * tau1 + 14. * tau2 * tau1)
               * sqrt (17.);
      else if (e == 1)
        return (1. / 17.) * sqrt (2.)
               * (-13. + 224. * tau2 - 588. * tau2 * tau2 + 28. * tau1 - 504. * tau2 * tau1
                  + 1344. * tau2 * tau2 * tau1)
               * sqrt (17.);
      else if (e == 2)
        return -(28. / 17.) * sqrt (2.)
               * (-13. + 64. * tau2 - 99. * tau2 * tau2 + 48. * tau2 * tau2 * tau2 + 26. * tau1 - 102. * tau2 * tau1
                  + 96. * tau2 * tau2 * tau1)
               * sqrt (17.);
      else
        return (1. / 17.) * sqrt (2.)
               * (-15. + 308. * tau2 - 1260. * tau2 * tau2 + 1344. * tau2 * tau2 * tau2 + 28. * tau1
                  - 504. * tau2 * tau1 + 1344. * tau2 * tau2 * tau1)
               * sqrt (17.);
    case 2:
      if (e == 0)
        return -(2. / 51.) * sqrt (2.)
               * (318. * tau2 - 594. * tau2 * tau2 + 384. * tau2 * tau2 * tau2 - 780. * tau2 * tau1
                  + 576. * tau2 * tau2 * tau1 + 576. * tau2 * tau1 * tau1 - 55. + 204. * tau1 - 204. * tau1 * tau1)
               * sqrt (255.);
      else if (e == 1)
        return (1. / 51.) * sqrt (2.)
               * (-144. * tau2 + 1284. * tau2 * tau2 - 1024. * tau2 * tau2 * tau2 - 456. * tau2 * tau1
                  - 960. * tau2 * tau2 * tau1 + 768. * tau2 * tau1 * tau1 - 3. + 60. * tau1 - 72. * tau1 * tau1)
               * sqrt (255.);
      else if (e == 2)
        return -(28. / 51.) * sqrt (2.)
               * (30. * tau2 - 39. * tau2 * tau2 + 16. * tau2 * tau2 * tau2 - 138. * tau2 * tau1
                  + 96. * tau2 * tau2 * tau1 + 96. * tau2 * tau1 * tau1 - 7. + 42. * tau1 - 42. * tau1 * tau1)
               * sqrt (255.);
      else
        return (1. / 51.) * sqrt (2.)
               * (252. * tau2 - 828. * tau2 * tau2 + 704. * tau2 * tau2 * tau2 - 1224. * tau2 * tau1
                  + 2496. * tau2 * tau2 * tau1 + 768. * tau2 * tau1 * tau1 - 15. + 84. * tau1 - 72. * tau1 * tau1)
               * sqrt (255.);
    case 3:
      if (e == 0)
        return (8. / 17.) * sqrt (2.) * sqrt (3.)
               * (-41. * tau2 * tau2 + 27. * tau2 + 32. * tau1 - 60. * tau1 * tau1 + 40. * tau1 * tau1 * tau1 - 6.
                  - 82. * tau2 * tau1 + 60. * tau2 * tau1 * tau1 + 60. * tau2 * tau2 * tau1 + 20. * tau2 * tau2 * tau2)
               * sqrt (119.);
      else if (e == 1)
        return -(1. / 51.) * sqrt (2.) * sqrt (3.)
               * (1632. * tau2 * tau2 - 1056. * tau2 - 72. * tau1 - 120. * tau1 * tau1 + 160. * tau1 * tau1 * tau1 + 47.
                  + 2976. * tau2 * tau1 - 1920. * tau2 * tau1 * tau1 - 2496. * tau2 * tau2 * tau1
                  - 640. * tau2 * tau2 * tau2)
               * sqrt (119.);
      else if (e == 2)
        return (64. / 51.) * sqrt (2.) * sqrt (3.)
               * (-3. * tau2 * tau2 + 3. * tau2 + 12. * tau1 - 30. * tau1 * tau1 + 20. * tau1 * tau1 * tau1 - 1.
                  - 24. * tau2 * tau1 + 30. * tau2 * tau1 * tau1 + 12. * tau2 * tau2 * tau1 + tau2 * tau2 * tau2)
               * sqrt (119.);
      else
        return -(1. / 51.) * sqrt (2.) * sqrt (3.)
               * (-360. * tau2 * tau2 + 168. * tau2 + 168. * tau1 - 360. * tau1 * tau1 + 160. * tau1 * tau1 * tau1 - 15.
                  - 1584. * tau2 * tau1 + 2400. * tau2 * tau1 * tau1 + 1824. * tau2 * tau2 * tau1
                  + 224. * tau2 * tau2 * tau2)
               * sqrt (119.);
    case 4:
      if (e == 0)
        return -(6. / 17.) * sqrt (2.) * sqrt (3.)
               * (-44. * tau2 * tau2 + 36. * tau2 + 16. * tau2 * tau2 * tau2 + 40. * tau1 - 40. * tau1 * tau1 - 9.
                  - 120. * tau2 * tau1 + 80. * tau2 * tau1 * tau1 + 80. * tau2 * tau2 * tau1)
               * sqrt (17.);
      else if (e == 1)
        return -(1. / 17.) * sqrt (2.) * sqrt (3.)
               * (-528. * tau2 * tau2 + 704. * tau2 + 128. * tau2 * tau2 * tau2 + 1480. * tau1 - 2280. * tau1 * tau1
                  - 305. - 2400. * tau2 * tau1 + 1920. * tau2 * tau1 * tau1 + 960. * tau2 * tau2 * tau1
                  + 1120. * tau1 * tau1 * tau1)
               * sqrt (17.);
      else if (e == 2)
        return 0.;
      else
        return (1. / 17.) * sqrt (2.) * sqrt (3.)
               * (-72. * tau2 * tau2 + 56. * tau2 + 32. * tau2 * tau2 * tau2 + 280. * tau1 - 1080. * tau1 * tau1 - 15.
                  - 720. * tau2 * tau1 + 1440. * tau2 * tau1 * tau1 + 480. * tau2 * tau2 * tau1
                  + 1120. * tau1 * tau1 * tau1)
               * sqrt (17.);
    case 5:
      if (e == 0)
        return -(2. / 30243. * (50. + 789. * tau2 - 7896. * tau2 * tau2 + 13104. * tau2 * tau2 * tau2)) * sqrt (60486.);
      else if (e == 1)
        return -(2. / 30243. * (-159. + 4551. * tau2 - 25452. * tau2 * tau2 + 36512. * tau2 * tau2 * tau2))
               * sqrt (60486.);
      else if (e == 2)
        return -(2. / 30243. * (-16466. + 75645. * tau2 - 112560. * tau2 * tau2 + 54320. * tau2 * tau2 * tau2))
               * sqrt (60486.);
      else
        return -(2. / 30243. * (-159. + 4551. * tau2 - 25452. * tau2 * tau2 + 36512. * tau2 * tau2 * tau2))
               * sqrt (60486.);
    case 6:
      if (e == 0)
        return -(2. / 9367.
                 * (-515. + 1030. * tau1 + 4281. * tau2 - 9478. * tau2 * tau2 - 7532. * tau2 * tau1
                    + 5712. * tau2 * tau2 * tau2 + 11424. * tau2 * tau2 * tau1))
               * sqrt (28101.);
      else if (e == 1)
        return (1. / 28101.
                * (-471. + 1208. * tau1 + 8338. * tau2 - 9660. * tau2 * tau2 - 27048. * tau2 * tau1
                   - 53312. * tau2 * tau2 * tau2 + 87360. * tau2 * tau2 * tau1))
               * sqrt (28101.);
      else if (e == 2)
        return (14. / 28101.
                * (-3971. + 7942. * tau1 + 17333. * tau2 - 24258. * tau2 * tau2 - 26724. * tau2 * tau1
                   + 10896. * tau2 * tau2 * tau2 + 21792. * tau2 * tau2 * tau1))
               * sqrt (28101.);
      else
        return (1. / 28101.
                * (-737. + 1208. * tau1 + 19918. * tau2 - 104748. * tau2 * tau2 - 27048. * tau2 * tau1
                   + 140672. * tau2 * tau2 * tau2 + 87360. * tau2 * tau2 * tau1))
               * sqrt (28101.);
    case 7:
      if (e == 0)
        return -(2. / 904640832569121.
                 * (-4234020. * tau1 - 35277354. * tau2 * tau2 + 44348976. * tau2 * tau2 * tau2 + 4234020. * tau1 * tau1
                    + 7172928. * tau2 * tau2 * tau1 + 7172928. * tau2 * tau1 * tau1 + 949516. - 2938908. * tau2 * tau1
                    + 5518869. * tau2))
               * sqrt (8892619384154459430.);
      else if (e == 1)
        return (1. / 904640832569121.
                * (3408564. * tau1 + 286030164. * tau2 * tau2 - 172101440. * tau2 * tau2 * tau2 - 5977440. * tau1 * tau1
                   - 310029888. * tau2 * tau2 * tau1 + 77308224. * tau2 * tau1 * tau1 + 857691. + 1095864. * tau2 * tau1
                   - 46352034. * tau2))
               * sqrt (8892619384154459430.);
      else if (e == 2)
        return (2. / 904640832569121.
                * (194216988. * tau1 - 146235138. * tau2 * tau2 + 57097040. * tau2 * tau2 * tau2
                   - 194216988. * tau1 * tau1 + 362631360. * tau2 * tau2 * tau1 + 362631360. * tau2 * tau1 * tau1
                   - 31318235. - 556848348. * tau2 * tau1 + 120405741. * tau2))
               * sqrt (8892619384154459430.);
      else
        return (1. / 904640832569121.
                * (8546316. * tau1 - 1711185. - 185689476. * tau2 * tau2 + 215236672. * tau2 * tau2 * tau2
                   - 5977440. * tau1 * tau1 + 464646336. * tau2 * tau2 * tau1 + 77308224. * tau2 * tau1 * tau1
                   - 167667192. * tau2 * tau1 + 40598370. * tau2))
               * sqrt (8892619384154459430.);
    case 8:
      if (e == 0)
        return (2. / 844880549868557.) * sqrt (3.)
               * (-11048795. + 52777270. * tau1 + 44847573. * tau2 - 92039040. * tau1 * tau1
                  + 61359360. * tau1 * tau1 * tau1 - 128956916. * tau2 * tau1 + 92039040. * tau2 * tau1 * tau1
                  - 59297674. * tau2 * tau2 + 25498896. * tau2 * tau2 * tau2 + 81677472. * tau2 * tau2 * tau1)
               * sqrt (5813623063645540717.);
      else if (e == 1)
        return (1. / 2534641649605671.) * sqrt (3.)
               * (11964793. - 33749584. * tau1 - 200649314. * tau2 + 3306000. * tau1 * tau1
                  + 22568960. * tau1 * tau1 * tau1 + 812678904. * tau2 * tau1 - 648504960. * tau2 * tau1 * tau1
                  - 224511180. * tau2 * tau2 + 364940096. * tau2 * tau2 * tau2 - 55120320. * tau2 * tau2 * tau1)
               * sqrt (5813623063645540717.);
      else if (e == 2)
        return (2. / 2534641649605671.) * sqrt (3.)
               * (-19878199. + 246403438. * tau1 + 58705037. * tau2 - 619941120. * tau1 * tau1
                  + 413294080. * tau1 * tau1 * tau1 - 490947756. * tau2 * tau1 + 619941120. * tau2 * tau1 * tau1
                  - 57578022. * tau2 * tau2 + 18751184. * tau2 * tau2 * tau2 + 244149408. * tau2 * tau2 * tau1)
               * sqrt (5813623063645540717.);
      else
        return (1. / 2534641649605671.) * sqrt (3.)
               * (-4090169. + 40569296. * tau1 + 77044666. * tau2 - 71012880. * tau1 * tau1
                  + 22568960. * tau1 * tau1 * tau1 - 626356776. * tau2 * tau1 + 716211840. * tau2 * tau1 * tau1
                  - 275712396. * tau2 * tau2 + 251013504. * tau2 * tau2 * tau2 + 1309596480. * tau2 * tau2 * tau1)
               * sqrt (5813623063645540717.);
    case 9:
      if (e == 0)
        return -(2. / 10765691411249447639.
                 * (161209100100. * tau1 + 186071273154. * tau2 - 283831332126. * tau2 * tau2
                    + 157399264224. * tau2 * tau2 * tau2 - 461612434260. * tau2 * tau1
                    + 300403334160. * tau2 * tau2 * tau1 - 42240490525. - 161209100100. * tau1 * tau1
                    + 300403334160. * tau2 * tau1 * tau1))
               * sqrt (150719679757492266946.);
      else if (e == 1)
        return (1. / 10765691411249447639.
                * (187982776440. * tau1 - 992699014632. * tau2 + 1708127525124. * tau2 * tau2
                   - 714127013504. * tau2 * tau2 * tau2 + 2439191566680. * tau2 * tau1
                   - 2260670213760. * tau2 * tau2 * tau1 - 13652221827. - 425611235640. * tau1 * tau1
                   + 258043932320. * tau1 * tau1 * tau1 - 1419282286080. * tau2 * tau1 * tau1))
               * sqrt (150719679757492266946.);
      else if (e == 2)
        return -(4. / 633275965367614567.
                 * (-5897464. + 85721730. * tau1 + 19354035. * tau2 - 16054755. * tau2 * tau2
                    + 2194440. * tau2 * tau2 * tau2 - 258888210. * tau2 * tau1 + 173166480. * tau2 * tau2 * tau1
                    - 85721730. * tau1 * tau1 + 173166480. * tau2 * tau1 * tau1))
               * sqrt (150719679757492266946.);
      else
        return -(1. / 10765691411249447639.
                 * (110892102120. * tau1 + 83681836152. * tau2 - 195350878164. * tau2 * tau2
                    + 130783018144. * tau2 * tau2 * tau2 - 1096414128120. * tau2 * tau1
                    + 1352026155360. * tau2 * tau2 * tau1 - 6763251293. - 348520561320. * tau1 * tau1
                    + 258043932320. * tau1 * tau1 * tau1 + 2193414083040. * tau2 * tau1 * tau1))
               * sqrt (150719679757492266946.);
    case 10:
      if (e == 0)
        return -(24. / 920240407198300427.) * sqrt (3.)
               * (127279662. + 204039544. * tau1 - 940664913. * tau2 + 1882607643. * tau2 * tau2
                  - 1375796604. * tau1 * tau1 - 1069222392. * tau2 * tau2 * tau2 + 709572766. * tau2 * tau1
                  - 1679845916. * tau2 * tau2 * tau1 + 1375796604. * tau2 * tau1 * tau1
                  + 917197736. * tau1 * tau1 * tau1)
               * sqrt (4601202035991502135.);
      else if (e == 1)
        return (2. / 2760721221594901281.) * sqrt (3.)
               * (-177394216769. + 763898820394. * tau1 + 426324583460. * tau2 - 336092076180. * tau2 * tau2
                  - 1062952041984. * tau1 * tau1 + 86762894400. * tau2 * tau2 * tau2 - 1261049609064. * tau2 * tau1
                  + 514694054112. * tau2 * tau2 * tau1 + 898559860800. * tau2 * tau1 * tau1
                  + 479569493648. * tau1 * tau1 * tau1)
               * sqrt (4601202035991502135.);
      else if (e == 2)
        return (4. / 2760721221594901281.) * sqrt (3.)
               * (-47762895. + 697977710. * tau1 + 128536021. * tau2 - 110649630. * tau2 * tau2
                  - 1807355760. * tau1 * tau1 + 29876504. * tau2 * tau2 * tau2 - 1366450092. * tau2 * tau1
                  + 662204928. * tau2 * tau2 * tau1 + 1807355760. * tau2 * tau1 * tau1
                  + 1204903840. * tau1 * tau1 * tau1)
               * sqrt (4601202035991502135.);
      else
        return (2. / 2760721221594901281.) * sqrt (3.)
               * (-3122055289. + 76703217370. * tau1 + 12868382174. * tau2 - 18288304356. * tau2 * tau2
                  - 375756438960. * tau1 * tau1 + 8940792560. * tau2 * tau2 * tau2 - 215442765384. * tau2 * tau1
                  + 156282813456. * tau2 * tau2 * tau1 + 540148620144. * tau2 * tau1 * tau1
                  + 479569493648. * tau1 * tau1 * tau1)
               * sqrt (4601202035991502135.);
    case 11:
      if (e == 0)
        return (1. / 905676852258502428609.) * sqrt (2.) * sqrt (3.)
               * (7544141429440. * tau2 * tau2 * tau2 + 144804769200. * tau1 - 144804769200. * tau1 * tau1
                  - 247573720553. + 39483400320. * tau2 * tau1 * tau1 + 39483400320. * tau2 * tau2 * tau1
                  + 2796129345102. * tau2 - 8729779876230. * tau2 * tau2 - 184288169520. * tau2 * tau1)
               * sqrt (100630761362055825401.);
      else if (e == 1)
        return -(1. / 905676852258502428609.) * sqrt (2.) * sqrt (3.)
               * (24266273541440. * tau2 * tau2 * tau2 + 20612007000. * tau1 - 59971373280. * tau1 * tau1
                  + 64071047040. * tau1 * tau1 * tau1 - 77084019997. - 943529408640. * tau2 * tau1 * tau1
                  + 1809294963840. * tau2 * tau2 * tau1 + 2345646250218. * tau2 - 16088530112970. * tau2 * tau2
                  + 597259943280. * tau2 * tau1)
               * sqrt (100630761362055825401.);
      else if (e == 2)
        return (1. / 905676852258502428609.) * sqrt (2.) * sqrt (3.)
               * (37820822392640. * tau2 * tau2 * tau2 + 2971563613200. * tau1 - 2971563613200. * tau1 * tau1
                  - 13847389199837. + 5548562772480. * tau2 * tau1 * tau1 + 5548562772480. * tau2 * tau2 * tau1
                  + 59417823663318. * tau2 - 82988273790870. * tau2 * tau2 - 8520126385680. * tau2 * tau1)
               * sqrt (100630761362055825401.);
      else
        return (1. / 905676852258502428609.) * sqrt (2.) * sqrt (3.)
               * (-21449378121920. * tau2 * tau2 * tau2 + 92882401560. * tau1 - 132241767840. * tau1 * tau1
                  + 64071047040. * tau1 * tau1 * tau1 + 52372339237. + 1135742549760. * tau2 * tau1 * tau1
                  + 3888566922240. * tau2 * tau2 * tau1 - 1906494383298. * tau2 + 12857194507290. * tau2 * tau2
                  - 1554282409680. * tau2 * tau1)
               * sqrt (100630761362055825401.);
    case 12:
      if (e == 0)
        return -(72. / 24280578051280066842136053269.) * sqrt (2.)
               * (-89913208265. * tau2 + 297997576. * tau1 + 346354066891. * tau2 * tau2 - 2597609196. * tau1 * tau1
                  + 1731739464. * tau1 * tau1 * tau1 + 177526804910. * tau2 * tau1 - 512583719676. * tau2 * tau2 * tau1
                  + 2597609196. * tau2 * tau1 * tau1 + 283936078. - 256724794704. * tau2 * tau2 * tau2)
               * sqrt (10381525354495561779357901615959985.);
      else if (e == 1)
        return -(4. / 24280578051280066842136053269.) * sqrt (2.)
               * (742622726374. * tau2 + 169713834857. * tau1 + 3170650107534. * tau2 * tau2 - 3545105382. * tau1 * tau1
                  + 3644377464. * tau1 * tau1 * tau1 - 4311394474968. * tau2 * tau1
                  + 15972451837104. * tau2 * tau2 * tau1 - 51883082160. * tau2 * tau1 * tau1 - 53427686327.
                  - 20656458047232. * tau2 * tau2 * tau2)
               * sqrt (10381525354495561779357901615959985.);
      else if (e == 2)
        return -(4. / 24280578051280066842136053269.) * sqrt (2.)
               * (-77050203333827. * tau2 - 37443894884206. * tau1 + 102499197109260. * tau2 * tau2
                  - 105392523600. * tau1 * tau1 + 70261682400. * tau1 * tau1 * tau1 + 116551119259848. * tau2 * tau1
                  - 88341882435072. * tau2 * tau2 * tau1 + 105392523600. * tau2 * tau1 * tau1 + 18739512862703.
                  - 44188506638136. * tau2 * tau2 * tau2)
               * sqrt (10381525354495561779357901615959985.);
      else
        return -(4. / 24280578051280066842136053269.) * sqrt (2.)
               * (3794211587239. * tau2 + 173556756485. * tau1 - 23565650610936. * tau2 * tau2
                  - 7388027010. * tau1 * tau1 + 3644377464. * tau1 * tau1 * tau1 - 4429936693308. * tau2 * tau1
                  + 16087151133816. * tau2 * tau2 * tau1 + 62816214552. * tau2 * tau1 * tau1 - 116385420612.
                  + 36684437343960. * tau2 * tau2 * tau2)
               * sqrt (10381525354495561779357901615959985.);
    case 13:
      if (e == 0)
        return (4. / 450190343850906324042093431419389.
                * (7723887049488. * tau2 + 5633936712900. * tau1 - 11310827442900. * tau2 * tau2
                   - 5633936712900. * tau1 * tau1 + 4307107230800. * tau2 * tau2 * tau2 - 20941320604140. * tau2 * tau1
                   + 15307383891240. * tau2 * tau2 * tau1 + 15307383891240. * tau2 * tau1 * tau1 - 1511837033851.))
               * sqrt (34560885477513793515221944312683818679790383.);
      else if (e == 1)
        return -(4. / 450190343850906324042093431419389.
                 * (2564480914314. * tau2 - 28384522317. * tau1 - 15633177382020. * tau2 * tau2
                    + 107501634540. * tau1 * tau1 + 6549688499104. * tau2 * tau2 * tau2 - 980844788268. * tau2 * tau1
                    + 19462977054288. * tau2 * tau2 * tau1 - 3349829668320. * tau2 * tau1 * tau1
                    + 62214683880. * tau1 * tau1 * tau1 - 67122037922.))
               * sqrt (34560885477513793515221944312683818679790383.);
      else if (e == 2)
        return -(4. / 450190343850906324042093431419389.
                 * (17036960891259. * tau2 + 30770051855220. * tau1 - 19841738599500. * tau2 * tau2
                    - 30770051855220. * tau1 * tau1 + 7487604817160. * tau2 * tau2 * tau2
                    - 82446561694980. * tau2 * tau1 + 51676509839760. * tau2 * tau2 * tau1
                    + 51676509839760. * tau2 * tau1 * tau1 - 4692563956424.))
               * sqrt (34560885477513793515221944312683818679790383.);
      else
        return (4. / 450190343850906324042093431419389.
                * (2139456340677. * tau2 + 373262798403. * tau1 - 11804449483356. * tau2 * tau2
                   - 294145686180. * tau1 * tau1 + 16325332907384. * tau2 * tau2 * tau2 - 8268795497268. * tau2 * tau1
                   + 26349280442568. * tau2 * tau2 * tau1 + 3536473719960. * tau2 * tau1 * tau1
                   + 62214683880. * tau1 * tau1 * tau1 - 74209758181.))
               * sqrt (34560885477513793515221944312683818679790383.);
    case 14:
      if (e == 0)
        return (8. / 6380720215578715256716408976310607022277.
                * (-865666978019157. * tau2 - 55556560994562. * tau1 - 302892267120612. * tau1 * tau1
                   - 1312841251683268. * tau2 * tau2 * tau2 + 201928178080408. * tau1 * tau1 * tau1
                   + 1372885127923140. * tau2 * tau1 - 2524718414326332. * tau2 * tau2 * tau1
                   + 302892267120612. * tau2 * tau1 * tau1 + 78260325017383. + 2100247904685042. * tau2 * tau2))
               * sqrt (38697087397201085046824763022410609242891822943821145.);
      else if (e == 1)
        return -(4. / 6380720215578715256716408976310607022277.
                 * (-522961534414587. * tau2 - 284427010918203. * tau1 + 168637635998520. * tau1 * tau1
                    + 6088129267565552. * tau2 * tau2 * tau2 + 76974385176312. * tau1 * tau1 * tau1
                    + 4891309274804220. * tau2 * tau1 + 6079189293644784. * tau2 * tau2 * tau1
                    - 4946764497496536. * tau2 * tau1 * tau1 + 72312277135198. - 7310482464324780. * tau2 * tau2))
               * sqrt (38697087397201085046824763022410609242891822943821145.);
      else if (e == 2)
        return (4. / 6380720215578715256716408976310607022277.
                * (1444448045208459. * tau2 + 7586203235228550. * tau1 - 19591382245714032. * tau1 * tau1
                   + 359309007353848. * tau2 * tau2 * tau2 + 13060921497142688. * tau1 * tau1 * tau1
                   - 14894075100902400. * tau2 * tau1 + 7249078763279040. * tau2 * tau2 * tau1
                   + 19591382245714032. * tau2 * tau1 * tau1 - 527871243328603. - 1275885809233704. * tau2 * tau2))
               * sqrt (38697087397201085046824763022410609242891822943821145.);
      else
        return -(4. / 6380720215578715256716408976310607022277.
                 * (862188173714676. * tau2 + 283771416607773. * tau1 - 399560791527456. * tau1 * tau1
                    + 5014798908752080. * tau2 * tau2 * tau2 + 76974385176312. * tau1 * tau1 * tau1
                    - 5801341303243764. * tau2 * tau1 + 16203641444166792. * tau2 * tau2 * tau1
                    + 5177687653025472. * tau2 * tau1 * tau1 - 33497287391827. - 4170487341036312. * tau2 * tau2))
               * sqrt (38697087397201085046824763022410609242891822943821145.);
    case 15:
      if (e == 0)
        return -(4. / 35516731602837412764420243765174837.
                 * (-12837488172413. - 123415494481620. * tau2 * tau1 + 73509208386120. * tau2 * tau1 * tau1
                    + 73509208386120. * tau2 * tau2 * tau1 + 49906286095500. * tau1 + 42690194479914. * tau2
                    - 49906286095500. * tau1 * tau1 - 29319191117400. * tau2 * tau2
                    - 16354402687400. * tau2 * tau2 * tau2))
               * sqrt (4587144914088926523642597305294551825737974887.);
      else if (e == 1)
        return -(4. / 35516731602837412764420243765174837.
                 * (1863564729991. + 293378451286824. * tau2 * tau1 - 192708078064200. * tau2 * tau1 * tau1
                    - 132611304892464. * tau2 * tau2 * tau1 + 14128387753800. * tau1 * tau1 * tau1
                    + 2683151627571. * tau1 - 103745743940397. * tau2 - 18078506976780. * tau1 * tau1
                    + 60644719016880. * tau2 * tau2 + 33991136549968. * tau2 * tau2 * tau2))
               * sqrt (4587144914088926523642597305294551825737974887.);
      else if (e == 2)
        return (4. / 35516731602837412764420243765174837.
                * (-3855190312. - 91426887169020. * tau2 * tau1 + 55766887370640. * tau2 * tau1 * tau1
                   + 55766887370640. * tau2 * tau2 * tau1 + 35659999798380. * tau1 - 3794741004603. * tau2
                   - 35659999798380. * tau1 * tau1 + 9706701686460. * tau2 * tau2
                   - 6040814783960. * tau2 * tau2 * tau2))
               * sqrt (4587144914088926523642597305294551825737974887.);
      else
        return (4. / 35516731602837412764420243765174837.
                * (-596597134582. - 140651017410816. * tau2 * tau1 + 235093241325600. * tau2 * tau1 * tau1
                   + 295190014497336. * tau2 * tau2 * tau1 + 14128387753800. * tau1 * tau1 * tau1
                   + 8911300935411. * tau1 + 11986671653184. * tau2 - 24306656284620. * tau1 * tau1
                   - 44377775250612. * tau2 * tau2 + 40234024375568. * tau2 * tau2 * tau2))
               * sqrt (4587144914088926523642597305294551825737974887.);
    case 16:
      if (e == 0)
        return -(2. / 4429162483638155085928905495293501.
                 * (-902027828018790. * tau2 * tau1 + 1114817925578880. * tau2 * tau1 * tau1
                    + 419714427181440. * tau2 * tau2 * tau1 - 1114817925578880. * tau1 * tau1
                    - 103462164810675. * tau2 * tau2 + 24054225994240. * tau2 * tau2 * tau2 + 467843784396552. * tau1
                    + 743211950385920. * tau1 * tau1 * tau1 - 48118904601796. + 127526843418231. * tau2))
               * sqrt (1977693186713807260499140827732350041750787.);
      else if (e == 1)
        return -(2. / 4429162483638155085928905495293501.
                 * (3861292811332650. * tau2 * tau1 - 2047919222304000. * tau2 * tau1 * tau1
                    - 4240917078976800. * tau2 * tau2 * tau1 - 1446732052196520. * tau1 * tau1
                    + 3506259687314565. * tau2 * tau2 - 1619202493724960. * tau2 * tau2 * tau2 + 892404533563320. * tau1
                    + 726050022346080. * tau1 * tau1 * tau1 - 163614571731043. - 1736827118543049. * tau2))
               * sqrt (1977693186713807260499140827732350041750787.);
      else if (e == 2)
        return (2. / 4429162483638155085928905495293501.
                * (589738585996. - 116842040777610. * tau2 * tau1 + 204619294393440. * tau2 * tau1 * tau1
                   + 46139163828000. * tau2 * tau2 * tau1 + 67026954292488. * tau1 - 10375149661671. * tau2
                   - 204619294393440. * tau1 * tau1 + 20819044893915. * tau2 * tau2
                   - 11033633818240. * tau2 * tau2 * tau2 + 136412862928960. * tau1 * tau1 * tau1))
               * sqrt (1977693186713807260499140827732350041750787.);
      else
        return -(2. / 4429162483638155085928905495293501.
                 * (-1697381662958790. * tau2 * tau1 + 4226069289342240. * tau2 * tau1 * tau1
                    + 2033071432669440. * tau2 * tau2 * tau1 - 731418014841720. * tau1 * tau1
                    - 231306256454835. * tau2 * tau2 + 152254659398240. * tau2 * tau2 * tau2 + 177090496208520. * tau1
                    + 726050022346080. * tau1 * tau1 * tau1 - 8107931981837. + 100544025722919. * tau2))
               * sqrt (1977693186713807260499140827732350041750787.);
    case 17:
      if (e == 0)
        return -(2. / 2347070638132331273272455.) * sqrt (5.)
               * (-5396852251249. - 81879933691710. * tau2 * tau1 + 56363129182080. * tau2 * tau1 * tau1
                  + 56363129182080. * tau2 * tau2 * tau1 + 25516804509630. * tau1 + 20019115103835. * tau2
                  - 25516804509630. * tau1 * tau1 - 17380798498455. * tau2 * tau2 - 1108238286400. * tau2 * tau2 * tau2)
               * sqrt (9023182675486518006136327.);
      else if (e == 1)
        return (2. / 2347070638132331273272455.) * sqrt (5.)
               * (-109424567344688. - 756000640994130. * tau2 * tau1 + 507350406064320. * tau2 * tau1 * tau1
                  + 318653554050720. * tau2 * tau2 * tau1 + 443580065894538. * tau1 + 274131673456605. * tau2
                  - 585513182972970. * tau1 * tau1 - 225817070108385. * tau2 * tau2
                  + 61185468551200. * tau2 * tau2 * tau2 + 252343928010720. * tau1 * tau1 * tau1)
               * sqrt (9023182675486518006136327.);
      else if (e == 2)
        return -(14. / 2347070638132331273272455.) * sqrt (5.)
               * (-429235621. - 726587567010. * tau2 * tau1 + 437189459040. * tau2 * tau1 * tau1
                  + 437189459040. * tau2 * tau2 * tau1 + 289398107970. * tau1 - 31687044699. * tau2
                  - 289398107970. * tau1 * tau1 + 81609115335. * tau2 * tau2 - 50560943360. * tau2 * tau2 * tau2)
               * sqrt (9023182675486518006136327.);
      else
        return -(2. / 2347070638132331273272455.) * sqrt (5.)
               * (-986243587600. - 84337030983870. * tau2 * tau1 + 249681377967840. * tau2 * tau1 * tau1
                  + 60984525954240. * tau2 * tau2 * tau1 + 29585483980758. * tau1 + 4104045453963. * tau2
                  - 171518601059190. * tau1 * tau1 - 5654913867015. * tau2 * tau2 + 2461607445920. * tau2 * tau2 * tau2
                  + 252343928010720. * tau1 * tau1 * tau1)
               * sqrt (9023182675486518006136327.);
    case 18:
      if (e == 0)
        return (1. / 4685971895550392990415.) * sqrt (2.)
               * (15979364504400. * tau2 * tau2 * tau2 - 943114772473. + 6026917395480. * tau2
                  - 16185801835185. * tau2 * tau2 + 2851410481710. * tau1 - 2851410481710. * tau1 * tau1
                  - 7853829170670. * tau2 * tau1 + 5002418688960. * tau2 * tau2 * tau1
                  + 5002418688960. * tau2 * tau1 * tau1)
               * sqrt (32801803268852750932905.);
      else if (e == 1)
        return (1. / 4685971895550392990415.) * sqrt (2.)
               * (-8590475828400. * tau2 * tau2 * tau2 - 27199413019. - 4934866738800. * tau2
                  + 10439985152895. * tau2 * tau2 + 657180239154. * tau1 - 1432137090510. * tau1 * tau1
                  + 825312928160. * tau1 * tau1 * tau1 + 10428888842610. * tau2 * tau1
                  - 9003435769440. * tau2 * tau2 * tau1 - 5851383884640. * tau2 * tau1 * tau1)
               * sqrt (32801803268852750932905.);
      else if (e == 2)
        return -(1. / 4685971895550392990415.) * sqrt (2.)
               * (12776635041040. * tau2 * tau2 * tau2 - 5500681931591. + 22207156652916. * tau2
                  - 29395021886415. * tau2 * tau2 + 3307928504370. * tau1 - 3307928504370. * tau1 * tau1
                  - 8256801806610. * tau2 * tau1 + 4948873302240. * tau2 * tau2 * tau1
                  + 4948873302240. * tau2 * tau1 * tau1)
               * sqrt (32801803268852750932905.);
      else
        return -(1. / 4685971895550392990415.) * sqrt (2.)
               * (6263736871760. * tau2 * tau2 * tau2 - 23156663785. + 626206623444. * tau2
                  - 3754230004095. * tau2 * tau2 + 268844842614. * tau1 - 1043801693970. * tau1 * tau1
                  + 825312928160. * tau1 * tau1 * tau1 - 3361482314610. * tau2 * tau1
                  + 5175270784320. * tau2 * tau2 * tau1 + 8327322669120. * tau2 * tau1 * tau1)
               * sqrt (32801803268852750932905.);
    case 19:
      if (e == 0)
        return (1. / 224288466335587227020110639.) * sqrt (2.)
               * (26401542816284. * tau2 + 7569453530758. * tau1 - 53912941735305. * tau2 * tau2
                  - 447340032960. * tau1 * tau1 + 31221569012240. * tau2 * tau2 * tau2
                  + 298226688640. * tau1 * tau1 * tau1 - 45680972134770. * tau2 * tau1
                  + 62592251368800. * tau2 * tau2 * tau1 + 447340032960. * tau2 * tau1 * tau1 - 3710170093219.)
               * sqrt (3810436754575291399844659645971.);
      else if (e == 1)
        return -(1. / 672865399006761681060331917.) * sqrt (2.)
               * (-16828358446012. * tau2 + 2703515000770. * tau1 - 1950569881365. * tau2 * tau2
                  - 6315542069040. * tau1 * tau1 + 86557884476880. * tau2 * tau2 * tau2
                  + 3602896084320. * tau1 * tau1 * tau1 + 51120651827190. * tau2 * tau1
                  - 70486646544960. * tau2 * tau2 * tau1 - 24140323097760. * tau2 * tau1 * tau1 - 199629664249.)
               * sqrt (3810436754575291399844659645971.);
      else if (e == 2)
        return (1. / 672865399006761681060331917.) * sqrt (2.)
               * (334299298799272. * tau2 + 200999473364474. * tau1 - 424257177681195. * tau2 * tau2
                  - 86599091838240. * tau1 * tau1 + 176024433591120. * tau2 * tau2 * tau2
                  + 57732727892160. * tau1 * tau1 * tau1 - 554198216072310. * tau2 * tau1
                  + 380915231128320. * tau2 * tau2 * tau1 + 86599091838240. * tau2 * tau1 * tau1 - 86066554709197.)
               * sqrt (3810436754575291399844659645971.);
      else
        return -(1. / 672865399006761681060331917.) * sqrt (2.)
               * (-9270851167768. * tau2 + 881119115650. * tau1 + 70784075874075. * tau2 * tau2
                  - 4493146183920. * tau1 * tau1 - 129301311839760. * tau2 * tau2 * tau2
                  + 3602896084320. * tau1 * tau1 * tau1 - 6146286736170. * tau2 * tau1
                  - 11397312096480. * tau2 * tau2 * tau1 + 34949011350720. * tau2 * tau1 * tau1 + 208760648199.)
               * sqrt (3810436754575291399844659645971.);
    case 20:
      if (e == 0)
        return (4. / 4891272024616371501303055622019.) * sqrt (2.)
               * (2916945992442. * tau1 - 2916945992442. * tau1 * tau1 - 1263543266489. + 11203578173532. * tau2
                  - 32224956773007. * tau2 * tau2 + 30463180770552. * tau2 * tau2 * tau2 - 12641311089594. * tau2 * tau1
                  + 9724365097152. * tau2 * tau2 * tau1 + 9724365097152. * tau2 * tau1 * tau1)
               * sqrt (338531287325211031860850338067150911335.);
      else if (e == 1)
        return (4. / 4891272024616371501303055622019.) * sqrt (2.)
               * (1151007682716. * tau1 - 2943675525354. * tau1 * tau1 + 1690378577200. * tau1 * tau1 * tau1
                  + 61957824709. - 12956432108676. * tau2 + 29505156081585. * tau2 * tau2
                  - 3451882567800. * tau2 * tau2 * tau2 + 26852854261950. * tau2 * tau1
                  - 42974094337392. * tau2 * tau2 * tau1 - 11578107511848. * tau2 * tau1 * tau1)
               * sqrt (338531287325211031860850338067150911335.);
      else if (e == 2)
        return (4. / 4891272024616371501303055622019.) * sqrt (2.)
               * (51646435564170. * tau1 - 51646435564170. * tau1 * tau1 - 8760244672247. + 31302616329774. * tau2
                  - 36358388734263. * tau2 * tau2 + 13817755030144. * tau2 * tau2 * tau2
                  - 132558664209354. * tau2 * tau1 + 80912228645184. * tau2 * tau2 * tau1
                  + 80912228645184. * tau2 * tau1 * tau1)
               * sqrt (338531287325211031860850338067150911335.);
      else
        return -(4. / 4891272024616371501303055622019.) * sqrt (2.)
               * (334792363608. * tau1 - 2127460206246. * tau1 * tau1 + 1690378577200. * tau1 * tau1 * tau1
                  + 40331440729. - 1983522277818. * tau2 + 15038117287815. * tau2 * tau2
                  - 26253725680544. * tau2 * tau2 * tau2 - 558281174238. * tau2 * tau1
                  - 14746743582096. * tau2 * tau2 * tau1 + 16649243243448. * tau2 * tau1 * tau1)
               * sqrt (338531287325211031860850338067150911335.);
    case 21:
      if (e == 0)
        return (4. / 8010895569642034387093225507772921797.) * sqrt (2.) * sqrt (3.)
               * (-544074453531273. * tau2 * tau2 + 430223607436900. * tau2 + 222581252051144. * tau2 * tau2 * tau2
                  + 511981397866982. * tau1 - 108730405956771. - 1232027574867138. * tau2 * tau1
                  + 739683090055728. * tau2 * tau2 * tau1 - 883561757860320. * tau1 * tau1
                  + 589041171906880. * tau1 * tau1 * tau1 + 883561757860320. * tau2 * tau1 * tau1)
               * sqrt (2206416960936351241282951190627144110360471628935.);
      else if (e == 1)
        return -(4. / 24032686708926103161279676523318765391.) * sqrt (2) * sqrt (3)
               * (-454480788017877. * tau2 * tau2 - 325425903991988. * tau2 + 418098131421528. * tau2 * tau2 * tau2
                  + 31792165898948. * tau1 + 633733829767. + 1128878245892526. * tau2 * tau1
                  + 339357636579840. * tau2 * tau2 * tau1 - 98098390802376. * tau1 * tau1
                  + 69600670060224. * tau1 * tau1 * tau1 - 868498944007416. * tau2 * tau1 * tau1)
               * sqrt (2206416960936351241282951190627144110360471628935.);
      else if (e == 2)
        return -(4. / 24032686708926103161279676523318765391.) * sqrt (2.) * sqrt (3.)
               * (145211743860843. * tau2 * tau2 - 76312238624402. * tau2 - 74523035155968. * tau2 * tau2 * tau2
                  + 674082719974322. * tau1 + 5623529919527. - 1229282142217002. * tau2 * tau1
                  + 536283709501440. * tau2 * tau2 * tau1 - 2055989339440128. * tau1 * tau1
                  + 1370659559626752. * tau1 * tau1 * tau1 + 2055989339440128. * tau2 * tau1 * tau1)
               * sqrt (2206416960936351241282951190627144110360471628935.);
      else
        return -(4. / 24032686708926103161279676523318765391.) * sqrt (2.) * sqrt (3.)
               * (-603700110062565. * tau2 * tau2 + 109443996581746. * tau2 + 859359119225952. * tau2 * tau2 * tau2
                  + 44397394474868. * tau1 - 3928178986563. - 829526880878898. * tau2 * tau1
                  + 2285157534775344. * tau2 * tau2 * tau1 - 110703619378296. * tau1 * tau1
                  + 69600670060224. * tau1 * tau1 * tau1 + 1077300954188088. * tau2 * tau1 * tau1)
               * sqrt (2206416960936351241282951190627144110360471628935.);
    case 22:
      if (e == 0)
        return (1. / 36794210061790373790848509118616973701.
                * (-810501145781685. * tau2 * tau2 + 547993219860726. * tau2 + 339256701212520. * tau2 * tau2 * tau2
                   + 399493881998430. * tau1 - 399493881998430. * tau1 * tau1 + 1204911521342400. * tau2 * tau1 * tau1
                   + 1204911521342400. * tau2 * tau2 * tau1 - 104211082717321. - 1604405403340830. * tau2 * tau1))
               * sqrt (2170343179656283287883941488614423306932736027111121.);
      else if (e == 1)
        return (1. / 36794210061790373790848509118616973701.
                * (-2433875842005. * tau2 * tau2 - 239334306942864. * tau2 + 107884877384280. * tau2 * tau2 * tau2
                   + 25050161456406. * tau1 - 65328395848230. * tau1 * tau1 + 42242443327040. * tau1 * tau1 * tau1
                   - 475351055700720. * tau2 * tau1 * tau1 - 93080754512160. * tau2 * tau2 * tau1 - 426028325995.
                   + 694645326664650. * tau2 * tau1))
               * sqrt (2170343179656283287883941488614423306932736027111121.);
      else if (e == 2)
        return (1. / 36794210061790373790848509118616973701.
                * (-226081054029765. * tau2 * tau2 + 188089757677170. * tau2 + 89659331725880. * tau2 * tau2 * tau2
                   + 193907952084750. * tau1 - 193907952084750. * tau1 * tau1 - 51450006715717.
                   - 481676057880750. * tau2 * tau1 + 287768105796000. * tau2 * tau1 * tau1
                   + 287768105796000. * tau2 * tau2 * tau1))
               * sqrt (2170343179656283287883941488614423306932736027111121.);
      else
        return -(1. / 36794210061790373790848509118616973701.
                 * (-221941088515515. * tau2 * tau2 + 41160735720000. * tau2 + 316627867131320. * tau2 * tau2 * tau2
                    + 21120699741066. * tau1 - 61398934132890. * tau1 * tau1 + 42242443327040. * tau1 * tau1 * tau1
                    + 602078385681840. * tau2 * tau1 * tau1 + 984348686870400. * tau2 * tau2 * tau1 - 1538180609221.
                    - 378854653002570. * tau2 * tau1))
               * sqrt (2170343179656283287883941488614423306932736027111121.);
    case 23:
      if (e == 0)
        return -(1. / 473145761369413066487006334443891328327.) * sqrt (3.)
               * (454765962881490. * tau2 + 593669146897674. * tau1 - 1153945315866048. * tau1 * tau1
                  + 769296877244032. * tau1 * tau1 * tau1 - 104510354137829. - 1469808094731354. * tau2 * tau1
                  + 1153945315866048. * tau2 * tau1 * tau1 - 728103048224517. * tau2 * tau2
                  + 377847439480856. * tau2 * tau2 * tau2 + 1140343317583728. * tau2 * tau2 * tau1)
               * sqrt (115935593470518835817766866162978250334523437784363455.);
      else if (e == 1)
        return (1. / 473145761369413066487006334443891328327.) * sqrt (3.)
               * (-374908740803664. * tau2 + 142540313543538. * tau1 - 227459974019856. * tau1 * tau1
                  + 113528693807904. * tau1 * tau1 * tau1 - 27349684692337. + 935860201978914. * tau2 * tau1
                  - 554277230276784. * tau2 * tau1 * tau1 + 453425348549205. * tau2 * tau2
                  - 23324339238248. * tau2 * tau2 * tau2 - 683569706048304. * tau2 * tau2 * tau1)
               * sqrt (115935593470518835817766866162978250334523437784363455.);
      else if (e == 2)
        return -(1. / 473145761369413066487006334443891328327.) * sqrt (3.)
               * (-252914338292682. * tau2 - 117930052625598. * tau1 - 44709478109664. * tau1 * tau1
                  + 29806318739776. * tau1 * tau1 * tau1 + 66416605997743. + 343189145850102. * tau2 * tau1
                  + 44709478109664. * tau2 * tau1 * tau1 + 315004838691411. * tau2 * tau2
                  - 128507106396472. * tau2 * tau2 * tau2 - 242111053423056. * tau2 * tau2 * tau1)
               * sqrt (115935593470518835817766866162978250334523437784363455.);
      else
        return (1. / 473145761369413066487006334443891328327.) * sqrt (3.)
               * (21532216029072. * tau2 + 28206446927538. * tau1 - 113126107403856. * tau1 * tau1
                  + 113528693807904. * tau1 * tau1 * tau1 - 1259348639249. - 398946473382366. * tau2 * tau1
                  + 894863311700496. * tau2 * tau1 * tau1 - 55676008479411. * tau2 * tau2
                  + 7560557274632. * tau2 * tau2 * tau2 + 765570835928976. * tau2 * tau2 * tau1)
               * sqrt (115935593470518835817766866162978250334523437784363455.);
    case 24:
      if (e == 0)
        return -(1. / 12008131465057998535947.) * sqrt (5.)
               * (-1236464420742. * tau2 * tau1 + 844896575904. * tau2 * tau1 * tau1
                  + 844896575904. * tau2 * tau2 * tau1 - 92201316317. + 377683669254. * tau2
                  - 473373785937. * tau2 * tau2 + 195277801704. * tau2 * tau2 * tau2 + 391567844838. * tau1
                  - 391567844838. * tau1 * tau1)
               * sqrt (77136233821044229928744879.);
      else if (e == 1)
        return -(1. / 12008131465057998535947.) * sqrt (5.)
               * (8953699554. * tau2 * tau1 + 64684666032. * tau2 * tau1 * tau1 - 469802454528. * tau2 * tau2 * tau1
                  - 118869058463. - 51848336400. * tau2 + 438289085391. * tau2 * tau2
                  - 266962421928. * tau2 * tau2 * tau2 + 484688476422. * tau1 - 642448007502. * tau1 * tau1
                  + 277725306880. * tau1 * tau1 * tau1)
               * sqrt (77136233821044229928744879.);
      else if (e == 2)
        return (1. / 12008131465057998535947.) * sqrt (5.)
               * (-45938028330. * tau2 * tau1 + 27915256128. * tau2 * tau1 * tau1 + 27915256128. * tau2 * tau2 * tau1
                  + 21460463777. - 85872627018. * tau2 + 112800364737. * tau2 * tau2 - 48666996472. * tau2 * tau2 * tau2
                  + 18022772202. * tau1 - 18022772202. * tau1 * tau1)
               * sqrt (77136233821044229928744879.);
      else
        return (1. / 12008131465057998535947.) * sqrt (5.)
               * (-243132794658. * tau2 * tau1 + 768491254608. * tau2 * tau1 * tau1 + 234004134048. * tau2 * tau2 * tau1
                  - 1096717337. + 11178352872. * tau2 - 20891512383. * tau2 * tau2 + 10200608248. * tau2 * tau2 * tau2
                  + 32968382058. * tau1 - 190727913138. * tau1 * tau1 + 277725306880. * tau1 * tau1 * tau1)
               * sqrt (77136233821044229928744879.);
    case 25:
      if (e == 0)
        return (1. / 1555339646295669787201.
                * (161849280390. * tau2 * tau1 + 42195622944. * tau2 * tau1 * tau1 - 238651009584. * tau2 * tau2 * tau1
                   + 221347956459. * tau2 * tau2 - 126358108616. * tau2 * tau2 * tau2 - 42195622944. * tau1 * tau1
                   - 18838110822. * tau1 + 16451659235. - 111441507078. * tau2 + 28130415296. * tau1 * tau1 * tau1))
               * sqrt (9902847527964529535108767.);
      else if (e == 1)
        return (1. / 1555339646295669787201.
                * (-1650197576610. * tau2 * tau1 + 1049759965968. * tau2 * tau1 * tau1
                   + 752444691024. * tau2 * tau2 * tau1 - 574802546469. * tau2 * tau2
                   + 171461098696. * tau2 * tau2 * tau2 - 1135971347472. * tau1 * tau1 + 899615172742. * tau1
                   - 233516057775. + 636822295568. * tau2 + 471006069600. * tau1 * tau1 * tau1))
               * sqrt (9902847527964529535108767.);
      else if (e == 2)
        return -(1. / 1555339646295669787201.
                 * (5001081450. * tau2 * tau1 + 7244473344. * tau2 * tau1 * tau1 - 4393240272. * tau2 * tau2 * tau1
                    + 1726856865. - 6642222038. * tau2 + 8319397533. * tau2 * tau2 - 3404032360. * tau2 * tau2 * tau2
                    - 1038889282. * tau1 - 7244473344. * tau1 * tau1 + 4829648896. * tau1 * tau1 * tau1))
               * sqrt (9902847527964529535108767.);
      else
        return (1. / 1555339646295669787201.
                * (-104771367330. * tau2 * tau1 + 363258242832. * tau2 * tau1 * tau1 + 65942967888. * tau2 * tau2 * tau1
                   - 5366650557. * tau2 * tau2 + 2229695960. * tau2 * tau2 * tau2 - 277046861328. * tau1 * tau1
                   + 40690686598. * tau1 - 1133837095. + 4306001672. * tau2 + 471006069600. * tau1 * tau1 * tau1))
               * sqrt (9902847527964529535108767.);
    case 26:
      if (e == 0)
        return -(2. / 2676102628445301.) * sqrt (2.) * sqrt (5.)
               * (-79707420. * tau2 * tau1 + 59898912. * tau2 * tau1 * tau1 + 59898912. * tau2 * tau2 * tau1 - 3337763.
                  + 2836659. * tau2 + 56291073. * tau2 * tau2 + 19808508. * tau1 - 19808508. * tau1 * tau1
                  - 99945333. * tau2 * tau2 * tau2)
               * sqrt (14155690870266160523.);
      else if (e == 1)
        return -(2. / 2676102628445301.) * sqrt (2.) * sqrt (5.)
               * (13981120. * tau1 * tau1 * tau1 - 38797740. * tau2 * tau1 + 23864640. * tau2 * tau1 * tau1
                  + 18024864. * tau2 * tau2 * tau1 - 7280351. + 15723903. * tau2 - 15723903. * tau2 * tau2
                  + 27473508. * tau1 - 34141452. * tau1 * tau1 + 7117227. * tau2 * tau2 * tau2)
               * sqrt (14155690870266160523.);
      else if (e == 2)
        return -(2. / 2676102628445301.) * sqrt (2.) * sqrt (5.)
               * (-68969916. * tau2 * tau1 + 41383968. * tau2 * tau1 * tau1 + 41383968. * tau2 * tau2 * tau1 + 8535709.
                  - 35729901. * tau2 + 48771201. * tau2 * tau2 + 27585948. * tau1 - 27585948. * tau1 * tau1
                  - 21726677. * tau2 * tau2 * tau2)
               * sqrt (14155690870266160523.);
      else
        return (2. / 2676102628445301.) * sqrt (2.) * sqrt (5.)
               * (13981120. * tau1 * tau1 * tau1 - 6672276. * tau2 * tau1 + 18078720. * tau2 * tau1 * tau1
                  + 12238944. * tau2 * tau2 * tau1 - 32825. + 343161. * tau2 - 1171329. * tau2 * tau2 + 1133964. * tau1
                  - 7801908. * tau1 * tau1 + 1024117. * tau2 * tau2 * tau2)
               * sqrt (14155690870266160523.);
    case 27:
      if (e == 0)
        return -(4. / 481435005651514309593.) * sqrt (7.)
               * (128335460. - 1019851173. * tau2 - 68609487. * tau1 + 2210359119. * tau2 * tau2
                  - 1318843406. * tau2 * tau2 * tau2 - 564184299. * tau1 * tau1 + 376122866. * tau1 * tau1 * tau1
                  - 2449625379. * tau2 * tau2 * tau1 + 564184299. * tau2 * tau1 * tau1 + 1406908560. * tau2 * tau1)
               * sqrt (2111308236157200531699814841.);
      else if (e == 1)
        return (4. / 481435005651514309593.) * sqrt (7.)
               * (-67398760. - 118444963. * tau2 + 274520323. * tau1 + 386770449. * tau2 * tau2
                  - 210605746. * tau2 * tau2 * tau2 - 362797173. * tau1 * tau1 + 156249390. * tau1 * tau1 * tau1
                  - 422918589. * tau2 * tau2 * tau1 - 83614923. * tau2 * tau1 * tau1 + 211588080. * tau2 * tau1)
               * sqrt (2111308236157200531699814841.);
      else if (e == 2)
        return (4. / 481435005651514309593.) * sqrt (7.)
               * (62426000. - 241814863. * tau2 - 53115257. * tau1 + 304845393. * tau2 * tau2
                  - 125456530. * tau2 * tau2 * tau2 - 215210229. * tau1 * tau1 + 143473486. * tau1 * tau1 * tau1
                  - 179176317. * tau2 * tau2 * tau1 + 215210229. * tau2 * tau1 * tau1 + 215304240. * tau2 * tau1)
               * sqrt (2111308236157200531699814841.);
      else
        return (4. / 481435005651514309593.) * sqrt (7.)
               * (-573780. + 8145953. * tau2 + 17674147. * tau1 - 25444623. * tau2 * tau2
                  + 27551470. * tau2 * tau2 * tau2 - 105950997. * tau1 * tau1 + 156249390. * tau1 * tau1 * tau1
                  + 213059427. * tau2 * tau2 * tau1 + 552363093. * tau2 * tau1 * tau1 - 167543760. * tau2 * tau1)
               * sqrt (2111308236157200531699814841.);
    case 28:
      if (e == 0)
        return (2. / 37417.) * sqrt (2.) * sqrt (7.)
               * (4929. * tau2 * tau2 * tau2 + 39174. * tau2 * tau2 * tau1 + 39174. * tau2 * tau1 * tau1 + 17214. * tau1
                  - 17214. * tau1 * tau1 - 4301. + 18051. * tau2 - 20967. * tau2 * tau2 - 56388. * tau2 * tau1)
               * sqrt (187085.);
      else if (e == 1)
        return -(2. / 37417.) * sqrt (2.) * sqrt (7.)
               * (351. * tau2 * tau2 * tau2 + 1280. * tau1 * tau1 * tau1 - 4038. * tau2 * tau2 * tau1
                  - 3750. * tau2 * tau1 * tau1 + 1986. * tau1 - 2802. * tau1 * tau1 - 455. - 2535. * tau2
                  + 2535. * tau2 * tau2 + 6324. * tau2 * tau1)
               * sqrt (187085.);
      else if (e == 2)
        return (26. / 37417.) * sqrt (2.) * sqrt (7.)
               * (-25. + 87. * tau2 + 150. * tau1 - 99. * tau2 * tau2 + 37. * tau2 * tau2 * tau2 - 150. * tau1 * tau1
                  + 222. * tau2 * tau2 * tau1 + 222. * tau2 * tau1 * tau1 - 372. * tau2 * tau1)
               * sqrt (187085.);
      else
        return (2. / 37417.) * sqrt (2.) * sqrt (7.)
               * (641. * tau2 * tau2 * tau2 + 1280. * tau1 * tau1 * tau1 + 7302. * tau2 * tau2 * tau1
                  + 7590. * tau2 * tau1 * tau1 + 222. * tau1 - 1038. * tau1 * tau1 - 9. + 183. * tau2
                  - 711. * tau2 * tau2 - 3252. * tau2 * tau1)
               * sqrt (187085.);
    case 29:
      if (e == 0)
        return (4. / 47457369.) * sqrt (7.)
               * (-7381. + 16722. * tau2 + 45879. * tau1 - 1326. * tau2 * tau2 - 8015. * tau2 * tau2 * tau2
                  - 93351. * tau1 * tau1 + 62234. * tau1 * tau1 * tau1 + 15087. * tau2 * tau2 * tau1
                  + 93351. * tau2 * tau1 * tau1 - 80916. * tau2 * tau1)
               * sqrt (13736271805.);
      else if (e == 1)
        return (4. / 47457369.) * sqrt (7.)
               * (-533. - 5642. * tau2 + 2531. * tau1 - 1326. * tau2 * tau2 + 3601. * tau2 * tau2 * tau2
                  - 3879. * tau1 * tau1 + 1914. * tau1 * tau1 * tau1 - 753. * tau2 * tau2 * tau1
                  - 11097. * tau2 * tau1 * tau1 + 16284. * tau2 * tau1)
               * sqrt (13736271805.);
      else if (e == 2)
        return -(52. / 47457369.) * sqrt (7.)
               * (-19. + 50. * tau2 + 325. * tau1 - 42. * tau2 * tau2 + 11. * tau2 * tau2 * tau2 - 861. * tau1 * tau1
                  + 574. * tau1 * tau1 * tau1 + 309. * tau2 * tau2 * tau1 + 861. * tau2 * tau1 * tau1
                  - 636. * tau2 * tau1)
               * sqrt (13736271805.);
      else
        return (4. / 47457369.) * sqrt (7.)
               * (-33. + 970. * tau2 + 515. * tau1 - 5694. * tau2 * tau2 + 8657. * tau2 * tau2 * tau2
                  - 1863. * tau1 * tau1 + 1914. * tau1 * tau1 * tau1 + 27183. * tau2 * tau2 * tau1
                  + 16839. * tau2 * tau1 * tau1 - 9636. * tau2 * tau1)
               * sqrt (13736271805.);
    default:;
    }
  default:;
  };
}

/* Here we define a function that calls the correct motherwavelet. */
double
muttermultiwavelet (int p, int i, double tau1, double tau2)
{
  if ((tau1 < 0.) || (tau1 > 1.) || (tau2 < 0.) || (tau2 > 1.) || (tau1 + tau2 > 1.))
    return 0.;
  int e;
  if ((tau1 <= 0.5) && (tau2 <= 0.5) && (tau1 + tau2 >= 0.5))
    e = 0;
  else if (tau1 > 0.5)
    e = 1;
  else if (tau2 > 0.5)
    e = 2;
  else
    e = 3;
  return muttermultiwavelets (p, i, tau1, tau2, e);
}

double
skalierungsfunktion_nextlevel (int i, double tau1, double tau2)
{
  if ((tau1 < 0.) || (tau1 > 1.) || (tau2 < 0.) || (tau2 > 1.) || (tau1 + tau2 > 1.))
    return 0.;
  return skalierungsfunktion (i, tau1, tau2);
}
//T8_EXTERN_C_END ();
