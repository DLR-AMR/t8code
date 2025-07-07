/**
 * @file
 * @brief Dunavant quadrature rule over the interior of a triangle in 2D 
 * See https://people.sc.fsu.edu/~jburkardt/m_src/triangle_dunavant_rule/triangle_dunavant_rule.html
 */

#pragma once

#ifdef T8_ENABLE_MRA

namespace t8_mra
{
int
dunavant_degree (int rule);
int
dunavant_order_num (int rule);
void
dunavant_rule (int rule, int order_num, double xy[], double w[]);
int
dunavant_rule_num (void);
int *
dunavant_suborder (int rule, int suborder_num);
int
dunavant_suborder_num (int rule);
void
dunavant_subrule (int rule, int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_01 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_02 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_03 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_04 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_05 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_06 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_07 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_08 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_09 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_10 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_11 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_12 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_13 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_14 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_15 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_16 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_17 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_18 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_19 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
dunavant_subrule_20 (int suborder_num, double suborder_xyz[], double suborder_w[]);
void
file_name_inc (char *file_name);
int
i4_max (int i1, int i2);
int
i4_min (int i1, int i2);
int
i4_modp (int i, int j);
int
i4_wrap (int ival, int ilo, int ihi);
double
r8_huge (void);
int
r8_nint (double x);
void
reference_to_physical_t3 (double t[], int n, double ref[], double phy[]);
int
s_len_trim (char *s);
void
timestamp (void);
char *
timestring (void);
double
triangle_area (double t[2 * 3]);
void
triangle_points_plot (char *file_name, double node_xy[], int node_show, int point_num, double point_xy[],
                      int point_show);
}  // namespace t8_mra

#endif
