#ifndef T8_DQUAD_BITS_H
#define T8_DQUAD_BITS_H

#include <t8.h>
#include <t8_schemes/t8_default/t8_default_quad/t8_dquad.h>

T8_EXTERN_C_BEGIN ();

void t8_dquad_parent (const t8_dquad_t * q, t8_dquad_t * r);

void t8_dquad_sibling (const t8_dquad_t * q, t8_dquad_t * r, int sibling_id);

int t8_dquad_compare (const void *v1, const void *v2);

int t8_dquad_is_extended (const t8_dquad_t * q);

int t8_dquad_is_parent (const t8_dquad_t * q, const t8_dquad_t * r);

void t8_dquad_childrenpv (const t8_dquad_t * q, t8_dquad_t * c[]);

int t8_dquad_child_id (const t8_dquad_t * q);

int t8_dquad_ancestor_id (const t8_dquad_t * q, int level);

int t8_dquad_is_familypv (t8_dquad_t * q[]);

void t8_dquad_set_morton (t8_dquad_t * quadrant, int level, uint64_t id);

uint64_t t8_dquad_linear_id (const t8_dquad_t * quadrant, int level);

void t8_dquad_first_descendant (const t8_dquad_t * q, t8_dquad_t * fd, int level);

void t8_dquad_nearest_common_ancestor (const t8_dquad_t * q1, const t8_dquad_t * q2, t8_dquad_t * r);

void t8_dquad_corner_descendant (const t8_dquad_t * q, t8_dquad_t * r, int c, int level);

void t8_dquad_face_neighbor (const t8_dquad_t * q, int face, t8_dquad_t * r);

int t8_dquad_is_inside_root (const t8_dquad_t * q);

void t8_dquad_print (int log_priority, const t8_dquad_t * q);

void t8_dquad_last_descendant (const t8_dquad_t * q, t8_dquad_t * ld, int level);

int t8_dquad_is_node (const t8_dquad_t * q, int inside);

int t8_dquad_is_inside_3x3 (const t8_dquad_t * q);

int t8_dquad_is_valid (const t8_dquad_t * q);

int t8_dquad_is_family (const t8_dquad_t * q0,
                          const t8_dquad_t * q1,
                          const t8_dquad_t * q2,
                          const t8_dquad_t * q3);

void
t8_dquad_children (const t8_dquad_t * q,
                         t8_dquad_t * c0, t8_dquad_t * c1,
                         t8_dquad_t * c2, t8_dquad_t * c3);


T8_EXTERN_C_END ();

#endif /* T8_DEFAULT_QUAD_BITS_H */
