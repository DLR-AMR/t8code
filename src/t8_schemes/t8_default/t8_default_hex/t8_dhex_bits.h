#ifndef T8_DHEX_BITS_H
#define T8_DHEX_BITS_H

#include <t8.h>
#include <t8_schemes/t8_default/t8_default_hex/t8_dhex.h>

T8_EXTERN_C_BEGIN ();

void t8_hex_parent (const t8_phex_t * q, t8_phex_t * r);

void t8_hex_sibling (const t8_phex_t * q, t8_phex_t * r, int sibling_id);

int t8_hex_compare (const void *v1, const void *v2);

int t8_hex_is_extended (const t8_phex_t * q);

int t8_hex_is_parent (const t8_phex_t * q, const t8_phex_t * r);

void t8_hex_childrenpv (const t8_phex_t * q, t8_phex_t * c[]);

int t8_hex_child_id (const t8_phex_t * q);

int t8_hex_ancestor_id (const t8_phex_t * q, int level);

int t8_hex_is_familypv (t8_phex_t * q[]);

void t8_hex_set_morton (t8_phex_t * quadrant, int level, uint64_t id);

uint64_t t8_hex_linear_id (const t8_phex_t * quadrant, int level);

void t8_hex_first_descendant (const t8_phex_t * q, t8_phex_t * fd, int level);

void t8_hex_nearest_common_ancestor (const t8_phex_t * q1, const t8_phex_t * q2, t8_phex_t * r);

void t8_hex_corner_descendant (const t8_phex_t * q, t8_phex_t * r, int c, int level);

void t8_hex_face_neighbor (const t8_phex_t * q, int face, t8_phex_t * r);

int t8_hex_is_inside_root (const t8_phex_t * q);

void t8_hex_print (int log_priority, const t8_phex_t * q);

void t8_hex_last_descendant (const t8_phex_t * q, t8_phex_t * ld, int level);

int t8_hex_is_node (const t8_phex_t * q, int inside);

int t8_hex_is_inside_3x3 (const t8_phex_t * q);

int t8_hex_is_valid (const t8_phex_t * q);

int t8_hex_is_family (    const t8_phex_t * q0,
                          const t8_phex_t * q1,
                          const t8_phex_t * q2,
                          const t8_phex_t * q3,
                          const t8_phex_t * q4,
                          const t8_phex_t * q5,
                          const t8_phex_t * q6,
                          const t8_phex_t * q7);

void
t8_hex_children (const t8_phex_t * q,
                         t8_phex_t * c0, t8_phex_t * c1,
                         t8_phex_t * c2, t8_phex_t * c3,
                         t8_phex_t * c4, t8_phex_t * c5,
                         t8_phex_t * c6, t8_phex_t * c7
);

T8_EXTERN_C_END ();

#endif
