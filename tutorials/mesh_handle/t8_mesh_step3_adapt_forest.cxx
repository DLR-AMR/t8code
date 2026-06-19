/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2015 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file t8_mesh_step3_adapt_forest.cxx
 * This is the same as general/t8_step3_adapt_forest.cxx but using the mesh handle interface instead of the forest
 * interface.
*/

#include <t8.h>
#include <mesh_handle/mesh.hxx>
#include <mesh_handle/competence_pack.hxx>
#include <mesh_handle/constructor_wrappers.hxx>
#include <mesh_handle/mesh_io.hxx>
#include <mesh_handle/concepts.hxx>
#include<t8_types/t8_vec.hxx>
#include <memory>
#include <span>

struct adapt_data {
    std::array<double, 3> midpoint;
    double refine_radius;
    double coarsen_radius;
};

template <t8_mesh_handle::T8MeshType TMeshClass>
int adapt_callback(const TMeshClass &mesh, std::span<const typename TMeshClass::element_class> elements, const adapt_data &adapt_data) {
    auto element_centroid = elements[0].get_centroid();
    double dist = t8_dist<t8_3D_vec, t8_3D_vec>(element_centroid, adapt_data.midpoint);
    if (dist < adapt_data.refine_radius) {
        return 1; // refine
    } else if ((elements.size() > 1) && (dist > adapt_data.coarsen_radius)) {
        return -1; // coarsen
    }
    return 0; // do nothing
}

template <t8_mesh_handle::T8MeshType TMeshClass>
std::unique_ptr<TMeshClass> build_mesh(sc_MPI_Comm comm, int level) {
    auto mesh = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<TMeshClass>(level, comm);
    struct adapt_data adapt_params = {
        {0.5, 0.5, 1.0},
        0.2,
        0.4
    };
    mesh->set_balance();
    mesh->set_partition();
    mesh->set_adapt(
        TMeshClass::template mesh_adapt_callback_wrapper<adapt_data>(
            adapt_callback<TMeshClass>, 
            adapt_params
        )
    );
    mesh->set_ghost();
    return mesh;
}

int main(int argc, char** argv) {
    int mpiret = sc_MPI_Init(&argc, &argv);
    SC_CHECK_MPI(mpiret);

    sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);

    t8_init (SC_LP_PRODUCTION);

    sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

    using mesh_type = t8_mesh_handle::mesh<>;

    int uniform_level = 2;
    auto mesh = build_mesh<mesh_type>(comm, uniform_level);

    t8_mesh_handle::write_mesh_to_vtk(mesh, "adapted_mesh.vtu");

    sc_finalize();
    return 0;
}
