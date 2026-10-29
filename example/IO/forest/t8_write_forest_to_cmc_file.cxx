#include <t8.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_partition.h>

#include "cmc.hxx"
#include "cmc_amr_mesh_compression.hxx"

#include <string>

namespace t8::example::cmc::mesh_compressor
{

constexpr bool kWriteVTKDataTestAMRMeshCompression = true;

/** An exemplary adaptation function to adapt the uniform forest.
 */
static t8_locidx_t
TestAdapt ([[maybe_unused]] t8_forest_t forest,
           t8_forest_t forest_from,
           t8_locidx_t which_tree,
           [[maybe_unused]] const t8_eclass_t tree_class,
           [[maybe_unused]] t8_locidx_t lelement_id,
           [[maybe_unused]] const t8_scheme_c * ts,
           const int is_family,
           [[maybe_unused]] const int num_elements,
           [[maybe_unused]] t8_element_t * elements[])
{

    std::vector<double> elem_midpoint(3);
    t8_forest_element_centroid (forest_from, which_tree, elements[0], elem_midpoint.data());

    const t8_gloidx_t gtree_id = t8_forest_global_tree_id (forest_from, which_tree);
    if ((gtree_id == 3 || gtree_id == 7 || gtree_id == 10) && is_family) {return -1;}
    if (gtree_id == 4 || gtree_id == 5) {return 1;}
    if (std::fabs(elem_midpoint[0] - 0.5) <= 0.2 && std::fabs(elem_midpoint[1] - 0.5) <= 0.2)
    {
        return 1;
    } else
    {
        return 0;
    }
}

/** Create an example mesh for the compression example.
 * \return The example mesh
 */
t8_forest_t
CreateInitMesh()
{
    /* Create a mesh */
    const sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
    t8_cmesh_t cmesh;
    t8_cmesh_init(&cmesh);
    t8_cmesh_new_hypercube_hybrid (cmesh, comm, 0);

    const t8_scheme *scheme = t8_scheme_new_default ();

    const int initial_level = 4;
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, initial_level, 0, comm);
    
    /* Partition For coarsening */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);

    /* Partition the forest for coarsening before the adaptation */
    constexpr int partition_for_coarsening = 1;
    t8_forest_set_partition(partitioned_forest, forest, partition_for_coarsening);
    t8_forest_commit(partitioned_forest);

    partitioned_forest = t8_forest_new_adapt(partitioned_forest, TestAdapt, 0, 0, NULL);

    t8_forest_t adapted_partitioned_forest;
    t8_forest_init(&adapted_partitioned_forest);

    /* Partition the forest (we can supply a not partitioned for coarsening forest to the compression) */
    constexpr int partition_for_coarsening2 = 0;
    t8_forest_set_partition(adapted_partitioned_forest, partitioned_forest, partition_for_coarsening2);
    t8_forest_commit(adapted_partitioned_forest);

    return adapted_partitioned_forest;
}

/** Supply the same cmesh and scheme as it has been used in CreateInitMesh().
 * \return The cmesh and the scheme from the example mesh as a pair.
 */
std::pair<t8_cmesh_t, const t8_scheme *>
GetInitCMeshAndScheme()
{
    t8_cmesh_t cmesh;
    t8_cmesh_init(&cmesh);
    t8_cmesh_new_hypercube_hybrid (cmesh, sc_MPI_COMM_WORLD, 0);
    return std::make_pair(cmesh, t8_scheme_new_default());
}

}

int
main(void)
{
    /**** Initialize all components ****/
    /* Initialize MPI. This has to happen before we initialize sc or t8code. */
    const int rv_mpi_init = sc_MPI_Init (NULL, NULL);
    /* Error check the MPI return value. */
    SC_CHECK_MPI (rv_mpi_init);
    /* Initialize the sc library, has to happen before we initialize t8code. */
    sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
    /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
    t8_init (SC_LP_PRODUCTION);
    /* Initialize cmc in it's minimal form (**optional**, if t8code is already initialized) */
    cmc::CmcInitialize(cmc::kMinimumInitialization);

    {
        /* Define a file name for the compressed data */
        const std::string file_name("t8_cmc_example_amr_mesh_compressor.cmc");
        
        /**** SETUP of an example mesh ****/
        /* Create an initial mesh */
        t8_forest_t mesh = t8::example::cmc::mesh_compressor::CreateInitMesh();
        if constexpr (t8::example::cmc::mesh_compressor::kWriteVTKDataTestAMRMeshCompression)
        {
            t8_forest_write_vtk_ext (mesh, "t8_cmc_example_initial_mesh", 1, 1, 1, 1, 0, 0, 0, 0, NULL);
        }
        /* Store the global number of elements */
        const uint64_t num_global_init_elems = static_cast<uint64_t>(t8_forest_get_global_num_leaf_elements(mesh));
        cmc::cmc_debug_msg("The initial mesh has ", num_global_init_elems, " elements.");
        /* Store the initial data offset */
        const uint64_t init_mesh_offset = static_cast<uint64_t>(t8_forest_get_first_local_leaf_element_id(mesh));
        /* Declare a forest pointer for the reconstructed mesh */
        t8_forest_t reconstructed_mesh;

        /**** COMPRESSION ****/
        {
        /** Construct a mesh compressor. The only parameter it takes is the mesh whose adaptivity structure
         *  we want to compress/store. The compression will be performed via the MPI-Communicator of the mesh.
         **/
        /* Define the mesh compressor */
        cmc::mesh::Compressor mesh_compressor(mesh);

        /** The following function performs the actual compression/serialization of the adaptivity structure.
         *  We are able to supply a parameter which controls whether we want to perform the compression 
         *  strictly level-wise (parameter value: true) or not (parameter value: false).
         *  The strictly leve-wise approach ensures that only elements of a certain refinement level are 
         *  compressed corresponding to the compression step. In particular, this means that during the 
         *  decompression, all elements of refinement level N+1 are constructed by the N-th decompression
         *  step. This allows a controllable progressive level-wise reconstruction of the mesh. This comes
         *  at a cost of larger output files.
         *  The default (parameter value: false) does not compress strictly-levelwise and coarsens the 
         *  mesh more aggressively, yielding smaller compressed files and faster runtimes at a cost
         *  of a less controllable progressive reconstruction.
         **/
        /* Compress */
        constexpr bool compress_strictly_level_wise = false;
        mesh_compressor.CompressMeshStructure(compress_strictly_level_wise);

        /** After the mesh structure has been compressed, we are able to write out the compressed file
         *  at the supplied file_name. Currently, the data is collected on single rank, to reduce
         *  fragmented small write accesses. As a default, all processes wait until the file has been
         *  written to disk, however, if permitted, some ranks may return earlier although the file is
         *  not fully written and closed.
         **/
        /* Write the encoded data */
        constexpr bool barrier_until_the_file_is_written = true;
        mesh_compressor.WriteCompressedMeshStructure(file_name, barrier_until_the_file_is_written);
        }
        /**** END OF COMPRESSION ****/

        /** The following function allows to print out some meta information from the
         *  compressed mesh file, e.g. number of global elements per compression level, 
         *  the initial MPI-communicator size, etc.
         **/
        /* Print out some information about the compressed mesh */
        cmc::mesh::PrintMeshCompressionInfo(file_name);

        /**** DECOMPRESSION ****/
        {
        /** In order to reconstruct the mesh correctly, we need to supply the intial cmesh and the scheme,
         *  since the compressed .cmc-file only holds the adaptivity structure of the supplied mesh **/
        /* Get the initial cmesh and scheme */
        auto [cmesh, scheme] = t8::example::cmc::mesh_compressor::GetInitCMeshAndScheme();

        /** The mesh decompressor takes the file with the encoded adaptivity structure, the cmesh and the scheme
         *  as well as the MPI-Communicator over which the forest will be reconstructed.
         *  The MPI-Communicator can differ from the one of the intial mesh, in particular, it is permitted
         *  to decompress/reconstruct the mesh with a different amount (more, less or equal) of MPI-processes
         *  than it has been compressed/serialized with.
         */
        /* Define the mesh decompressor */
        cmc::mesh::Decompressor mesh_decompressor(file_name, MPI_COMM_WORLD, cmesh, scheme);

        /** This function call performs the progressive reconstruction of the mesh from the root level
         *  up to the leaf level of the initial forest.
         *  Depending on how the compression has been performed, the decompression is also either strictly
         *  level-wise or not.
         */
        /* Decompress the mesh */
        reconstructed_mesh = mesh_decompressor.ReconstructMesh();

        const uint64_t num_global_reconstructed_elements = static_cast<uint64_t>(t8_forest_get_global_num_leaf_elements(reconstructed_mesh));
        cmc::cmc_debug_msg("The reconstructed mesh has ", num_global_reconstructed_elements, " elements.");
        }
        /**** END OF DECOMPRESSION ****/

        if constexpr (t8::example::cmc::mesh_compressor::kWriteVTKDataTestAMRMeshCompression)
        {
            t8_forest_write_vtk_ext (reconstructed_mesh, "t8_cmc_example_reconstructed_mesh_", 1, 1, 1, 1, 0, 0, 0, 0, NULL);
        }

        /* Destroy all forests */
        t8_forest_unref(&mesh);
        t8_forest_unref(&reconstructed_mesh);
    }

    /**** Finalize all components ****/
    cmc::CmcFinalize();
    sc_finalize ();
    const int rv_mpi_fin = sc_MPI_Finalize ();
    SC_CHECK_MPI (rv_mpi_fin);

    return 0;
}