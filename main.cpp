#include <Eigen/Dense>
#include <cinolib/meshes/meshes.h>
#include <cinolib/vector_field.h>
#include <cinolib/scalar_field.h>
#include <cinolib/gradient.h>
#include <cinolib/heat_flow.h>

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// ALGORITHM PARAMETERS
#define LAMBDA           0.1
#define SMOOTHING_PASSES 5

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

using namespace cinolib;

template <class Mesh>
ScalarField smooth_discrete_hyper_surface(Mesh & m)
{
    // STEP ONE: compute heat flow
    std::vector<uint> heat_sources;
    for(uint vid=0; vid<m.num_verts(); ++vid)
    {
        bool has_A = false;
        bool has_B = false;
        for(uint pid : m.adj_v2p(vid))
        {
            if (m.poly_data(pid).label == 0) has_A = true; else
            if (m.poly_data(pid).label == 1) has_B = true;
        }
        if (has_A && has_B) heat_sources.push_back(vid);
    }
    double t = std::pow(m.edge_avg_length(),2.0);
    ScalarField u = heat_flow(m, heat_sources, t, COTANGENT);
    u.normalize_in_01();
    u.copy_to_mesh(m);
    u.serialize("u.txt");
    std::cout << "heat flow computed (see file u.txt)" << std::endl;

    VectorField field = VectorField(m.num_polys());
    field = gradient_matrix(m) * u;
    field.serialize("u_gradient.txt");
    std::cout << "u gradient computed (see file u_gradient.txt)" << std::endl;

    // STEP TWO: flip the gradient of one of the regions
    for(uint pid=0; pid<m.num_polys(); ++pid)
    {
        if (m.poly_data(pid).label == 1) field.set(pid, -field.vec_at(pid));
    }
    field.normalize();
    field.serialize("X.txt");
    std::cout << "X field generated (see file X.txt)" << std::endl;

    // STEP THREE: smooth the resulting gradient
    for(uint i=0; i<SMOOTHING_PASSES; ++i)
    for(uint pid=0; pid<m.num_polys(); ++pid)
    {
        vec3d avg_g = field.vec_at(pid);
        for(uint nbr : m.adj_p2p(pid))
        {
            avg_g += field.vec_at(nbr);
        }
        avg_g /= static_cast<double>(m.adj_p2p(pid).size()+1);
        avg_g.normalize();
        field.set(pid,avg_g);
    }
    field.normalize();
    field.serialize("X_prime.txt");
    std::cout << "smoothed X field generated (see file X_prime.txt)" << std::endl;

    // STEP FOUR: find the scalar field corresponding to it
    int type = UNIFORM;
    if (m.mesh_type()==TRIMESH || m.mesh_type()==TETMESH) type = COTANGENT; // use cotangent weights for tris and tets
    std::vector<Eigen::Triplet<double>> entries = laplacian_matrix_entries(m, type);
    Eigen::VectorXd div = gradient_matrix(m).transpose() * field;
    Eigen::SparseMatrix<double> L(m.num_verts()+heat_sources.size(), m.num_verts());
    Eigen::VectorXd rhs(m.num_verts()+heat_sources.size());
    for(uint vid=0; vid<m.num_verts(); ++vid) rhs[vid] = div[vid];
    for(uint i=0; i<heat_sources.size(); ++i)
    {
        uint vid = heat_sources.at(i);
        entries.push_back(Entry(m.num_verts()+i, vid, LAMBDA));
        rhs[m.num_verts()+i] = 0.0;
    }
    L.setFromTriplets(entries.begin(), entries.end());
    ScalarField phi;
    solve_least_squares(-L, rhs, phi);
    phi.copy_to_mesh(m);
    phi.normalize_in_01();
    return phi;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        std::cout << "\nThis is a demo implementation of the hyper surface relaxation described in:" << std::endl;
        std::cout << "A Heat Flow Based Relaxation Scheme for n Dimensional Discrete Hyper Surfaces" << std::endl;
        std::cout << "Marco Livesu                                                                 " << std::endl;
        std::cout << "Computers and Graphics, 2018                                                 " << std::endl;
        std::cout << "\nusage:\n\thyper_surface_smoothing mesh labeling\n                          " << std::endl;
        std::cout << "\tmesh     : a triangle mesh (both OBJ and OFF format are supported)         " << std::endl;
        std::cout << "\tlabeling : a bipartition of the mesh in the form of a text file            " << std::endl;
        std::cout << "\t           having one line per triangle, valued 0 or 1.\n                  " << std::endl;
        return -1;
    }

    // load both mesh and labeling and copy labeling onto mesh
    Trimesh<> m(argv[1]);
    ScalarField sf(argv[2]);
    assert(sf.size() == m.num_polys());
    for(uint pid=0; pid<m.num_polys(); ++pid) m.poly_data(pid).label = sf[pid];

    // Hyper Surface Relaxation.
    // Notice that this function is templated, and will work with any mesh
    // supported by the CinoLib. All you have to do is to change the parsing
    // of the input argument argv[1], substituting Trimesh<> with Quadmesh<>,
    // or Polygonmesh<> or Tetmesh<> or Hexmesh<> or Polyhedralmesh<>.
    //
    ScalarField res = smooth_discrete_hyper_surface(m);
    res.serialize("res.txt");
    std::cout << "Output scalar field computed (see file res.txt). The smoothed" << std::endl;
    std::cout << "boundary corresponds to the zero level set of such field.\n  " << std::endl;

    return 0;
}
