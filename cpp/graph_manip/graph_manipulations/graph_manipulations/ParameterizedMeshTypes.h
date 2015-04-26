#ifndef PARAMETERIZED_MESH_TYPES_H
#define PARAMETERIZED_MESH_TYPES_H
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
//#include <CGAL/Parameterization_mesh_patch_3.h>
#include "CGAL/Parameterization_mesh_patch_3.h"

template < class Polyhedron_ >
struct ParameterizedMeshTypes
{
    typedef Polyhedron_ PolyhedronType;
    typedef CGAL::Parameterization_polyhedron_adaptor_3< PolyhedronType > 
        AdaptorType;
    typedef CGAL::Parameterization_mesh_patch_3< AdaptorType >
        PatchType;
};

/**
Provides types related to the ParameterizedMesh.

The base template assumes we are using a PolyhedronPatch, ie.
a PolyhedronAdaptor wrapped in a MeshPatch ( see ParameterizedMeshTypes ).

Not sure if you see more layers of wrapping in practice, eg. when you need
to introduce multiple seams on more complicated topologies? TODO: figure out.
*/
template < class ParameterizedMesh >
struct ParameterizedMeshTraits
{
    typedef ParameterizedMesh PatchType;
    typedef typename PatchType::Adaptor AdaptorType;
    typedef typename AdaptorType::Polyhedron PolyhedronType;
    typedef typename PolyhedronType::Kernel Kernel;
};
#endif // PARAMETERIZED_MESH_TYPES_H
