/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Söderström (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#include "SimpleDecimationMesh.h"

void SimpleDecimationMesh::computeCollapse(EdgeCollapse *collapse) {
  // The new vertex position is implicitly stored as the
  // position halfway along the edge. The cost is computed as
  // the vertex-to-vertex distance between the new vertex
  // and the old vertices at the edge's endpoints
  const size_t index0 = mEdges[collapse->halfEdge].vert;
  const size_t index1 = mEdges[mEdges[collapse->halfEdge].pair].vert;

  //Verts
  const Vector3<float> &v0 = mVerts[index0].pos;
  const Vector3<float> &v1 = mVerts[index1].pos;

 /* std::vector<size_t> verts1 = FindNeighborVertices(index0);
  std::vector<size_t> verts2 = FindNeighborVertices(index1);
  const Vector3<float> n = (VertexNormal(index0) + VertexNormal(index1)).Normalize();
  float sum = 100;
  for (size_t i = 0; i < verts1.size(); i++)
  {
      sum -= n * VertexNormal(verts1[i]);
  }
  for (size_t i = 0; i < verts2.size(); i++) {
      sum -= n * VertexNormal(verts2[i]);
  }*/

 // sum = sum < 0 ? 0 : sum;
  const Vector3<float> n = (VertexNormal(index0) + VertexNormal(index1)).Normalize();
  float sum = 1.0-VertexNormal(index0)*VertexNormal(index1);
  collapse->position = (v0 + v1) * 0.5;
  collapse->cost = sum;  //(collapse->position - v0).Length();
  //std::cout << "Orginal " << collapse->cost << " My " << sum << std::endl;
}
