/*************************************************************************************************
*
* Modeling and animation (TNM079) 2007
* Code base for lab assignments. Copyright:
*   Gunnar Johansson (gunnar.johansson@itn.liu.se)
*   Ken Museth (ken.museth@itn.liu.se)
*   Michael Bang Nielsen (bang@daimi.au.dk)
*   Ola Nilsson (ola.nilsson@itn.liu.se)
*   Andreas Sderstrm (andreas.soderstrom@itn.liu.se)
*
*************************************************************************************************/
#include "QuadricDecimationMesh.h"

const QuadricDecimationMesh::VisualizationMode
QuadricDecimationMesh::QuadricIsoSurfaces =
NewVisualizationMode("Quadric Iso Surfaces");

void QuadricDecimationMesh::Initialize() {
	// Allocate memory for the quadric array
	size_t numVerts = mVerts.size();
	mQuadrics.reserve(numVerts);
	std::streamsize width = std::cerr.precision(); // store stream precision
	for (size_t i = 0; i < numVerts; i++) {

		// Compute quadric for vertex i here
		mQuadrics.push_back(createQuadricForVert(i));

		// Calculate initial error, should be numerically close to 0

		Vector3<float> v0 = mVerts[i].pos;
		Vector4<float> v(v0[0], v0[1], v0[2], 1);
		Matrix4x4<float> m = mQuadrics.back();

		float error = v * (m * v);
		// std::cerr << std::scientific << std::setprecision(2) << error << " ";
	}
	std::cerr << std::setprecision(width) << std::fixed; // reset stream precision

														 // Run the initialize for the parent class to initialize the edge collapses
	DecimationMesh::Initialize();
}

/*! \lab2 Implement the computeCollapse here */
/*!
* \param[in,out] collapse The edge collapse object to (re-)compute,
* DecimationMesh::EdgeCollapse
*/
void QuadricDecimationMesh::computeCollapse(EdgeCollapse *collapse) {
	// Compute collapse->position and collapse->cost here
	// based on the quadrics at the edge endpoints

	//Get the vertex position
	Vector3<float> p1 = v(e(collapse->halfEdge).vert).pos;
	Vector3<float> p2 = v(e(e(collapse->halfEdge).next).vert).pos;
	Vector3<float> temp = (p1 + p2) / 2;

	//Get Qbar
	Matrix4x4<float> Qbar = mQuadrics[e(collapse->halfEdge).vert] + mQuadrics[e(e(collapse->halfEdge).next).vert];

	//Create Qhat
	Matrix4x4<float> Qhat = Qbar;
	Qhat(3, 0) = 0;
	Qhat(3, 1) = 0;
	Qhat(3, 2) = 0;
	Qhat(3, 3) = 1;

	//Get the zero-vector
	Vector4<float> Zero = { { 0 },{ 0 },{ 0 },{ 1 } };

	//Create the four different Vbars
	Vector4<float> v1 = Qhat.Inverse()*Zero;
	Vector4<float> v2(temp[0], temp[1], temp[2], 1.0f);
	Vector4<float> v3(p1[0], p1[1], p1[2], 1.0f);
	Vector4<float> v4(p2[0], p2[1], p2[2], 1.0f);
	//std::cerr << "computeCollapse in QuadricDecimationMesh not implemented.\n";

	//Calculate the different costs
	float c1 = Qbar(0, 0)*v1[0] * v1[0] + 2 * Qbar(0, 1)*v1[0] * v1[1] + 2 * Qbar(0, 2)*v1[0] * v1[2] + 2 * Qbar(0, 3)*v1[0] + Qbar(1, 1)*v1[1] * v1[1] + 2 * Qbar(1, 2)*v1[1] * v1[2] + 2 * Qbar(1, 3)*v1[1] + Qbar(2, 2)*v1[2] * v1[2] + 2 * Qbar(2, 3)*v1[2] + Qbar(3, 3);
	float c2 = Qbar(0, 0)*v2[0] * v2[0] + 2 * Qbar(0, 1)*v2[0] * v2[1] + 2 * Qbar(0, 2)*v2[0] * v2[2] + 2 * Qbar(0, 3)*v2[0] + Qbar(1, 1)*v2[1] * v2[1] + 2 * Qbar(1, 2)*v2[1] * v2[2] + 2 * Qbar(1, 3)*v2[1] + Qbar(2, 2)*v2[2] * v2[2] + 2 * Qbar(2, 3)*v2[2] + Qbar(3, 3);
	float c3 = Qbar(0, 0)*v3[0] * v3[0] + 2 * Qbar(0, 1)*v3[0] * v3[1] + 2 * Qbar(0, 2)*v3[0] * v3[2] + 2 * Qbar(0, 3)*v3[0] + Qbar(1, 1)*v3[1] * v3[1] + 2 * Qbar(1, 2)*v3[1] * v3[2] + 2 * Qbar(1, 3)*v3[1] + Qbar(2, 2)*v3[2] * v3[2] + 2 * Qbar(2, 3)*v3[2] + Qbar(3, 3);
	float c4 = Qbar(0, 0)*v4[0] * v4[0] + 2 * Qbar(0, 1)*v4[0] * v4[1] + 2 * Qbar(0, 2)*v4[0] * v4[2] + 2 * Qbar(0, 3)*v4[0] + Qbar(1, 1)*v4[1] * v4[1] + 2 * Qbar(1, 2)*v4[1] * v4[2] + 2 * Qbar(1, 3)*v4[1] + Qbar(2, 2)*v4[2] * v4[2] + 2 * Qbar(2, 3)*v4[2] + Qbar(3, 3);

	//Find the smallest cost	
	if (c1 < c2 && c1 < c3 && c1 < c4)
	{
		collapse->cost = c1;
		collapse->position[0] = v1[0];
		collapse->position[1] = v1[1];
		collapse->position[2] = v1[2];
	}
	else if (c2 < c3 && c2 < c4)
	{
		collapse->cost = c2;
		collapse->position = temp;
	}
	else if (c3 < c4)
	{
		collapse->cost = c3;
		collapse->position = p1;
	}
	else
	{
		collapse->cost = c4;
		collapse->position = p2;
	}


}

/*! After each edge collapse the vertex properties need to be updated */
void QuadricDecimationMesh::updateVertexProperties(size_t ind) {
	DecimationMesh::updateVertexProperties(ind);
	mQuadrics[ind] = createQuadricForVert(ind);
}

/*!
* \param[in] indx vertex index, points into HalfEdgeMesh::mVerts
*/
Matrix4x4<float>
QuadricDecimationMesh::createQuadricForVert(size_t indx) const {
	float q[4][4] = { { 0, 0, 0, 0 },{ 0, 0, 0, 0 },{ 0, 0, 0, 0 },{ 0, 0, 0, 0 } };
	Matrix4x4<float> Q(q);

	// The quadric for a vertex is the sum of all the quadrics for the adjacent
	// faces Tip: Matrix4x4 has an operator +=
	Vertex v = mVerts[indx];
	std::vector<size_t> faces = FindNeighborFaces(indx);

	for (size_t i = 0; i < faces.size(); i++)
	{
		Q += createQuadricForFace(faces[i]);
	}

	return Q;
}

/*!
* \param[in] indx face index, points into HalfEdgeMesh::mFaces
*/
Matrix4x4<float>
QuadricDecimationMesh::createQuadricForFace(size_t indx) const {

	// Calculate the quadric (outer product of plane parameters) for a face
	// here using the formula from Garland and Heckbert
	Vector3<float> n = mFaces[indx].normal;

	float a = n[0];
	float b = n[1];
	float c = n[2];

	Vector3<float> pos = mVerts[mEdges[mFaces[indx].edge].vert].pos;

	float d = -1 * (a * pos[0] + b * pos[1] + c * pos[2]);

	float Kpi[4][4] = { { a*a, a*b, a*c, a*d } ,{ a*b, b*b, b*c, b*d } ,{ a*c, b*c, c*c, c*d } ,{ a*d, b*d, c*d, d*d } };

	return Matrix4x4<float>(Kpi);
}

void QuadricDecimationMesh::Render() {
	DecimationMesh::Render();

	glEnable(GL_LIGHTING);
	glMatrixMode(GL_MODELVIEW);

	if (mVisualizationMode == QuadricIsoSurfaces) {
		// Apply transform
		glPushMatrix(); // Push modelview matrix onto stack

						// Implement the quadric visualization here
		std::cout << "Quadric visualization not implemented" << std::endl;

		// Restore modelview matrix
		glPopMatrix();
	}
}
