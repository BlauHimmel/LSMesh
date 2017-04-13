#ifndef LS_MESH
#define LS_MESH

#include "GeoScene.h"
#include <vector>
#include <gl\glut.h>
#include <gl\GL.h>
#include <gl\GLU.h>
#include <gl\GLAUX.H>

class CLSMesh 
{
public:
	CLSMesh(int nVertices, int nEdges, int nTriangles,
		CGeoScene::CVertex* pVertices, CGeoScene::CEdge* pEdges,
		CGeoScene::CTriangle* pTriangles, std::vector<int>& controlPoints) :m_nVertices(nVertices),
		m_nEdges(nEdges), m_nTriangles(nTriangles)
	{
		m_pVertices = pVertices;
		m_pEdges = pEdges;
		m_pTriangles = pTriangles;
		m_ControlPoints = controlPoints;
	}
	CLSMesh(CLSMesh& other);
	~CLSMesh();

public:
	float* createMatrixA(float* &fpMatrix);
	float* createMatrixB(float* &fpMatrix);
	void calculate();

private:

	///////////////////////////////
	//    �����ȴ洢����ؾ������   //
	///////////////////////////////

	float* transpose(float* matrix, int i, int j);
	float getMatrix(float* matrix, int ld, int i, int j);
	void setMatrix(float* matrix, int ld, int i, int j, float data);
	void print(const char* filename, float* m, int ii, int jj);

private:

	/////////////////////
	//	    ԭʼ����    //
	/////////////////////

	int m_nVertices;
	CGeoScene::CVertex* m_pVertices;
	int m_nEdges;
	CGeoScene::CEdge* m_pEdges;
	int m_nTriangles;
	CGeoScene::CTriangle* m_pTriangles;

private:
	/////////////////////
	//	    ���Ƶ�      //
	/////////////////////

	std::vector<int> m_ControlPoints;

private:
	/////////////////////
	//	    �ؽ���      //
	/////////////////////
	std::vector<CPoint3D*> m_VerticesResult;
};

#endif