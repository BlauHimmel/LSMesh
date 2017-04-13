#include "stdafx.h"
#include "LSMesh.h"
#include <mkl.h>
#include <fstream>



CLSMesh::CLSMesh(CLSMesh& other)
{
	
}

CLSMesh::~CLSMesh()
{
	for (int i = 0; i < m_VerticesResult.size(); i++)
	{
		delete m_VerticesResult[i];
	}
}

float* CLSMesh::createMatrixA(float* &fpMatrix)
{
	memset(fpMatrix, 0, sizeof(float) * m_nVertices * (m_nVertices + m_ControlPoints.size()));
	for (int i = 0; i < m_nVertices; i++)
	{
		setMatrix(fpMatrix, m_nVertices, i, i, 1.0f);
	}
	for (int i = 0; i < m_nVertices; i++)
	{
		int nNeighbor = m_pVertices[i].m_nNeighbor;
		for (int j = 0; j < nNeighbor; j++)
		{
			setMatrix(fpMatrix, m_nVertices, i, m_pVertices[i].m_npNeighborVertexIndices[j], -1.0f / nNeighbor);
		}
	}
	for (int i = 0; i < m_ControlPoints.size(); i++)
	{
		setMatrix(fpMatrix, m_nVertices, m_nVertices + i, m_ControlPoints[i], 1.0f);
	}
	return fpMatrix;
}

float* CLSMesh::createMatrixB(float* &fpMatrix)
{
	memset(fpMatrix, 0, sizeof(float) * (m_nVertices + m_ControlPoints.size()) * 3);
	for (int i = 0; i < m_ControlPoints.size(); i++)
	{
		setMatrix(fpMatrix, 3, m_nVertices + i, 0, m_pVertices[m_ControlPoints[i]].m_vPosition.x);
		setMatrix(fpMatrix, 3, m_nVertices + i, 1, m_pVertices[m_ControlPoints[i]].m_vPosition.y);
		setMatrix(fpMatrix, 3, m_nVertices + i, 2, m_pVertices[m_ControlPoints[i]].m_vPosition.z);
	}
	return fpMatrix;
}

void CLSMesh::print(const char* filename, float* m, int ii, int jj)
{
	std::ofstream fout(filename, std::ios::out);
	for (int i = 0; i < ii; i++)
	{
		for (int j = 0; j < jj; j++)
		{
			fout << getMatrix(m, jj, i, j) << "\t\t";
		}
		fout << std::endl;
	}
}

void CLSMesh::calculate()
{
	float* fpMatrixA = new float[(m_nVertices + m_ControlPoints.size()) * m_nVertices];		//(n+m)*n
	float* fpMatrixB = new float[(m_nVertices + m_ControlPoints.size()) * 3];			//(n+m)*3
	float* fpMatrixATA = new float[m_nVertices * m_nVertices];	//n*n
	float* fpMatrixATB = new float[m_nVertices * 3];			//n*3
	memset(fpMatrixATA, 0, sizeof(float) * m_nVertices * m_nVertices);
	memset(fpMatrixATB, 0, sizeof(float) * m_nVertices * 3);

	createMatrixA(fpMatrixA);
	createMatrixB(fpMatrixB);

	float* fpMatrixAT = transpose(fpMatrixA, m_nVertices + m_ControlPoints.size(), m_nVertices);
	
	//AT:n*(n+m) A:(n+m)*n C:n*n AT*A
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m_nVertices, m_nVertices, m_nVertices + m_ControlPoints.size(), 1.0f, fpMatrixAT, m_nVertices + m_ControlPoints.size(),
		fpMatrixA, m_nVertices, 0.0f, fpMatrixATA, m_nVertices);
	//AT:n*(n+m) B:(n+m)*3 C:n*3 AT*B
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m_nVertices, 3, m_nVertices + m_ControlPoints.size(), 1.0f, fpMatrixAT, m_nVertices + m_ControlPoints.size(),
		fpMatrixB, 3, 0.0f, fpMatrixATB, 3);
	int* ipiv = new int[m_nVertices * m_nVertices];
	int result = LAPACKE_sgesv(LAPACK_ROW_MAJOR, m_nVertices, 3, fpMatrixATA, m_nVertices, ipiv, fpMatrixATB, 3);

	if (result == 0)
	{
		for (int i = 0; i < m_nVertices; i++)
		{
			//CPoint3D* point = new CPoint3D(getMatrix(fpMatrixATB, 3, i, 0), getMatrix(fpMatrixATB, 3, i, 1), getMatrix(fpMatrixATB, 3, i, 2));
			//m_VerticesResult.push_back(point);
			m_pVertices[i].m_vPosition.x = getMatrix(fpMatrixATB, 3, i, 0);
			m_pVertices[i].m_vPosition.y = getMatrix(fpMatrixATB, 3, i, 1);
			m_pVertices[i].m_vPosition.z = getMatrix(fpMatrixATB, 3, i, 2);
		}
	}
	else
	{
		printf("------------------------------------------------------------------");
		printf("error occur during the process of solving the linear equation:");
		printf("error code:%d", result);
		printf("If info = -i, parameter i had an illegal value.");
		printf("If info = i,U(i,i)(computed in double precision for mixed precision");
		printf("subroutines) is exactly zero. The factorization has been completed,");
		printf("but the factor U is exactly singular, so the solution could not be ");
		printf("computed.");
		printf("------------------------------------------------------------------");
		delete[] ipiv;
		delete[] fpMatrixA;
		delete[] fpMatrixB;
		delete[] fpMatrixATA;
		delete[] fpMatrixATB;
		assert("ERROR:CANNOT SOLVE EQUATION!");
	}
	delete[] ipiv;
	delete[] fpMatrixA;
	delete[] fpMatrixB;
	delete[] fpMatrixATA;
	delete[] fpMatrixATB;
	delete[] fpMatrixAT;
}

float CLSMesh::getMatrix(float* matrix, int ld, int i, int j)
{
	assert(ld >= 1);
	assert(i >= 0 && j >= 0);
	return matrix[i * ld + j];
}

void CLSMesh::setMatrix(float* matrix, int ld, int i, int j, float data)
{
	assert(ld >= 1);
	assert(i >= 0 && j >= 0);
	matrix[i * ld + j] = data;
}

float* CLSMesh::transpose(float* matrix, int i, int j)
{
	float* trans = new float[i * j];
	memset(trans, 0, sizeof(float) * i * j);
	for (int _i = 0; _i < i; _i++)
	{
		for (int _j = 0; _j < j; _j++)
		{
			float tmp = getMatrix(matrix, j, _i, _j);
			setMatrix(trans, i, _j, _i, tmp);
		}
	}
	return trans;
}
