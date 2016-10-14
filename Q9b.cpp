#include<iostream>
#include<fstream>
#include<stdio.h>
#include<vector>
#include<math.h>
#include"tnt.h"
#include"jama_eig.h"

using namespace std;

void main()
{
	int n, m, count, t;
	float temp, sum;

	cout << "\nEnter the number of Vectors: ";
	cin >> n;

	cout << "\nEnter the length of each Vector: ";
	cin >> m;

	float temp1[6][4] = {
		{-1,-2,-1,0 },
		{ 2,1,3,2 },
		{ 1,2,0,3 },
		{ 2,3,1,1 },
		{ -1,2,3,1 },
		{ 0,1,-1,-2 }
	};

	Array2D<float> TVector(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			TVector[i][j] = temp1[i][j];
		}
	}
	/*
	for (int i = 0; i < n; i++)
	{
	cout << "\nEnter Vector " << i << " : ";
	for (int j = 0; j < m; j++)
	{
	cin >> temp;
	TVector[j][i] = temp;
	}
	}*/

	//Normalizing the matrix
	/*
	for (int i = 0; i < m; i++)
	{
	sum = 0;
	count = 0;
	for (int j = 0; j < n; j++)
	{
	sum += TVector[i][j];
	count++;
	}

	for (int j = 0; j < n; j++)
	{
	TVector[i][j] -= (sum/count);
	}
	}
	*/

	//cout << "\nAfter normalizing:\n";
	cout << "\nA Matrix: \n" << TVector;

	Array2D<float> ATranspose(n, m);

	//A Transpose
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ATranspose[j][i] = TVector[i][j];
		}
	}

	//Compute L
	Array2D<float> L(n, n, 0.0);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < m; k++)
			{
				L[i][j] += ATranspose[i][k] * TVector[k][j];
			}
		}
	}


	// Create EigenValue object
	JAMA::Eigenvalue<float> eigs(L);
	Array1D<float> Eigenvalues(n);
	Array2D<float> Eigenvectors(n, n);

	//Get real Eigenvalues
	eigs.getRealEigenvalues(Eigenvalues);
	cout << "\nEigen values: \n" << Eigenvalues;

	//Get Eigenvectors
	eigs.getV(Eigenvectors);
	cout << Eigenvectors;

	Array2D<float> SortedEigenv(n, n);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			SortedEigenv[i][n - 1 - j] = Eigenvectors[i][j];
		}
	}
	Eigenvectors = SortedEigenv;
	cout << "\nEigenVectors\n" << Eigenvectors;

	//Get eigenvectors of C
	Array2D<float> CVector(m, n);
	CVector = TNT::matmult(TVector, Eigenvectors);
	cout << "\nC Eigen Vectors: \n" << CVector;

	//Get unit Eigenvectors by normalizing each Eigenvector
	for (int i = 0; i < n; i++)
	{
		sum = 0;
		for (int j = 0; j < m; j++)
		{
			sum = sum + pow(CVector[j][i], 2);
		}

		for (int j = 0; j < m; j++)
		{
			CVector[j][i] /= pow(sum, 0.5);
		}
	}

	Array2D<float> ScoringMatrix(n, n);

	ScoringMatrix = TNT::matmult(ATranspose, CVector);

	Array2D<float>TScoring(n, n);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			TScoring[j][i] = ScoringMatrix[i][j];
		}
	}

	cout << "\n\nScoring Matrix:\n" << TScoring;

	Array2D<float>SignificantScoring(n - 1, n);

	for (int i = 0; i < n - 1; i++)
	{
		for (int j = 0; j < n; j++)
		{
			SignificantScoring[i][j] = TScoring[i][j];
		}
	}

	cout << "\nMost significant Scoring Matrix: \n";
	cout << SignificantScoring;

}
/*

Enter the number of Vectors: 4

Enter the length of each Vector: 6

A Matrix:
6 4
-1 -2 -1 0
2 1 3 2
1 2 0 3
2 3 1 1
-1 2 3 1
0 1 -1 -2

Eigen values:
4
4.30397
8.97479
10.6452
50.076

4 4
-0.779274 0.379268 -0.368967 0.335786
0.414693 -0.284748 -0.634452 0.586873
-0.301719 -0.559617 0.552523 0.538993
0.360178 0.679635 0.395035 0.50231

EigenVectors
4 4
0.335786 -0.368967 0.379268 -0.779274
0.586873 -0.634452 -0.284748 0.414693
0.538993 0.552523 -0.559617 -0.301719
0.50231 0.395035 0.679635 0.360178

C Eigen Vectors:
6 4
-2.04852 1.08535 0.749844 0.251606
3.88004 1.07525 0.154207 -1.32865
3.01646 -0.452765 1.84868 1.13065
3.47349 -1.69373 0.0243109 -0.256008
2.95725 1.15267 -1.94798 1.06368
-0.95674 -1.97704 -1.0844 -0.0039438


Scoring Matrix:
4 4
2.37617 4.15297 3.81415 3.55457
-1.20383 -2.07003 1.80272 1.28888
1.13621 -0.853046 -1.6765 2.03605
-1.61668 0.860323 -0.625945 0.747225

Most significant Scoring Matrix:
3 4
2.37617 4.15297 3.81415 3.55457
-1.20383 -2.07003 1.80272 1.28888
1.13621 -0.853046 -1.6765 2.03605
Press any key to continue . . .
*/