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
	int n, m, count,t;
	float temp, sum;

	cout << "\nEnter the number of Vectors: ";
	cin >> n;

	cout << "\nEnter the length of each Vector: ";
	cin >> m;

	float temp1[8][4] = {
		{2,-2,-1,3},
		{-1,3,3,-1},
		{0,2,3,0},
		{1,3,1,3},
		{1,0,-1,2},
		{-3,2,4,-1},
		{5,-1,5,3},
		{2,1,2,0}
	};

	Array2D<float> TVector(m,n);
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
	Array1D<float> Eigenvalues (n);
	Array2D<float> Eigenvectors(n,n);

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
			SortedEigenv[i][n-1-j] = Eigenvectors[i][j];
		}
	}
	Eigenvectors = SortedEigenv;
	cout << "\nEigenVectors\n" << Eigenvectors;
	
	//Get eigenvectors of C
	Array2D<float> CVector(m,n);	
	CVector = TNT::matmult(TVector, Eigenvectors);
	cout << "\nC Eigen Vectors: \n" << CVector;

	//Get unit Eigenvectors by normalizing each Eigenvector
	for (int i = 0; i < n; i++)
	{
		sum = 0;
		for (int j = 0; j < m; j++)
		{
			sum = sum + pow(CVector[j][i],2);
		}

		for (int j = 0; j < m; j++)
		{
			CVector[j][i] /= pow(sum,0.5);
		}
	}
	
	Array2D<float> ScoringMatrix(n,n);

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

	Array2D<float>SignificantScoring(n-1, n);

	for (int i = 0; i < n-1; i++)
	{
		for (int j = 0; j < n; j++)
		{
			SignificantScoring[i][j] = TScoring[i][j];
		}
	}

	cout << "\nMost significant Scoring Matrix: \n";
	cout << SignificantScoring;

	float temp2[4][8] = {
		{1,5,1,5,5,1,1,3},
		{-2,3,2,3,0,2,-1,1},
		{2,-3,2,3,0,0,2,-1},
		{2,-2,2,2,-1,1,2,2}
	};

	Array2D<float> TestVector(4,8);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			TestVector[i][j] = temp2[i][j];
		}
	}

	Array2D<float>EigenSpace(n,n);

	EigenSpace = TNT::matmult(TestVector, CVector);

	//cout << "\nEigenspace: \n";
	//cout << EigenSpace;

	Array2D<float>TEigenSpace(n,n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n-1; j++)
		{
			TEigenSpace[j][i] = EigenSpace[i][j];
		}
	}

	
	Array2D<float>Scores(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int k = 0; k < n; k++)
		{
			sum = 0;
			for (int j = 0; j < n - 1; j++)
			{
				sum += pow(TEigenSpace[j][i] - SignificantScoring[j][k], 2);
			}
			Scores[i][k] = pow(sum, 0.5);
		}
	}

	//cout << "Scoring matrix: \n";
	//cout << Scores;

	for (int i = 0; i < n; i++)
	{
		float min = Scores[i][0];

		cout << "\nEuclidean Distances for X-" << i + 1 << endl;
		for (int j = 0; j < n; j++)
		{
			if (Scores[i][j] < min)
			{
				min = Scores[i][j];
			}
			cout << "d" << j+1 << " : " << Scores[i][j] << endl;
			cout << endl;
		}
		cout << "The score is the minimum distance, hence Score(X" << i + 1 << ") = " << min;
		cout << endl;
	}
}
/*

Enter the number of Vectors: 4

Enter the length of each Vector: 8

A Matrix:
8 4
2 -2 -1 3
-1 3 3 -1
0 2 3 0
1 3 1 3
1 0 -1 2
-3 2 4 -1
5 -1 5 3
2 1 2 0

Eigen values:
4
5.18611
16.2895
71.7458
82.7786

4 4
-0.646744 0.130932 -0.615787 -0.430563
-0.447686 -0.699107 0.497495 -0.251643
0.256227 0.390526 0.374114 -0.801172
0.561825 -0.584461 -0.483056 -0.330779

EigenVectors:
4 4
-0.430563 -0.615787 0.130932 -0.646744
-0.251643 0.497495 -0.699107 -0.447686
-0.801172 0.374114 0.390526 0.256227
-0.330779 -0.483056 -0.584461 0.561825

C Eigen Vectors:
8 4
-0.549007 -4.04985 -0.483831 1.03113
-2.3971 3.71367 -0.472214 -0.48946
-2.9068 2.11733 -0.226636 -0.126693
-2.979 -0.198357 -3.32925 -0.0481021
-0.29095 -1.95601 -1.42852 0.220679
-2.0855 4.82186 0.355556 1.50794
-6.89937 -3.15503 1.55302 0.180574
-2.71511 0.0141476 0.343809 -1.22872


Scoring Matrix:
4 4
-3.91738 -2.28952 -7.28928 -3.00952
-5.2159 4.21393 3.16885 -4.09163
0.528445 -2.82161 1.57617 -2.3589
-1.47283 -1.01952 0.583506 1.27945

Most significant Scoring Matrix:
3 4
-3.91738 -2.28952 -7.28928 -3.00952
-5.2159 4.21393 3.16885 -4.09163
0.528445 -2.82161 1.57617 -2.3589

Euclidean Distances for X-1
d0 : 9.00731

d1 : 5.49504

d2 : 8.06999

d3 : 6.57204

The score is the minimum distance, hence Score(X1) = 5.49504

Euclidean Distances for X-2
d0 : 10.1388

d1 : 0

d2 : 6.7402

d3 : 8.34954

The score is the minimum distance, hence Score(X2) = 0

Euclidean Distances for X-3
d0 : 3.91648

d1 : 6.88105

d2 : 8.40813

d3 : 1.81294

The score is the minimum distance, hence Score(X3) = 1.81294

Euclidean Distances for X-4
d0 : 4.06446

d1 : 6.11548

d2 : 6.38345

d3 : 3.40625

The score is the minimum distance, hence Score(X4) = 3.40625
Press any key to continue . . .
*/



