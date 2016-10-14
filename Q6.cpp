#include<iostream>
#include<stdio.h>
#include<vector>
#include<math.h>

using namespace std;

void MatrixMult(vector<float>, vector<vector<float>>, vector<float> &);
void EuclideanDistance(vector<float>, vector<float>, float &);

void main()
{
	int length,m,n,t;
	float A[5][5], temp;

	cout << "\nEnter the length of vectors: ";
	cin >> length;

	vector<vector<float>> EigenV(3, vector<float>(length));
	
	
	cout << "\nEnter the number of rows and columns of the scoring matrix: ";
	cin >> m >> n;

	cout << "\nEnter the scoring matrix: \n";
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cin >> A[i][j];
		}
	}
	
	// Get Unit Eigenvectors
	for (int i = 0; i < m; i++)
	{
		cout << "\nEnter the values for Unit Eigenvector " << i << " :";
		for (int j = 0; j < length; j++)
		{
			cin >> temp;
			EigenV[i][j] = temp;
		}
	}

	cout << "\nEnter the number of Test vectors: ";
	cin >> t;

	vector<vector<float>> EigenSpace(t, vector<float>(length, 0));

	vector<vector<float>> TestV(t, vector<float>(length));

	// Get Test Vectors
	for (int i = 0; i < t; i++)
	{
		cout << "\nEnter the values for Test Vector " << i << " :";
		for (int j = 0; j < length; j++)
		{
			cin >> temp;
			TestV[i][j] = temp;
		}
	}

	//Calculate Projections
	for (int i = 0; i < t; i++)
	{
		MatrixMult(TestV[i], EigenV, EigenSpace[i]);
	}

	//Calculate scores

	vector<vector<float>> distance(t, vector<float>(n, 0));
	vector<vector<float>> TempColumn(n, vector<float>(m, 0));

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			TempColumn[i][j] = A[j][i];
		}
	}

	for (int i = 0; i< t; i++)
	{
		for (int j = 0; j < n; j++)
		{
			EuclideanDistance(TempColumn[j], EigenSpace[i], distance[i][j]);
		}
	}
	 

	//Calculate minimum distance and display the result.

	for (int i = 0; i < t; i++)
	{
		float min = distance[i][0];

		cout << "\nEuclidean Distances for X-" << i + 1 << endl;
		for (int j = 0; j < n; j++)
		{
			if (distance[i][j] < min)
			{
				min = distance[i][j];
			}
			cout << "d" << j << " : " << distance[i][j] << endl;
			cout << endl;
		}
		cout << "The score is the minimum distance, hence Score(X" << i+1 << ") = " << min;
		cout << endl;
	}

}

void MatrixMult(vector<float> A, vector<vector<float>> B, vector<float> & C)
{
		for (int j = 0; j < B.size(); j++)
		{
			for (int k = 0; k < A.size(); k++)
			{	
				C[j] += A[k] * B[j][k];
			}
		}
}

void EuclideanDistance(vector<float> A, vector<float> B, float & C)
{
	float sum = 0;

	for (int i = 0; i < A.size(); i++)
	{
		sum += pow(A[i] - B[i], 2);
	}

	C = pow(sum, 0.5);
	 
}

/*
Output:

Enter the length of vectors: 6

Enter the number of rows and columns of the scoring matrix: 3 4

Enter the scoring matrix:
3.40 4.98 3.37 4.56
-1.21 0.75 -2.88 2.22
-1.47 -0.48 1.27 1.00

Enter the values for Unit Eigenvector 0 :0.43 0.47 0.29 0.44 0.52 0.23

Enter the values for Unit Eigenvector 1 :0.04 0.50 -0.37 -0.71 0.32 0.02

Enter the values for Unit Eigenvector 2 :-0.24 -0.11 0.79 -0.44 0.13 0.30

Enter the number of Test vectors: 4

Enter the values for Test Vector 0 :2 3 1 0 3 2

Enter the values for Test Vector 1 :-4 -5 0 3 1 -2

Enter the values for Test Vector 2 :2 3 0 1 3 2

Enter the values for Test Vector 3 :3 2 1 0 3 2

Euclidean Distances for X-1
d0 : 4.36376

d1 : 2.09621

d2 : 5.24044

d3 : 0.0374165

The score is the minimum distance, hence Score(X1) = 0.0374165

Euclidean Distances for X-2
d0 : 7.0281

d1 : 9.3025

d2 : 6.46398

d3 : 9.97466

The score is the minimum distance, hence Score(X2) = 6.46398

Euclidean Distances for X-3
d0 : 3.56643

d1 : 1.16846

d2 : 5.17233

d3 : 1.31871

The score is the minimum distance, hence Score(X3) = 1.16846

Euclidean Distances for X-4
d0 : 3.92394

d1 : 1.71348

d2 : 4.79486

d3 : 0.49689

The score is the minimum distance, hence Score(X4) = 0.49689
Press any key to continue . . .
*/