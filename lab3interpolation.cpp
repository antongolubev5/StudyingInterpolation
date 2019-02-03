// lab3interpolation.cpp: ���������� ����� ����� ��� ����������� ����������.
//
//��� �� 3 ��������� �������� ������� �������� � ����� ����, ��� ����� ������� ���������� ����� ��� ����� ���(������)++++++++++++++++++++
//�������������� ������� �=const � y=� y=x^3+++++++++++++++++++++++++++++++++++++
//������� ������� ��� ����������� ������������ n=3,10,100 �� ����������� �����, �� ����������� �����++++++++++++++++++++
//����� �� ������� ��� ������ ������������ - ��������++++++++++++++++
//������� ������� �������, � ����� ����������� � ���� ����� ����� ���, ����� �� ����������� ��� �� ��������+++++++++++++++(2 �������, ������ ������� ������)
//������� ����� �� ����������� ����� ���������� �������� ������ ������������� ++++++++++++++++++++++++(��� ������������ ����� ���������)
//�������� � �������� ������ ������������� � ������� �� ������� � ��������� ��� ���������� ��� ������������ � ����� 2.2
//������ �����������++++++++++++++++++++++

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include <malloc.h>
#include <iomanip>
#include <cmath>
#include <bitset>   // ������������ ���� ������� �����
#include <iomanip>  // ��� ������������ setw()
typedef unsigned char      u8;
typedef unsigned short int u16;
typedef unsigned int       u32;
typedef unsigned long long u64;
using namespace std;

double* gaussmethod(double** extMatrix, int n)
{
	double *solution;
	solution = new double[n];
	int imax;
	double maxvalue = 0;

	for (int cnt = 0; cnt < n; cnt++)
	{
		solution[cnt] = 0;
	}
	//double det = determinant(extMatrix, n);
	//if (abs(det) < 1e-30)
	//{
	//cout << "������������ ����� 0. �� ���������� ������������� �������." << endl;
	//solution = nullptr;
	//}
	//else
	//{

	for (int i = 0; i < n - 1; i++)//���� �� �������, ������� ���������� �� �����������
	{
		//����� ���� �������� �� i-�� �������
		maxvalue = 0;
		for (int il = i; il < n; il++)
		{
			if (maxvalue < abs(extMatrix[il][i]))
			{
				maxvalue = abs(extMatrix[il][i]);
				imax = il;
			}
		}

		if (maxvalue < 1e-10)
		{
			cout << "�� ���������� ������������� �������." << endl;
			return nullptr;
		}

		if (imax != i)
		{
			double* buf = extMatrix[imax];
			extMatrix[imax] = extMatrix[i];
			extMatrix[i] = buf;
		}

		//extMatrix[i][n] = extMatrix[i][n] / extMatrix[i][i];
		double aii = extMatrix[i][i];

		if (abs(aii) < 1e-10)
		{
			cout << "�� ���������� ������������� �������. ��������� ������ ������������������� ������� - �������" << endl;
			return nullptr;
		}

		for (int j = i; j <= n; j++)//���� �� ��������� ��������, ������� ���������� �� �����������  �� i+1???
		{
			extMatrix[i][j] = extMatrix[i][j] / aii;
		}

		for (int ii = i + 1; ii < n; ii++)//��������� �� ���������� ����� i-�� ������
		{
			double a_ii_i = extMatrix[ii][i];
			for (int jj = i; jj <= n; jj++)
			{
				extMatrix[ii][jj] -= a_ii_i * extMatrix[i][jj];
			}
		}
	}
	//��������� ������ ������
	double	 aii = extMatrix[n - 1][n - 1];
	if (abs(aii) < 1e-10)
	{
		cout << "�� ���������� ������������� �������. ��������� ������ ������������������� ������� - �������" << endl;
		return nullptr;
	}
	for (int j = n - 1; j <= n; j++)//���� �� ��������� ��������, ������� ���������� �� �����������  �� i+1???
	{
		extMatrix[n - 1][j] = extMatrix[n - 1][j] / aii;
	}
	//printMatrix(extMatrix, n, n + 1, true);
	//�������� ���

	double sum = 0;
	for (int i = n - 1; i >= 0; i--)
	{
		sum = 0;
		for (int j = i + 1; j < n; j++) //��������� ��� ����� ������� ����������  ���������� �� ������������ ������� ������
		{
			sum += solution[j] * extMatrix[i][j];
		}
		solution[i] = extMatrix[i][n] - sum;//�������� �� ������ ����� 
	}

	//printMatrix(extMatrix, n);//������ ������������������� (��� ��������)
	return solution;
}

void matrix_destroyer(double** ary, int n)
{
	if (ary != nullptr)
	{
		for (int i = 0; i < n; i++) {
			delete[] ary[i];
		}
		delete[] ary;
		ary = nullptr;
	}
}

void printMatrix(double** extMatrix, int k, int m)
{
	cout << endl;
	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < m; j++)
		{
			
			cout << /*setprecision(5) << fixed <<*/ extMatrix[i][j] << " ";//simb;

		}
		cout << endl;
	}
}

void printVector(double* Vector, int k)
{
	cout << endl;
	for (int i = 0; i < k; i++)
	{
		cout << Vector[i] /*<< " " << i+1 */ << endl;
	}

}

double function1(double x) 
{
	return 1;
}

double function2(double x)
{
	return 1/(1+10*x*x);
}

double function3(double x)
{
	return 1/(atan(1+10*x*x));
}

double function4(double x)
{
	double h = 4 * x*x*x + 2 * x*x - 4 * x + 2;
	double g = asin(1/(5+x-x*x));
	 
	return (pow(h, sqrt(2)) + g - 5);
}

double function5(double x)
{
	return exp(x);
}

double** UniformNet(double a, double b, int n, double(*ff) (double))
{
    double **result;
	result = new double*[2];
	double xi, yi;

	for (int ii = 0; ii<2; ii++)
		result[ii] = new double[n];

	double h = (b - a) / ((double)n - 1.0);//��� �����

	for (int i = 0; i < n; i++)
	{
		xi = a + h*i;
		yi = ff(xi);
		//cout << xi << ' ' << yi << endl;
	
		result[0][i] = xi;
		result[1][i] = yi;
	}

	return result;
}

double** ChebyshevNet(double a, double b, int n, double(*ff) (double))
{
  
	double Pi = acos(0)*2;
	double **result;
	result = new double*[2];
	double xi, yi;

	for (int ii = 0; ii<2; ii++)
		result[ii] = new double[n];

	double h = (b - a) / ((double)n - 1.0);//��� �����

	result[0][0] = a;
	result[1][0] = ff(a);
	result[0][n-1] = b;
	result[1][n-1] = ff(b);
	for (int i = 1; i < n-1; i++)
	{
		xi = (a+b)*0.5 + (b-a)*0.5*cos(((2*i+1)*Pi)/(2*(n+1)));
		yi = ff(xi);
		
		result[0][n-1-i] = xi;
		result[1][n-1-i] = yi;
	}

	return result;
}

int GetOnesNumber(int x)
{
	bitset<1024> number;
	number = x;
	int result = 0;
	
	for (int i = 0; i < 1024; i++) {
		result += number[i];
	}
	return result;
}

u8 CountOnes4(u64 n) {
	n = ((n >> 1) & 0x5555555555555555llu) + (n & 0x5555555555555555llu);
	n = ((n >> 2) & 0x3333333333333333llu) + (n & 0x3333333333333333llu);
	n = ((n >> 4) & 0x0F0F0F0F0F0F0F0Fllu) + (n & 0x0F0F0F0F0F0F0F0Fllu);
	n = ((n >> 8) & 0x00FF00FF00FF00FFllu) + (n & 0x00FF00FF00FF00FFllu);
	n = ((n >> 16) & 0x0000FFFF0000FFFFllu) + (n & 0x0000FFFF0000FFFFllu);
	n = ((n >> 32) & 0x00000000FFFFFFFFllu) + (n & 0x00000000FFFFFFFFllu);
	return n;
}

int CountOnes0(u64 m) {
	u64 n = m;
	int res = 0;
	while (n) {
		res += n & 1;
		n >>= 1;
	}
	return res;
}

double* LagrangePolynom(double** Net, int n)
{
	bitset<64> bits;
	
	double y = 0;
	double ck;
	double p = 1;
	int k = 3;
	double* arrayp;
	arrayp = new double[n-1];
	double* arrayz;
	arrayz = new double[n];
	double coefmltpled = 1,coefKJ=0;
	
	for (int ii = 0; ii < n; ii++)
	{
		if(ii<n-1)arrayp[ii] = 0.0; 
		arrayz[ii] = 1.0;
	}

	for (int k = 0; k < n; k++)//�������� ������� �� �������� ��������
	{
		for (int node = 0; node < n; node++)// ����������� ��� �� 
		{
			if (node != k) { arrayz[k] *= (Net[0][k]-Net[0][node]); }
		}

		for (int j = 0; j < n - 1; j++)//������� �������� ��� ��
		{
			coefKJ = 0;
			for (u64 i = 0; i < pow(2, n); i++)//���� �� ����������� � ����� (��������� ������������ ��������)
			{//�������� � ������ ���� - ���� �� �������� ������� �, 0 -���������
				bits = i; coefmltpled = 1;
				if (CountOnes0(i) == j && bits[k]==0) //���������� ��������� ����� =������� �, �-� ���� �� �������� (�� ������� ��� ��)
				{
					
					for (int m = 0; m < n; m++)//���� �� �����- �����
					{
						if (bits[m] == 0&& m!=k) { coefmltpled *= -Net[0][m]; } //����������� ��� j-� ������� 
		
					}
					//������������ ������������ �� n-j �������� ��� ������� ��
					coefKJ += coefmltpled;
				}
			}
           arrayp[j] += coefKJ*Net[1][k]/ arrayz[k];//���������� �� �� ����� ������������� ��=������������ ��� j-� ������� � � �������� ��������
		}

	}
	
	return arrayp;
}

double* LagrangePolynomNew(double** Net, int n)
{
	double y = 0;
	double c = 1, cl = 0, coeff=1;
	double p = 1;
	int k = 3;
	double* arrayp;
	arrayp = new double[n - 1];
	double* arrayz;
	arrayz = new double[n];
	double* coeffl;
	coeffl = new double[n-1];
	double coefmltpled = 1, coefKJ = 0;
	int* pos;
	int* nodebit = new int[n];
	pos = new int[n];
	for (int bb = 0; bb < n; bb++)
	{
		nodebit[bb] = 1;
	}
	for (int ii = 0; ii < n; ii++)
	{
		if (ii<n - 1)arrayp[ii] = 0.0;
		arrayz[ii] = 1.0;
		pos[ii] = 0;
		coeffl[ii] = 0.0;
	}
		
	for (int k = 0; k < n; k++)
	{
		for (int node = 0; node < n; node++)// ����������� ��� �� 
		{
			if (node != k) {
				arrayz[k] /= (Net[0][k] - Net[0][node]);
			}
		}
		arrayz[k] *= Net[1][k];
	}
		//����������  ���������� ����� �������� ��������
		for (int k = 0; k < n; k++)
		{
			c = 1;
			for (int node = 0; node < n; node++)// ����������� ��� �� 
			{
				if (node != k) { c *= -Net[0][node]; }
			}
			//cl += c*Net[1][k] / arrayz[k];
			cl += c* arrayz[k];
		}
		coeffl[0] = cl;

		for (int j = 1; j < n; j++)//������� �������� ��� ����� ��������
		{
			for (int ii = 0; ii < j; ii++)
			{
				pos[ii] = ii;
			}
			//��� ��������� ���������� ����� ��������� �����������
			bool ans = false; coeff = 1;
			for (int bb = 0; bb < n; bb++)
			{
				nodebit[bb] = 1;
			}
			for (int iii = 0; iii < j; iii++)
			{
				nodebit[pos[iii]] = 0;
			}
			for (int cknumber = 0; cknumber < n; cknumber++)
			{
				
				if (nodebit[cknumber] > 0)
				{ coeff = 1;
					for (int freechlennumb = 0; freechlennumb < n; freechlennumb++)
					{
						if (cknumber != freechlennumb&&nodebit[freechlennumb] > 0)
						{
							coeff *= -Net[0][freechlennumb];
						}
						
					}
                      //coeffl[j] += coeff*Net[1][cknumber] / arrayz[cknumber];
					coeffl[j] += coeff* arrayz[cknumber];
				}
				
			}

						

			while (pos[0]<n-j && j>0)
			//���� �� ����������� j ����� (��������� ����������� ��� j-� �������)
			{
				
				
				for (int i = 1; i < j+1; i++)
				{
					
					if (pos[j - i] <= n - i && (i==1 || pos[j - i+1]==n-i+1))
					{
						
						pos[j - i] ++;//j - i ���� �� j ���������� ������
						if (i > 1) { 
							for (int ii = i-1; ii > 0; ii--)
							{
								pos[j - ii] = pos[j - ii-1] + 1;//�������� ��� ����� ������ �� j - i (�� j) ���� � ���� �������� ������
							}
						}
					
				// ��� ������� ���������� ����������� ������������
					coeff = 1;
					for (int bb = 0; bb < n; bb++)
					{nodebit[bb] = 1;}
					for (int iii = 0; iii < j; iii++)
					{	nodebit[pos[iii]] = 0;}
					for (int cknumber = 0; cknumber < n; cknumber++)
					{

						if (nodebit[cknumber] > 0)
						{
							coeff = 1;
							for (int freechlennumb = 0; freechlennumb < n; freechlennumb++)
							{
								if (cknumber != freechlennumb&&nodebit[freechlennumb] > 0)
								{
									coeff *= -Net[0][freechlennumb];
								}

							}
							//coeffl[j] += coeff*Net[1][cknumber] / arrayz[cknumber];
							coeffl[j] += coeff* arrayz[cknumber];
						}

					}//cknumber


					}//if

					

				}//i

				
			}//while

			
		}//j

	return coeffl;
}

double LagrangeInterpolation(double x, double* coeffs, int n)
{//������ �����, c���� � ������� � � ������� ����� ������� �������� ������� � �������� �����
	double y = 0;
	double xk = 1;
	double ck;

	for (int k = 0; k < n; k++)
	{
				y += coeffs[k]* xk;
				xk *= x;
	}

	return y;
}

double LagrangeInterpolation(double x, double** Net, int n)
{//������ �����, c���� � ������� � � ������� ����� ������� �������� ������� � �������� �����
	double y=0;
	double ck;

	for (int k = 0; k < n; k++)
	{
		ck = 1;
		for (int i = 0; i < n; i++)
		{
			if (i != k) { ck *= (x - Net[0][i]) / (Net[0][k] - Net[0][i]); }
		}
		y += ck*Net[1][k];
	}

	return y;
}

double SplineInterpolation(double x, double** Net, int m)
{
	int n = m - 1;
	double y = 0;
	double *a;
	a = new double[n];
	double *b;
	b = new double[n];
	double *c;
	c = new double[n];
	double *d;
	d = new double[n];

	double *g;
	g = new double[n];
	double *h;
	h = new double[n];

	double **ExtMatrix;//����������� ������� ��� ������������� �[i]
	ExtMatrix = new double*[n+1];
	for (int ii = 0; ii<n+1; ii++)
		ExtMatrix[ii] = new double[n+2];

	for (int i = 0; i < n+1; i++)//��������� ����������� ������� ������
	{
		for (int j = 0; j < n+2; j++)
		{
			ExtMatrix[i][j] = 0;
		}
	}

	double **NetNew;
	NetNew = new double*[2];
	
	for (int ii = 0; ii<2; ii++)
		NetNew[ii] = new double[n];

	for (int i = 0; i < 2; i++)//��������� ����������� ������� ������
	{
		for (int j = 0; j < n; j++)
		{
			NetNew[i][j] = Net[i][j+1];
		}
	}

	for (int i = 0; i < n; i++)
	{
		if (i != 0) { h[i] = NetNew[0][i] - NetNew[0][i - 1]; }
		else { h[i] = NetNew[0][0] - Net[0][0]; }
	}
	

	for (int i = 0; i < n; i++)
	{
		if (i != 0) { g[i] = (NetNew[1][i] - NetNew[1][i - 1])/h[i]; }
		else { g[i] = (NetNew[1][0] - Net[1][0])/h[0]; }
	}

	ExtMatrix[0][n+1] = ExtMatrix[n][n+1] = 0;

	for (int i = 1; i < n; i++)//��������� ������ ������� ����������� ������� �������������� g[i]
	{
		ExtMatrix[i][n+1] = 3 * (g[i] - g[i-1]);	
	}

	ExtMatrix[0][0] = ExtMatrix[n][n] = 1;

	for (int i = 1; i < n; i++)//���������  ���������� ����� ������� ��������������
	{
		ExtMatrix[i][i - 1] = h[i-1];
		ExtMatrix[i][i] = 2*(h[i-1]+h[i]);
		ExtMatrix[i][i + 1] = h[i];
	}

	/*cout << "������� ��: " << endl;
	printMatrix(ExtMatrix, n + 1, n + 2);*/

	double* cext = gaussmethod(ExtMatrix, n+1);
	//printVector(cext, n+1);

	for (int i = 0; i < n; i++)
	{
		b[i] = g[i] - (cext[i + 1] + 2 * cext[i])*h[i] / 3;
		d[i] = (cext[i + 1] - cext[i]) / (3 * h[i]);
		a[i] = Net[1][i];//�� ������
		c[i] = cext[i];
	}

	/*cout << "������� ����� : " << endl;
	printMatrix(ExtMatrix, n + 1, n + 2);
	cout << "cext: " << endl;
	printVector(cext, n);
	cout << "h: " << endl;
	printVector(h, n-1);
	cout << "a: " << endl;
	printVector(a, n-1);
	cout << "b: " << endl;
	printVector(b, n-1);
	cout << "c: " << endl;
	printVector(c, n-1);
	cout << "d: " << endl;
	printVector(d, n-1);
	cout << endl;*/

	//���������� �������
	int k = n-1;
	while ((x<Net[0][k])&&(k>0))
	{
		k--;
	}
	double deltax = x - Net[0][k];
	y = a[k] + b[k] * deltax + c[k] * deltax*deltax + d[k] * deltax*deltax*deltax;
	return y;
}

double** EstimateDerive(double **Net, int n, int m)
{
	//n-����, m -������� �����������, �������� ������� ����������� �� ������������ ������� y'=yi-yi-1 / xi - xi-1
	double **CurrentDeriveNet;
	CurrentDeriveNet = new double*[2];

	for (int ii = 0; ii<2; ii++)
		CurrentDeriveNet[ii] = new double[n];

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < n; j++)
		{
			CurrentDeriveNet[i][j] = Net[i][j];
		}
	}

	/*cout << "currentderivenet" << endl;
	printMatrix(CurrentDeriveNet, 2, n);*/
	if (m == 1)

	{
		for (int i = 0; i < n; i++)
		{
			if (i > 0 && i < n - 1)
			{
				CurrentDeriveNet[1][i] = (Net[1][i + 1] - Net[1][i - 1]) / (Net[0][i + 1] - Net[0][i - 1]);
			}
			else if (i == 0)//����� �����
			{
				//CurrentDeriveNet[1][i] = (Net[1][i + 1] - Net[1][i]) / (Net[0][i + 1] - Net[0][i]);
				CurrentDeriveNet[1][i] = (Net[1][i + 2] - Net[1][i]) / (Net[0][i + 2] - Net[0][i]);
				
			}
			else {//������ �����
				  //CurrentDeriveNet[1][i] = (Net[1][i] - Net[1][i - 1]) / (Net[0][i] - Net[0][i - 1]);
				CurrentDeriveNet[1][i] = (Net[1][i] - Net[1][i - 2]) / (Net[0][i] - Net[0][i - 2]);
				
			}
		}
		/*for (int i = 0; i < 2; i++)
		{
			CurrentDeriveNet[i]++;
		}*/
	}
	else

	{
		for (int i = 0; i < n; i++)
		{
			if (i > 0 && i < n - 1)
			{
				CurrentDeriveNet[1][i] = (Net[1][i + 1] - Net[1][i - 1]) / (Net[0][i + 1] - Net[0][i - 1]);
			}
			else if (i == 0)//����� �����
			{
				//CurrentDeriveNet[1][i] = (Net[1][i + 1] - Net[1][i]) / (Net[0][i + 1] - Net[0][i]);
				CurrentDeriveNet[1][i] = (Net[1][i + 2] - Net[1][i ]) / (Net[0][i + 2] - Net[0][i]);
			}
			else {//������ �����
				//CurrentDeriveNet[1][i] = (Net[1][i] - Net[1][i - 1]) / (Net[0][i] - Net[0][i - 1]);
				CurrentDeriveNet[1][i] = (Net[1][i ] - Net[1][i - 2]) / (Net[0][i] - Net[0][i - 2]);
			}
		}

		/*cout << "derived currentderivenet" << endl;
		printMatrix(CurrentDeriveNet, 2, n);*/

		/*for (int i = 0; i < 2; i++)
		{
			CurrentDeriveNet[i]++;
		}*/

		CurrentDeriveNet = EstimateDerive(CurrentDeriveNet, n, m - 1);
	}

	return CurrentDeriveNet;
}

double** EstimateDerive(double **Net, int n, int m, double(*ff) (double))
{
	//n-����, m -������� �����������, �������� ������� ����������� �� ������������ ������� y'=yi-yi-1 / xi - xi-1
	int n3 = 3 * n;
	double a3 = 3 * Net[0][0] - 2 * Net[0][n - 1];
	double b3 = 3 * Net[0][n - 1] - 2 * Net[0][0];
	double **CurrentDeriveNet = ChebyshevNet(a3, b3, n3, ff);
	double **ResultNet = new double*[2];

	for (int ii = 0; ii<2; ii++)
		ResultNet[ii] = new double[n];

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ResultNet[i][j] = Net[i][j];
		}
	}
	/*cout << "currentderivenet" << endl;
	printMatrix(CurrentDeriveNet, 2, n);*/
	if (m == 1)

	{
		for (int i = 0; i < n3; i++)
		{
			if (i > 0 && i < n3 - 1)
			{
				CurrentDeriveNet[1][i] = (Net[1][i + 1] - Net[1][i - 1]) / (Net[0][i + 1] - Net[0][i - 1]);
			}
			else if (i == 0)//����� �����
			{
				//CurrentDeriveNet[1][i] = (Net[1][i + 1] - Net[1][i]) / (Net[0][i + 1] - Net[0][i]);
				//CurrentDeriveNet[1][i] = (Net[1][i + 2] - Net[1][i]) / (Net[0][i + 2] - Net[0][i]);
				CurrentDeriveNet[1][i] = 0;
			}
			else {//������ �����
				  //CurrentDeriveNet[1][i] = (Net[1][i] - Net[1][i - 1]) / (Net[0][i] - Net[0][i - 1]);
				  //CurrentDeriveNet[1][i] = (Net[1][i] - Net[1][i - 2]) / (Net[0][i] - Net[0][i - 2]);
				CurrentDeriveNet[1][i] = 0;
			}
		}
		while (CurrentDeriveNet[0][0]<Net[0][0])
		{
			n3--; n3--;//left & right elimination
		for (int i = 0; i < 2; i++)
		{
			CurrentDeriveNet[i]++;
		}
		}
		for (int i = 0; i < n; i++)
		{
			
			ResultNet[1][i] = CurrentDeriveNet[1][i];
		}
		return (ResultNet);
	}
	else

	{
		for (int i = 0; i < n; i++)
		{
			if (i > 0 && i < n - 1)
			{
				CurrentDeriveNet[1][i] = (Net[1][i + 1] - Net[1][i - 1]) / (Net[0][i + 1] - Net[0][i - 1]);
			}
			else if (i == 0)//����� �����
			{
				//CurrentDeriveNet[1][i] = (Net[1][i + 1] - Net[1][i]) / (Net[0][i + 1] - Net[0][i]);
				CurrentDeriveNet[1][i] = (Net[1][i + 2] - Net[1][i]) / (Net[0][i + 2] - Net[0][i]);
			}
			else {//������ �����
				  //CurrentDeriveNet[1][i] = (Net[1][i] - Net[1][i - 1]) / (Net[0][i] - Net[0][i - 1]);
				CurrentDeriveNet[1][i] = (Net[1][i] - Net[1][i - 2]) / (Net[0][i] - Net[0][i - 2]);
			}
		}

		/*cout << "derived currentderivenet" << endl;
		printMatrix(CurrentDeriveNet, 2, n);*/

		for (int i = 0; i < 2; i++)
		{
			CurrentDeriveNet[i]++;
		}

		CurrentDeriveNet = EstimateDerive(CurrentDeriveNet, n - 2, m - 1);
	}
	for (int i = 0; i < n; i++)
	{

		ResultNet[1][i] = CurrentDeriveNet[1][i];
	}
	
	return  (ResultNet);
}

bool LagrangeEstimateError(double** Net, int n, double eps)
{
	double kol;
	//eps < (Mn+1)/(n+1)!*|Wn+1(x)|
	//Mn+1 - ������������ ����������� �� �������
	//Wn+1(x)=� k=0..n  (x-xk)
	double M;
	double hmax = 0;
	double h;
	//���� hmax
	for (int i = 0; i < n; i++)
	{
		h = Net[0][i+1] - Net[0][i];
		if (h > hmax) hmax = h;
	}
	//���� ��������� ����������� m-�� �������
	double** NNet = EstimateDerive(Net, n, n + 1);

	double max = -1e+50;
	for (int i = 0; i < n; i++)
	{
		if (max < abs(NNet[1][i])) { max = abs(NNet[1][i]); }
	}
	
	if (eps < max / (n + 1)*pow(hmax, n + 1)) { return false; }
	else return true;
}

bool LagrangeEstimateError(double** Net, int n, double eps, double(*ff) (double))
{
	double kol;
	//eps < (Mn+1)/(n+1)!*|Wn+1(x)|
	//Mn+1 - ������������ ����������� �� �������
	//Wn+1(x)=� k=0..n  (x-xk)
	double M;
	double hmax = 0;
	double h;
	//���� hmax
	for (int i = 0; i < n; i++)
	{
		h = Net[0][i + 1] - Net[0][i];
		if (h > hmax) hmax = h;
	}
	//���� ��������� ����������� m-�� �������
	double** NNet = EstimateDerive(Net, n, n + 1,ff);

	double max = -1e+50;
	for (int i = 0; i < n; i++)
	{
		if (max < abs(NNet[1][i])) { max = abs(NNet[1][i]); }
	}

	if (eps < max / (n + 1)*pow(hmax, n + 1)) { return false; }
	else return true;
}

int main()
{
	setlocale(LC_ALL, "Russian");
	double a = -1;
	double b = 1;
	int nkol = 50;
	int n = 30;//���������� ����� �����
	int i;
	int ntest = 50;//���������� ����� ��� ���������� �������
	double x = 0.7;
	double h, yi;
	double **unifnet;
	double **chebnet;

	unifnet = UniformNet(a, b, n, function1);
	cout << "����������� �����: " << endl;
	//�������� ����������� �����
	printMatrix(unifnet, 2, n );
	cout <<  endl;

	chebnet = ChebyshevNet(a, b, n, function1);
	cout << "����������� �����: " << endl;
	//�������� ����������� �����
	printMatrix(chebnet, 2, n);
	cout << endl;

    //��������� �������� ��� 1 ������� �^2
	/*double test1 = LagrangeInterpolation(x, unifnet, n);
	cout << "Te��: " << test1 << endl;*/

	////��������� ������� ��� 1 ������� �^2
 //   double test = SplineInterpolation(x, unifnet, n);
	//cout << "Te��: " << test << endl;

	//���������� �������

	//����������� �����
	double xi = unifnet[0][0];
	double delta = (unifnet[0][n-1] - unifnet[0][0]) / (ntest-1);

	//����������� �����
	/*double xi = chebnet[0][0];
	double delta = (chebnet[0][n - 1] - chebnet[0][0]) / (ntest - 1);*/

    //double* coeffs = LagrangePolynomNew(unifnet, n);
	 //double* coeffs = LagrangePolynomNew(chebnet, n);
	/*cout << "test" << endl;
	printVector(coeffs,n);*/

	ofstream fout;
	fout.open("MassiveDots.txt"); // ��������� ������ � ������

	for (int i = 0; i < ntest; i++)
	{
		xi = a + delta*i;
		//yi = LagrangeInterpolation(xi, coeffs, n);
		
		//yi = LagrangeInterpolation(xi, unifnet, n);
		yi = LagrangeInterpolation(xi, chebnet, n);
		//yi = SplineInterpolation(xi, unifnet, n);
		//yi = SplineInterpolation(xi, chebnet, n);

		//���������� � ���� �������� � �������� ������ ��������	
		fout << xi << " " << yi << endl; // ������ ������ � ����
	}
	fout.close(); // ��������� ����

	cout << "���� ����������� ����� ����� ����������� ����� ��� ��������� eps (���������� ������������):" << endl;
	int nn = 50;
	double eps = 1e-4;
	//while (!LagrangeEstimateError(ChebyshevNet(a, b, nn, function1), nn, eps, function1) && nn<1e+3)
	while (!LagrangeEstimateError(ChebyshevNet(a, b, nn, function1), nn, eps) && nn<1e+3)
	{
		nn++;
	}
	cout << nn << endl;

	//��������� ����������� ��� ������� ������������
	double result;
	double xx,pogr;
	double deltamax = -1;
	//double* coeffs = LagrangePolynomNew(chebnet, n);
	//printVector(coeffs, n);

	for (int i = 0; i < n; i++)
	{
		xx = unifnet[0][i];
		//pogr = LagrangeInterpolation(xtest, coeffs, n);
		pogr = SplineInterpolation(xx, unifnet, n);
		result = abs(unifnet[1][i] - pogr);
		if (abs(result) > deltamax) {deltamax = result;	}
	}
	cout << "������������ �������� ����������� = " << deltamax << endl;

	//������ ������������ ����������
	//e^x, [0,2], h =0.2, ������������ ������������ �������� �� ����������� �����, ������� � ����� x=2.2 � �������� � ������������� �������
	int nexp = 11;
	double atest = 0;
	double btest = 2;
	double xtest = 2.2;
	double htest = 0.2;
	double **unifnettest;

	unifnettest = UniformNet(atest, btest, nexp, function5);
	cout << endl;
	cout << "����������� �����: " << endl;
	//�������� ����������� �����
	printMatrix(unifnettest, 2, nexp);
	cout << endl;
	double test1 = LagrangeInterpolation(xtest, unifnettest, n);
	cout << "��������� ���������������� � ����� x = 2.2: " << test1 << endl;
	double  accurvalue = function5(xtest);
	cout << "������ �������� ������� � ����� x = 2.2: " << accurvalue << endl;
	cout << endl;
	cout << "����������� ����������� ������������ � ����� x = 2.2: " << abs(accurvalue-test1) << endl;
	//������� ����������� ������������ �� ������� �� ���������
	double theorerror;
	theorerror = 0.5*pow(htest,n+1)*(n+2)*(n+3)*function5(xtest);
	cout << "������������� ����������� ������������ � ����� x = 2.2: " << theorerror << endl;
	cout << endl;

	//������ ����������� � ����� ������������� �����

	double resulttest;
	double xxtest, pogrtest;
	double deltamaxtest = -1;
	//double* coeffs = LagrangePolynomNew(chebnet, n);
	//printVector(coeffs, n);

	//������� �����
	//xi = unifnet[0][0];
	delta = (unifnet[0][n - 1] - unifnet[0][0]) / (ntest - 1);

	for (int i = 0; i < ntest; i++) {

	xxtest = unifnet[0][0] + delta*i;
	pogrtest = SplineInterpolation(xxtest,unifnet, n);
	resulttest = abs(function1(xxtest) - pogrtest);
	if (abs(resulttest) > deltamaxtest) { deltamaxtest = resulttest; }
	cout << "x: " << xxtest << " error: " << resulttest << endl;

	}
	cout << endl;
	cout << "������������ �������� ����������� � ����� ������������� ����� = " << deltamaxtest << endl;



    matrix_destroyer(unifnet, 2);
	matrix_destroyer(chebnet, 2);
	system("pause");
    return 0;
}




//// Lab2slau.cpp: ���������� ����� ����� ��� ����������� ����������.
////
//
//#include "stdafx.h"
//#include <fstream>
//#include <iostream>
//#include <malloc.h>
//
//#include <cmath>
//using namespace std;
//typedef double mytype;
//
//mytype** readfromfile(string filename, int &n)
//{
//	char buff[500];
//	mytype element;
//	int count = 0;
//	ifstream fin(filename); // ������� ���� ��� ������		
//
//	while (!fin.eof())//������������ ����� �����
//	{
//		fin.getline(buff, 500);
//		count++;
//	}
//	//�������� ������ ��� ������
//	int i;
//	mytype **extMatrix;
//	extMatrix = new mytype*[count];
//	for (i = 0; i<count; i++)
//		//a[i], a[i] �������� � ��������� ���� double
//		extMatrix[i] = new mytype[count + 1];
//
//	fin.close(); // ��������� ����
//	fin.open(filename); // ������� ���� ��� ������
//	int elcnt = 0;
//	while (!fin.eof())
//	{
//		fin >> element; // ������� ��������� �������
//		div_t ij = div(elcnt, count + 1);
//		extMatrix[ij.quot][ij.rem] = element;
//		elcnt++;
//	}
//	n = count;
//	fin.close(); // ��������� ����
//	return extMatrix;
//}
//
//void matrix_destroyer(mytype** ary, int n)
//{
//	if (ary != nullptr)
//	{
//		for (int i = 0; i < n; i++) {
//			delete[] ary[i];
//		}
//		delete[] ary;
//		ary = nullptr;
//	}
//}
//
//void vector_destroyer(mytype* vec, int n)
//{
//	if (vec != nullptr)
//	{
//
//		delete[] vec;
//		vec = nullptr;
//	}
//}
//
//mytype** readfromscreen(int &k)
//{
//	int i, j, raz;
//	cout << "������� ����������� �������: ";
//	cin >> raz;
//	k = raz;
//	mytype **extMatrix;
//	extMatrix = new mytype*[raz];
//
//	for (i = 0; i<raz; i++)
//		extMatrix[i] = new mytype[raz + 1];
//
//	for (i = 0; i < raz; i++)
//	{
//		for (j = 0; j < raz + 1; j++)
//		{
//			cin >> extMatrix[i][j];
//		}
//	}
//	return extMatrix;
//}
//
//mytype** readfromscreensquare(int &k)
//{
//	int i, j, raz;
//	cout << "������� ����������� ������� A: ";
//	cin >> raz;
//	k = raz;
//	mytype **extMatrix;
//	extMatrix = new mytype*[raz];
//
//	for (i = 0; i<raz; i++)
//		//a[i], a[i] �������� � ��������� ���� double
//		extMatrix[i] = new mytype[raz];
//
//	for (i = 0; i < raz; i++)
//	{
//		for (j = 0; j < raz; j++)
//		{
//			cin >> extMatrix[i][j];
//		}
//	}
//	return extMatrix;
//}
//
//void printMatrix(mytype** extMatrix, int k, int m, bool extended)
//{
//	cout << endl;
//	for (int i = 0; i < k; i++)
//	{
//		for (int j = 0; j < m; j++)
//		{
//			char simb;
//			if ((j == k - 1) && extended) { simb = '='; }
//			else { simb = ' '; }
//			cout << extMatrix[i][j] << simb;
//
//		}
//		cout << endl;
//	}
//
//}
//
//void printVector(mytype* Vector, int k)
//{
//	cout << endl;
//	for (int i = 0; i < k; i++)
//	{
//		cout << Vector[i] /*<< " " << i+1 */ << endl;
//	}
//
//}
//
//mytype** multiplyMatrix(mytype** Matrix1, mytype** Matrix2, int n)
//{
//	mytype **result;
//	result = new mytype*[n];
//	mytype s = 0;
//	for (int ii = 0; ii<n; ii++)
//		result[ii] = new mytype[n];
//
//	for (int i = 0; i < n; i++)
//	{
//		for (int l = 0; l < n; l++)
//		{
//			s = 0;
//
//			for (int j = 0; j < n; j++)
//			{
//				s += Matrix1[i][j] * Matrix2[j][l];
//			}
//			result[i][l] = s;
//		}
//	}
//	return result;
//}
//
//mytype* multiplyMatrixVector(mytype** Matrix, mytype* Vector, int n)
//{
//	mytype *result;
//	result = new mytype[n];
//	mytype s = 0;
//
//	for (int i = 0; i < n; i++)
//	{
//		s = 0;
//
//		for (int j = 0; j < n; j++)
//		{
//			s += Matrix[i][j] * Vector[j];
//		}
//		result[i] = s;
//
//	}
//	return result;
//}
//
//mytype** multiplyMatrixNumber(mytype** Matrix, mytype Number, int n)
//{
//	mytype **result;
//	result = new mytype*[n];
//	for (int ii = 0; ii<n; ii++)
//		result[ii] = new mytype[n];
//
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			result[i][j] = Matrix[i][j] * Number;
//		}
//	}
//	return result;
//}
//
//mytype* multiplyVectorNumber(mytype* Vector, mytype Number, int n)
//{
//	mytype *result;
//	result = new mytype[n];
//
//	for (int i = 0; i < n; i++)
//	{
//		result[i] = Vector[i] * Number;
//	}
//	return result;
//}
//
//mytype* substractVector(mytype* Vector1, mytype* Vector2, int n)
//{
//	mytype *result;
//	result = new mytype[n];
//
//	for (int i = 0; i < n; i++)
//	{
//		result[i] = Vector1[i] - Vector2[i];
//	}
//
//	return result;
//}
//
//mytype** substractMatrix(mytype** Matrix1, mytype** Matrix2, int n)
//{
//	mytype **result;
//	result = new mytype*[n];
//
//	for (int ii = 0; ii<n; ii++)
//		result[ii] = new mytype[n];
//
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			result[i][j] = Matrix1[i][j] - Matrix2[i][j];
//		}
//	}
//	return result;
//}
//
//mytype** sumMatrix(mytype** Matrix1, mytype** Matrix2, int n)
//{
//	mytype **result;
//	result = new mytype*[n];
//
//	for (int ii = 0; ii<n; ii++)
//		result[ii] = new mytype[n];
//
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			result[i][j] = Matrix1[i][j] + Matrix2[i][j];
//		}
//	}
//	return result;
//}
//
//mytype** transposeMatrix(mytype** Matrix, int n)
//{
//	mytype **result;
//	result = new mytype*[n];
//
//	for (int ii = 0; ii<n; ii++)
//		result[ii] = new mytype[n];
//
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = i; j < n; j++)
//		{
//			result[i][j] = Matrix[j][i];
//			result[j][i] = Matrix[i][j];
//		}
//	}
//
//	return result;
//}
//
//void unitMatrix(mytype** Matrix, int n)
//{
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			if (i == j) { Matrix[i][j] = 1; }
//			else  Matrix[i][j] = 0;
//		}
//	}
//}
//
//mytype discrepancy(mytype** extMatrix, mytype* solution, int n)
//{
//	mytype *b1;
//	b1 = new mytype[n];
//	mytype result = 0;
//
//	b1 = multiplyMatrixVector(extMatrix, solution, n);
//
//	for (int i = 0; i < n - 1; i++)
//	{
//		result += pow((b1[i] - extMatrix[i][n]), 2);
//	}
//	result = sqrt(result);
//	return result;
//}
//
//mytype normVectorUnit(mytype* Vector, int n)
//{
//	mytype max = 0;
//	for (int i = 0; i < n; i++)
//	{
//		max += abs(Vector[i]);
//	}
//	return max;
//}
//
//mytype normVectorInfinity(mytype* Vector, int n)
//{
//	mytype max = 0;
//	for (int i = 0; i < n; i++)
//	{
//		if (abs(Vector[i]) > max) max = abs(Vector[i]);
//	}
//	return max;
//}
//
//mytype normMatrixUnit(mytype** Matrix, int n)//���� �������
//{
//	mytype max = 0;
//	mytype *Vector;
//	Vector = new mytype[n];
//	for (int i = 0; i < n; i++)
//	{
//		Vector[i] = 0;
//	}
//
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			Vector[j] += Matrix[i][j];
//		}
//	}
//	max = normVectorInfinity(Vector, n);
//	return max;
//}
//
//mytype normMatrixInfinity(mytype** Matrix, int n)//���� ������
//{
//	mytype max = 0;
//	mytype *Vector;
//	Vector = new mytype[n];
//	for (int i = 0; i < n; i++)
//	{
//		Vector[i] = 0;
//	}
//
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			Vector[i] += Matrix[i][j];
//		}
//	}
//	max = normVectorInfinity(Vector, n);
//	return max;
//}
//
//mytype* gaussmethod(mytype** extMatrix, int n)
//{
//	mytype *solution;
//	solution = new mytype[n];
//	int imax;
//	mytype maxvalue = 0;
//
//	for (int cnt = 0; cnt < n; cnt++)
//	{
//		solution[cnt] = 0;
//	}
//	//double det = determinant(extMatrix, n);
//	//if (abs(det) < 1e-30)
//	//{
//	//cout << "������������ ����� 0. �� ���������� ������������� �������." << endl;
//	//solution = nullptr;
//	//}
//	//else
//	//{
//
//	for (int i = 0; i < n - 1; i++)//���� �� �������, ������� ���������� �� �����������
//	{
//		//����� ���� �������� �� i-�� �������
//		maxvalue = 0;
//		for (int il = i; il < n; il++)
//		{
//			if (maxvalue < abs(extMatrix[il][i]))
//			{
//				maxvalue = abs(extMatrix[il][i]);
//				imax = il;
//			}
//		}
//
//		if (maxvalue < 1e-10)
//		{
//			cout << "�� ���������� ������������� �������." << endl;
//			return nullptr;
//		}
//
//		if (imax != i)
//		{
//			mytype* buf = extMatrix[imax];
//			extMatrix[imax] = extMatrix[i];
//			extMatrix[i] = buf;
//		}
//
//		//extMatrix[i][n] = extMatrix[i][n] / extMatrix[i][i];
//		mytype aii = extMatrix[i][i];
//
//		if (abs(aii) < 1e-10)
//		{
//			cout << "�� ���������� ������������� �������. ��������� ������ ������������������� ������� - �������" << endl;
//			return nullptr;
//		}
//
//		for (int j = i; j <= n; j++)//���� �� ��������� ��������, ������� ���������� �� �����������  �� i+1???
//		{
//			extMatrix[i][j] = extMatrix[i][j] / aii;
//		}
//
//		for (int ii = i + 1; ii < n; ii++)//��������� �� ���������� ����� i-�� ������
//		{
//			mytype a_ii_i = extMatrix[ii][i];
//			for (int jj = i; jj <= n; jj++)
//			{
//				extMatrix[ii][jj] -= a_ii_i * extMatrix[i][jj];
//			}
//		}
//	}
//	//��������� ������ ������
//	mytype	 aii = extMatrix[n - 1][n - 1];
//	if (abs(aii) < 1e-10)
//	{
//		cout << "�� ���������� ������������� �������. ��������� ������ ������������������� ������� - �������" << endl;
//		return nullptr;
//	}
//	for (int j = n - 1; j <= n; j++)//���� �� ��������� ��������, ������� ���������� �� �����������  �� i+1???
//	{
//		extMatrix[n - 1][j] = extMatrix[n - 1][j] / aii;
//	}
//	//printMatrix(extMatrix, n, n + 1, true);
//	//�������� ���
//
//	mytype sum = 0;
//	for (int i = n - 1; i >= 0; i--)
//	{
//		sum = 0;
//		for (int j = i + 1; j < n; j++) //��������� ��� ����� ������� ����������  ���������� �� ������������ ������� ������
//		{
//			sum += solution[j] * extMatrix[i][j];
//		}
//		solution[i] = extMatrix[i][n] - sum;//�������� �� ������ ����� 
//	}
//
//	//printMatrix(extMatrix, n);//������ ������������������� (��� ��������)
//	return solution;
//}
//
//mytype iterationnumbertheor(mytype eps, mytype* solution, mytype* x0, mytype normC, int n)
//{
//	mytype k;
//	mytype normsubstract = normVectorUnit((substractVector(solution, x0, n)), n);
//	k = (log(eps / normsubstract)) / (log(normC));
//	return k;
//}
//
//mytype estimateepssimpleiterationaprior(mytype* solution, mytype* x0, mytype* xk, mytype normC, int n, int k)
//{
//	mytype est;
//	mytype normsubstract = normVectorUnit((substractVector(x0, solution, n)), n);
//	est = pow(normC, k) * normsubstract;
//	return est;
//}
//
//mytype estimateepssimpleiterationaposter(mytype* solution, mytype* x0, mytype* xk, mytype normC, int n)
//{
//	mytype est;
//	mytype normsubstract = normVectorUnit((substractVector(solution, xk, n)), n);
//	est = (normC / (1 - normC)) *normsubstract;
//	return est;
//}
//
//mytype estimateepsrelaxationaprior(mytype* solution, mytype* x0, mytype* xk, mytype normCu, mytype normCl, int n, int k)
//{
//	mytype est;
//	mytype normsubstract = normVectorUnit((substractVector(x0, xk, n)), n);
//	mytype q = normCu / (1 - normCl);
//	est = pow(q, k) * normsubstract;
//	return est;
//}
//
//mytype estimateepsrelaxationaposter(mytype* solution, mytype* x0, mytype* xk, mytype normC, mytype normCu, int n)
//{
//	mytype est;
//	mytype normsubstract = normVectorUnit((substractVector(solution, xk, n)), n);
//	est = (normCu / (1 - normC)) *normsubstract;
//	return est;
//}
//
//mytype* simpleiteration(mytype** extMatrix, mytype* x0, mytype eps, int n)
//
//{// ��� 1-� ������� eps=1�-10 tau= 0.000512 -� 1�-6 �������� �� 2 � ����� 200 ��� (����� 1)
// //mytype *x;//true solution
// /*x = new mytype[n];
// x[0] = 5;
// x[1] = -7;
// x[2] = 12;
// x[3] = 4;*/
//	mytype *vectokill;
//	mytype *solution;
//	solution = new mytype[n];
//	mytype *b;
//	b = new mytype[n];
//	mytype *y;
//	y = new mytype[n];
//	mytype **A;
//	A = new mytype*[n];
//	mytype **boofA = nullptr;
//	boofA = new mytype*[n];
//	mytype **C;
//	C = new mytype*[n];
//	mytype **E;
//	E = new mytype*[n];
//	mytype **memtokill;
//	mytype iterparam = 1e-9;
//	mytype normCmin = 1e+60;
//	mytype normC, iterparammin;
//
//	//diag norm
//	mytype diagii = 0;
//	for (int str = 0; str < n; str++)
//	{
//		diagii = extMatrix[str][str];
//		for (int col = 0; col < n + 1; col++)
//		{
//			extMatrix[str][col] = extMatrix[str][col] / diagii;
//		}
//	}
//	cout << "��������������� �������:" << endl;
//	printMatrix(extMatrix, n, n + 1, true);
//
//	for (int ii = 0; ii < n; ii++)
//	{
//		A[ii] = new mytype[n];
//		C[ii] = new mytype[n];
//		E[ii] = new mytype[n];
//	}
//
//	for (int i = 0; i < n; i++)//��������� ������� � ������
//	{
//		for (int j = 0; j < n; j++)
//		{
//			E[i][j] = 0;
//		}
//	}
//
//	unitMatrix(E, n);//�������������� ��������� ������� �
//
//	for (int i = 0; i < n; i++)//����������� ��������� ���� �� �������  � ������ 
//	{
//		b[i] = extMatrix[i][n];
//	}
//
//	for (int i = 0; i < n; i++)//��������� �������� �� �������  � ���������� ������� 
//	{
//		for (int j = 0; j < n; j++)
//		{
//			A[i][j] = extMatrix[i][j];
//		}
//	}
//	boofA = nullptr;
//	//������� ����� ������������ ��������, ��� ������� ����� ������� � ����� ����������
//	for (int i = 1; i < 200; i++)
//	{
//		//��������� ������� � �� ������� �=-(tau*A-E)
//		//matrix_destroyer(mytype** ary, int n)
//		memtokill = boofA;
//		boofA = multiplyMatrixNumber(A, iterparam, n);
//		matrix_destroyer(memtokill, n);
//		memtokill = boofA;
//		boofA = substractMatrix(boofA, E, n);
//		matrix_destroyer(memtokill, n);
//		memtokill = C;
//		C = multiplyMatrixNumber(boofA, -1, n);
//		matrix_destroyer(memtokill, n);
//		/*normC = normMatrixInfinity(C, n);*/
//		normC = normMatrixUnit(C, n);
//		if (normC < normCmin)
//		{
//			normCmin = normC;
//			iterparammin = iterparam;
//		}
//
//		/*printMatrix(C, n, n, false);
//		cout << endl;*/
//		iterparam *= 2.0;
//	}
//
//	//iterparammin = 0.1;
//	//������� ������� � � ������ �, ������ �� ���������� ������������� ���������
//	memtokill = boofA;
//	boofA = multiplyMatrixNumber(A, iterparammin, n);
//
//	boofA = substractMatrix(boofA, E, n);
//	C = multiplyMatrixNumber(boofA, -1, n);
//	matrix_destroyer(memtokill, n);
//	y = multiplyVectorNumber(b, iterparammin, n);
//
//	cout << endl;
//	cout << "����� inf ������� �:" << normCmin << endl;
//	cout << "����� unit ������� �:" << normMatrixUnit(C, n) << endl;
//	eps *= abs((1 - normCmin) / normCmin);
//	cout << "tau = " << iterparammin << endl;
//	cout << endl;
//	cout << "������� �: " << endl;
//	printMatrix(C, n, n, false);
//	cout << endl;
//	cout << "������ �: " << endl;
//	printVector(y, n);
//
//	for (int i = 0; i < n; i++)// solution = x0
//	{
//		solution[i] = x0[i];
//	}
//
//	mytype* xk = multiplyVectorNumber(solution, 10, n);
//	xk[0] += 10;
//	mytype* minusy = multiplyVectorNumber(y, -1, n);// y = -y
//	int cnt = 0;
//	while (normVectorUnit(substractVector(solution, xk, n), n) > eps&&cnt<1e7)
//	{
//		for (int w = 0; w < n; w++)
//		{
//			xk[w] = solution[w];
//		}
//
//		vectokill = solution;
//		solution = multiplyMatrixVector(C, solution, n);
//		vector_destroyer(vectokill, n);
//		vectokill = solution;
//		solution = substractVector(solution, minusy, n);
//		vector_destroyer(vectokill, n);
//		/*cout << "iteration = " << cnt << endl;
//		cout <<"xk+1 = " << endl;
//		printVector(solution, n);
//		cout << "xk = " << endl;
//		printVector(xk, n);
//		cout << "xk+1-xk " << endl;
//		cout << normVectorUnit(substractVector(solution, xk, n), n) << endl;
//		cout << "xk+1-xtrue " << endl;
//		//cout << normVectorUnit(substractVector(solution, x, n), n) << endl;*/
//		cnt++;
//	}
//	cout << endl;
//	cout << "���������� ��������: " << cnt << endl;
//	cout << endl;
//	mytype kol = round(iterationnumbertheor(eps, solution, x0, normCmin, n));
//	cout << "������������� ������ ���������� ��������: " << kol << endl;
//	cout << endl;
//	mytype est = estimateepssimpleiterationaprior(solution, x0, xk, normCmin, n, cnt);
//	cout << "��������� ����������� ������� = " << est << endl;
//	cout << endl;
//
//	est = estimateepssimpleiterationaposter(solution, x0, xk, normCmin, n);
//	cout << "������������� ����������� ������� = " << est << endl;
//	return solution;
//}
//
//mytype* jacobimethod(mytype** extMatrix, mytype* x0, mytype eps, int n)
//{
//	mytype normC;
//	mytype *vectokill;
//	mytype *solution;
//	solution = new mytype[n];
//	mytype *b;
//	b = new mytype[n];
//	mytype *y;
//	y = new mytype[n];
//
//	mytype **A;
//	A = new mytype*[n];
//	mytype **boofA = nullptr;
//	boofA = new mytype*[n];
//	mytype **C;
//	C = new mytype*[n];
//	mytype **E;
//	E = new mytype*[n];
//	mytype **L;
//	L = new mytype*[n];
//	mytype **D;
//	D = new mytype*[n];
//	mytype **U;
//	U = new mytype*[n];
//
//	for (int ii = 0; ii < n; ii++)
//	{
//		A[ii] = new mytype[n];
//		C[ii] = new mytype[n];
//		E[ii] = new mytype[n];
//		L[ii] = new mytype[n];
//		D[ii] = new mytype[n];
//		U[ii] = new mytype[n];
//	}
//	//�������������� ������� �, ������� �� ����������� ������� ������� � � ������� b ��������������
//	for (int i = 0; i < n; i++)//��������� ����� ������� ������
//	{
//		for (int j = 0; j < n; j++)
//		{
//			E[i][j] = 0;
//			L[i][j] = 0;
//			D[i][j] = 0;
//			U[i][j] = 0;
//		}
//	}
//
//	unitMatrix(E, n);//�������������� ��������� ������� �
//
//	for (int i = 0; i < n; i++)//����������� ��������� ���� �� �������  � ������ 
//	{
//		b[i] = extMatrix[i][n];
//	}
//
//	for (int i = 0; i < n; i++)//��������� �������� �� �������  � ���������� ������� 
//	{
//		for (int j = 0; j < n; j++)
//		{
//			A[i][j] = extMatrix[i][j];
//		}
//	}
//	boofA = nullptr;
//	//������������ ������� � � ���� ���� ������: ����������������� L, ������������ D  � ���������������� U
//
//	for (int i = 0; i < n; i++)//���������������� L
//	{
//		for (int j = 0; j < n; j++)
//		{
//			if (i > j) { L[i][j] = A[i][j]; }
//		}
//	}
//
//	/*cout << "������� L: " << endl;
//	printMatrix(L, n, n, false);
//	cout << endl;*/
//
//	for (int i = 0; i < n; i++)//������������ D
//	{
//		for (int j = 0; j < n; j++)
//		{
//			if (i == j) { D[i][j] = A[i][j]; }
//		}
//	}
//
//	/*cout << "������� D: " << endl;
//	printMatrix(D, n, n, false);
//	cout << endl;*/
//
//	for (int i = 0; i < n; i++)//����������������� U
//	{
//		for (int j = 0; j < n; j++)
//		{
//			if (i < j) { U[i][j] = A[i][j]; }
//		}
//	}
//
//	/*cout << "������� U: " << endl;
//	printMatrix(U, n, n, false);
//	cout << endl;*/
//
//	//��������� ������� � ����������
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			if (i != j) { C[i][j] = -(A[i][j] / A[i][i]); }
//			else (C[i][j] = 0);
//		}
//	}
//
//	cout << "������� �: " << endl;
//	printMatrix(C, n, n, false);
//	cout << endl;
//	normC = normMatrixUnit(C, n);
//	cout << "����� unit ������� �:" << normC << endl;
//	cout << "����� inf ������� �:" << normMatrixUnit(C, n) << endl;
//	cout << endl;
//	eps *= abs((1 - normC) / normC);
//
//	//��������� ������� y ����������
//	for (int i = 0; i < n; i++)
//	{
//		y[i] = b[i] / A[i][i];
//	}
//
//	cout << "������ �: " << endl;
//	printVector(y, n);
//
//	for (int i = 0; i < n; i++)// solution = x0
//	{
//		solution[i] = x0[i];
//	}
//
//	mytype* xk = multiplyVectorNumber(solution, 10, n);
//	xk[0] += 10;
//	mytype* minusy = multiplyVectorNumber(y, -1, n);// y = -y
//	int cnt = 0;
//
//	while (normVectorUnit(substractVector(solution, xk, n), n) > eps&&cnt<1e7)
//	{
//		for (int w = 0; w < n; w++)
//		{
//			xk[w] = solution[w];
//		}
//
//		vectokill = solution;
//		solution = multiplyMatrixVector(C, solution, n);
//		vector_destroyer(vectokill, n);
//		vectokill = solution;
//		solution = substractVector(solution, minusy, n);
//		vector_destroyer(vectokill, n);
//		/*cout << "iteration = " << cnt << endl;
//		cout <<"xk+1 = " << endl;
//		printVector(solution, n);
//		cout << "xk = " << endl;
//		printVector(xk, n);
//		cout << "xk+1-xk " << endl;
//		cout << normVectorUnit(substractVector(solution, xk, n), n) << endl;
//		cout << "xk+1-xtrue " << endl;
//		//cout << normVectorUnit(substractVector(solution, x, n), n) << endl;*/
//		cnt++;
//	}
//	cout << endl;
//	cout << "���������� ��������: " << cnt << endl;
//	cout << endl;
//	mytype kol = round(iterationnumbertheor(eps, solution, x0, normC, n));
//	cout << "������������� ������ ���������� ��������: " << kol << endl;
//	cout << endl;
//
//	mytype est = estimateepssimpleiterationaprior(solution, x0, xk, normC, n, cnt);
//	cout << "��������� ����������� ������� = " << est << endl;
//	cout << endl;
//
//	est = estimateepssimpleiterationaposter(solution, x0, xk, normC, n);
//	cout << "������������� ����������� ������� = " << est << endl;
//
//	return solution;
//}
//
//mytype* relaxationmethoddiag(mytype* vecta, mytype* vectb, mytype* vectc, mytype* vectd, mytype* x0, mytype eps, int n)
//{
//	mytype *solution;
//	solution = new mytype[n];
//	mytype relaxparam = 1;
//
//	for (int i = 0; i < n; i++)// solution = x0
//	{
//		solution[i] = x0[i];
//	}
//
//	mytype* xk = multiplyVectorNumber(solution, 10, n);
//	xk[0] += 10;//xk!=x0
//	int cnt = 0;
//
//	while (normVectorUnit(substractVector(solution, xk, n), n) > eps&&cnt<1e7)
//	{
//		for (int w = 0; w < n; w++)
//		{
//			xk[w] = solution[w];
//		}
//
//		solution[0] = (1 - relaxparam)*xk[0] - relaxparam*vectc[0] / vectb[0] * xk[1] + relaxparam*vectd[0] / vectb[0];
//
//		for (int i = 1; i < n - 1; i++)
//		{
//			solution[i] = -relaxparam*vecta[i] / vectb[i] * solution[i - 1] + (1 - relaxparam)*xk[i]
//				- relaxparam*vectc[i] / vectb[i] * xk[i + 1] + relaxparam*vectd[i] / vectb[i];
//		}
//
//		solution[n - 1] = -relaxparam*vecta[n - 2] / vectb[n - 1] * solution[n - 2] + (1 - relaxparam)*xk[n - 1]
//			+ relaxparam*vectd[n - 1] / vectb[n - 1];
//
//		cnt++;
//		//printVector(solution, n);
//	}
//
//	cout << endl;
//	cout << "���������� �������� = " << cnt << endl;
//
//	return solution;
//}
//
//mytype** reverseMatrix(mytype** Matrix, int n)
//{
//	mytype **extMatrix;
//	extMatrix = new mytype*[n];
//	mytype **result;
//	result = new mytype*[n];
//	mytype *ordinary;
//	ordinary = new mytype[n];
//
//	for (int ii = 0; ii<n; ii++)
//		result[ii] = new mytype[n];
//	for (int ii = 0; ii<n; ii++)
//		extMatrix[ii] = new mytype[n + 1];
//
//	for (int j = 0; j < n; j++)//��������� ������ ������� extMatrix ������
//	{
//		extMatrix[j][n] = 0;
//	}
//
//	for (int i = 0; i < n; i++)
//	{
//		for (int ii = 0; ii < n; ii++)//��������� �������� �� ������� � ����������� �������
//		{
//			for (int jj = 0; jj < n; jj++)
//			{
//				extMatrix[ii][jj] = Matrix[ii][jj];
//			}
//		}
//		for (int jj = 0; jj < n; jj++)//��������� ������ ������� extMatrix ������
//		{
//			extMatrix[jj][n] = 0;
//		}
//
//		extMatrix[i][n] = 1;
//
//		ordinary = gaussmethod(extMatrix, n);
//		if (ordinary == nullptr)
//		{
//			cout << "�� ���������� �������� �������." << endl;
//			return nullptr;
//		}
//		for (int j = 0; j < n; j++)
//		{
//			result[j][i] = ordinary[j];
//		}
//
//	}
//
//	return result;
//}
//
//mytype* relaxationmethodgeneral(mytype** extMatrix, mytype* x0, mytype eps, int n)
//{
//	mytype relaxparam = 1;
//	mytype normCu, normC, normCl;
//	mytype *solution;
//	solution = new mytype[n];
//	mytype **A;
//	A = new mytype*[n];
//	mytype **E;
//	E = new mytype*[n];
//	mytype *b;
//	b = new mytype[n];
//	mytype **L;
//	L = new mytype*[n];
//	mytype **D;
//	D = new mytype*[n];
//	mytype **U;
//	U = new mytype*[n];
//	mytype **C;
//	C = new mytype*[n];
//	mytype **Cu;
//	Cu = new mytype*[n];
//	mytype **Cl;
//	Cl = new mytype*[n];
//	mytype **boofMatrix1;
//	boofMatrix1 = new mytype*[n];
//	mytype **boofMatrix2;
//	boofMatrix2 = new mytype*[n];
//	mytype **boofMatrix3;
//	boofMatrix3 = new mytype*[n];
//
//	for (int ii = 0; ii < n; ii++)
//	{
//		A[ii] = new mytype[n];
//		E[ii] = new mytype[n];
//		C[ii] = new mytype[n];
//		Cu[ii] = new mytype[n];
//		Cl[ii] = new mytype[n];
//		L[ii] = new mytype[n];
//		D[ii] = new mytype[n];
//		U[ii] = new mytype[n];
//		boofMatrix1[ii] = new mytype[n];
//		boofMatrix2[ii] = new mytype[n];
//		boofMatrix3[ii] = new mytype[n];
//	}
//
//	unitMatrix(E, n);//�������������� ��������� ������� �
//
//					 //�������� ����������� ������� �� ���������� � �������
//	for (int i = 0; i < n; i++)//��������� �������� �� �������  � ���������� ������� 
//	{
//		for (int j = 0; j < n; j++)
//		{
//			A[i][j] = extMatrix[i][j];
//		}
//	}
//
//	for (int i = 0; i < n; i++)//����������� ��������� ���� �� �������  � ������ 
//	{
//		b[i] = extMatrix[i][n];
//	}
//
//	for (int i = 0; i < n; i++)//��������� ����� ������� ������
//	{
//		for (int j = 0; j < n; j++)
//		{
//			L[i][j] = 0;
//			D[i][j] = 0;
//			U[i][j] = 0;
//			Cu[i][j] = 0;
//			Cl[i][j] = 0;
//			boofMatrix1[i][j] = 0;
//			boofMatrix2[i][j] = 0;
//			boofMatrix3[i][j] = 0;
//		}
//	}
//
//	//������������ ������� � � ���� ���� ������: ���������������� L, ������������ D  � ����������������� U
//
//	for (int i = 0; i < n; i++)//���������������� L
//	{
//		for (int j = 0; j < n; j++)
//		{
//			if (i > j) { L[i][j] = A[i][j]; }
//		}
//	}
//
//	cout << "������� L: " << endl;
//	printMatrix(L, n, n, false);
//	cout << endl;
//
//	for (int i = 0; i < n; i++)//������������ D
//	{
//		for (int j = 0; j < n; j++)
//		{
//			if (i == j) { D[i][j] = A[i][j]; }
//		}
//	}
//
//	cout << "������� D: " << endl;
//	printMatrix(D, n, n, false);
//	cout << endl;
//
//	for (int i = 0; i < n; i++)//����������������� U
//	{
//		for (int j = 0; j < n; j++)
//		{
//			if (i < j) { U[i][j] = A[i][j]; }
//		}
//	}
//
//	cout << "������� U: " << endl;
//	printMatrix(U, n, n, false);
//	cout << endl;
//	//������� ������� � �� ���������� ���������
//	boofMatrix1 = reverseMatrix(D, n);
//	boofMatrix1 = multiplyMatrixNumber(boofMatrix1, relaxparam, n);
//	boofMatrix1 = multiplyMatrix(boofMatrix1, L, n);
//	boofMatrix1 = sumMatrix(E, boofMatrix1, n);
//	boofMatrix2 = multiplyMatrixNumber(boofMatrix2, 1 - relaxparam, n);
//	boofMatrix2 = multiplyMatrix(boofMatrix2, E, n);
//	boofMatrix3 = reverseMatrix(D, n);
//	boofMatrix3 = multiplyMatrixNumber(boofMatrix3, relaxparam, n);
//	boofMatrix3 = multiplyMatrix(boofMatrix3, U, n);
//	boofMatrix2 = substractMatrix(boofMatrix2, boofMatrix3, n);
//	boofMatrix1 = reverseMatrix(boofMatrix1, n);
//	boofMatrix3 = multiplyMatrix(boofMatrix3, boofMatrix1, n);
//
//	normC = normMatrixUnit(boofMatrix3, n);
//	//�������� ����������������� ������� �u
//	for (int i = 0; i < n; i++)//����������������� U
//	{
//		for (int j = 0; j < n; j++)
//		{
//			if (i < j) { Cu[i][j] = boofMatrix3[i][j]; }
//		}
//	}
//
//	normCu = normMatrixUnit(Cu, n);
//
//	//�������� ���������������� ������� �l
//	for (int i = 0; i < n; i++)//���������������� L
//	{
//		for (int j = 0; j < n; j++)
//		{
//			if (i > j) { Cl[i][j] = boofMatrix3[i][j]; }
//		}
//	}
//
//	normCl = normMatrixUnit(Cu, n);
//
//	for (int i = 0; i < n; i++)// solution = x0
//	{
//		solution[i] = x0[i];
//	}
//
//	mytype* xk = multiplyVectorNumber(solution, 10, n);
//	xk[0] += 10;//xk!=x0
//	int cnt = 0;
//	mytype sum1, sum2;
//	eps *= abs((1 - normC) / normCu);
//
//	while (normVectorUnit(substractVector(solution, xk, n), n) > eps&&cnt<1e7)
//	{
//		for (int w = 0; w < n; w++)
//		{
//			xk[w] = solution[w];
//		}
//
//		//solution[0] = (1 - relaxparam)*xk[0] - relaxparam*vectc[0] / vectb[0] * xk[1] + relaxparam*vectd[0] / vectb[0];
//
//		for (int i = 0; i < n; i++)
//		{
//			sum1 = sum2 = 0;
//			for (int j = 0; j < n; j++)
//			{
//				if (j < i) { sum1 += A[i][j] / A[i][i] * solution[j]; }
//				if (j > i) { sum2 += A[i][j] / A[i][i] * xk[j]; }
//			}
//			solution[i] = -relaxparam*sum1 + (1 - relaxparam)*xk[i] - relaxparam*sum2 + relaxparam*b[i] / A[i][i];
//		}
//
//		//solution[n - 1] = -relaxparam*vecta[n - 2] / vectb[n - 1] * solution[n - 2] + (1 - relaxparam)*xk[n - 1]
//		//+ relaxparam*vectd[n - 1] / vectb[n - 1];
//
//		cnt++;
//		//printVector(solution, n);
//	}
//
//	cout << "���������� �������� = " << cnt << endl;
//	cout << endl;
//	mytype kol = round(iterationnumbertheor(eps, solution, x0, normC, n));
//	cout << "������������� ������ ���������� ��������: " << kol << endl;
//	cout << endl;
//
//	mytype est = estimateepsrelaxationaprior(solution, x0, xk, normCu, normCl, n, cnt);
//	cout << "��������� ����������� ������� = " << est << endl;
//	cout << endl;
//
//	est = estimateepsrelaxationaposter(solution, x0, xk, normC, normCu, n);
//	cout << "������������� ����������� ������� = " << est << endl;
//	return solution;
//}
//
//int main()
//{
//	setlocale(LC_ALL, "rus");
//	int k;
//	mytype eps = 1e-10;
//
//	mytype** extMatrix = readfromfile("test1.txt", k);
//	mytype *x0;
//	x0 = new mytype[k];
//	for (int i = 0; i < k; i++)
//	{
//		x0[i] = 1;
//	}
//
//	/*x0[0] = 5;
//	x0[1] = -7;
//	x0[2] = 12;
//	x0[3] = 4;*/
//
//	x0[0] = 1;
//	x0[1] = 1;
//	x0[2] = 1;
//	x0[3] = 1;
//
//	cout << "������ ���������� �����������: " << endl;
//	printVector(x0, k);;
//	cout << endl;
//
//	cout << "������� ���� ������� ������� ��������: " << endl;
//	printMatrix(extMatrix, k, k + 1, true);
//	cout << endl;
//	mytype *solution = simpleiteration(extMatrix, x0, eps, k);
//
//	cout << endl;
//	cout << "�������: " << endl;
//	printVector(solution, k);
//	cout << endl;
//
//	mytype dis�rep = discrepancy(extMatrix, solution, k);
//	cout << "����� ������� ������� = " << dis�rep << endl;
//	cout << endl;
//
//	/*cout << "������� ���� ������� �����: " << endl;
//	printMatrix(extMatrix, k, k + 1, true);
//	cout << endl;
//	mytype *solution = jacobimethod(extMatrix, x0, eps, k);
//
//	cout << endl;
//	cout << "�������: " << endl;
//	printVector(solution, k);
//	cout << endl;
//
//	mytype dis�rep = discrepancy(extMatrix, solution, k);
//	cout << "����� ������� ������� = " << dis�rep << endl;
//	cout << endl;*/
//
//	//��������� 4 ������� ��� ���������������� ������� (��� ������� ������� � ����������)
//	k = 205;//����������� �� �������
//
//	mytype *vecta;//������������
//	vecta = new mytype[k - 1];
//	mytype *vectb;//���������
//	vectb = new mytype[k];
//	mytype *vectc;//������������
//	vectc = new mytype[k - 1];
//	mytype *vectd;//������ �����
//	vectd = new mytype[k];
//
//	mytype  *x0large;//������ ���������� �����������
//	x0large = new mytype[k];
//	for (int i = 0; i < k; i++)
//	{
//		x0large[i] = 1;
//		//x0large[i] = 2 - ((i + 1) % 2);
//	}
//
//	for (int i = 0; i < k - 1; i++)
//	{
//		vecta[i] = vectc[i] = 1;
//	}
//
//
//	for (int i = 0; i < k; i++)
//	{
//		vectb[i] = 4;
//
//		if (i == 0)
//		{
//			vectd[i] = 6;
//		}
//		else  if (i == k - 1)
//		{
//			vectd[i] = 9 - 3 * ((k % 2));
//		}
//		else
//		{
//			vectd[i] = 10 - 2 * ((i + 1) % 2);
//		}
//	}
//
//	/* printVector(vecta, k-1);
//	cout << endl;
//	printVector(vectb, k);
//	cout << endl;
//	printVector(vectc, k - 1);
//	cout << endl;
//	printVector(vectd, k);
//	cout << endl;*/
//
//	/*cout << endl;
//	cout << "������� ���������������� ���� ������� ����������� ������� ����������: " << endl;
//	mytype *solution = relaxationmethoddiag(vecta, vectb, vectc, vectd, x0large, eps, k);
//
//	cout << endl;
//	cout << "�������: " << endl;
//	printVector(solution, k);
//	cout << endl;*/
//
//	/*cout << "������� ���� ������� ����������: " << endl;
//	extMatrix = readfromfile("test1.txt", k);
//	printMatrix(extMatrix, k, k + 1, true);
//	cout << endl;
//	mytype *solution = relaxationmethodgeneral(extMatrix, x0, eps, k);
//
//	cout << "�������: " << endl;
//	printVector(solution, k);
//	cout << endl;*/
//	system("pause");
//	return 0;
//}


