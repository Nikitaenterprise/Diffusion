
#include "stdafx.h"

int main()
{
	int Ni = 100, Nj = 200, Nk = 400, l = 1000, t = 3000;
	double Da = 5, Db = 1, Dab = Da/Db;
	double h = 2 * static_cast<double>(l) / static_cast<double>(Ni), tau = static_cast<double>(t) / static_cast<double>(Nk), gamma = tau / pow(h, 2);
	double ***C, ***D;

	std::cout << gamma << std::endl;
	
	//Massives init and set to 0
	C = new double**[Ni];
	D = new double**[Ni];
	for (int i = 0; i < Ni; i++)
	{
		C[i] = new double*[Nj];
		D[i] = new double*[Nj];
		for (int j = 0; j < Nj; j++)
		{
			C[i][j] = new double[Nk];
			D[i][j] = new double[Nk];
			for (int k = 0; k < Nk; k++)
			{
				C[i][j][k] = 0;
				D[i][j][k] = 0;
			}
		}
	}
	//Setting up borders
	for (int i = Ni / 2; i < Ni; i++)
	{
		for (int j = 0; j < Nj; j++)
		{
			C[i][j][0] = 1;
		}
	}
	for (int i = 0; i < Ni; i++)
	{
		for (int j = Nj / 2; j < Nj; j++)
		{
			C[i][j][0] = 1;
		}
	}

	for (int i = 0; i < Ni; i++)
	{
		for (int j = 0; j < Nj; j++)
		{
			for (int k = 1; k < Nk; k++)
			{
				C[Ni - 1][j][k] = 1;
				C[i][Nj - 1][k] = 1;
			}
		}
	}

	std::ofstream out0("C.txt");
	for (int i = 0; i < Ni; i++)
	{
		for (int j = 0; j < Nj; j++)
		{
			out0 << i << "\t" << j << "\t" << C[i][j][0] << std::endl;
		}
	}
	out0.close();
	
	for (int k = 1; k < Nk; k++)
	{
		for (int i = 1; i < Ni; i++)
		{
			for (int j = 1; j < Nj; j++)
			{
				D[i][j][k] = Dab / (C[i][j][k] * (Dab - 1) + 1 / Dab);
			}
		}
	}

	//Start
	std::ofstream out1("D_ij.txt");
	std::ofstream out2("C_ij.txt");
	try
	{
		for (int k = 1; k < Nk; k++)
		{
			//std::cout << D[Ni/2][Nj/2][k] << "\t" << C[Ni / 2][Nj / 2][k] << std::endl;
			for (int i = 2; i < Ni - 2; i += 2)
			{
				for (int j = 2; j < Nj - 2; j += 2)
				{
					//std::cout << "C_s = " << C[i][j][k];
					C[i][j][k] = C[i][j][k - 1] +
						gamma / 2 *
						(D[i + 1][j][k - 1] *
						(C[i + 2][j][k - 1] - C[i][j][k - 1]) -
							D[i - 1][j][k - 1] *
							(C[i][j][k - 1] - C[i - 2][j][k - 1])) /
							(2 * pow(h, 2))

						+

						gamma / 2 *
						(D[i][j + 1][k - 1] *
						(C[i][j + 2][k - 1] - C[i][j][k - 1]) -
							D[i][j - 1][k - 1] *
							(C[i][j][k - 1] - C[i][j - 2][k - 1])) /
							(2 * pow(h, 2));
					//std::cout << "\tC_e = " << C[i][j][k] << std::endl;
					D[i][j][k] = Dab / (C[i][j][k] * (Dab - 1) + 1 / Dab);
					if (k == Nk - 2)
					{
						out1 << i << "\t" << j << "\t" << D[i][j][k] << std::endl;
						out2 << i << "\t" << j << "\t" << C[i][j][k] << std::endl;
					}
				}
				//Watching for gamma
				if (gamma > 0.25) throw 1;

				for (int i = 1; i < Ni - 1; i += 2)
				{
					for (int j = 1; j < Nj - 1; j += 2)
					{
						D[i][j][k] = (D[i - 1][j][k] + D[i + 1][j][k]) / 2 + (D[i][j - 1][k] + D[i][j + 1][k]) / 2;
					}
				}
			}
		}
	}
	catch (int a) { if (a == 1) std::cout << "Gamma = " << gamma << std::endl; }
	catch (...) { std::cout << "Something wrong" << std::endl; }
	out1.close();
	out2.close();
	system("pause");
    return 0;
}
