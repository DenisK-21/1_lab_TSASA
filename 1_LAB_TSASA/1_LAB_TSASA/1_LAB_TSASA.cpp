#include <iostream>
#include<cmath>
#include<vector>
#include <iomanip>

double Function(const double& x)
{
	//double y = -sqrt(x) * sin(x) + 2;
	double y = cos(x) * tanh(x);
	return y;
}
double fun_min(std::vector<std::pair<double, double>>& Сoordinate_point, int  N)
{
	std::pair<double, double> min = Сoordinate_point[0];
	for (int i = 0; i < N; i++)
	{
		if (Сoordinate_point[i].second < min.second)
			min = Сoordinate_point[i];
	}
	return min.first;
}
int main()
{
	//optimal passive search
	const double a = 1.5, b = 4;
	int N = 1;
	double Search_accuracy = 1;
	std::cout << "  Optimal passive search  " << std::endl;
	std::cout << "_________________________________" << std::endl;
	std::cout << "| Number of |\t" << "Value  x\t|" << std::endl;
	std::cout << "| points N  |\t" << "at a minimum  \t|" << std::endl;
	std::cout << "---------------------------------" << std::endl;
	while (Search_accuracy > 0.1)
	{
		std::vector < std::pair<double, double>> Сoordinate_point(N);
		for (size_t k = 1; k <= N; k++)
		{
			Сoordinate_point[k - 1].first = (b - a) / (N + 1) * k + a;
			Сoordinate_point[k - 1].second = Function(Сoordinate_point[k - 1].first);
		}
		Search_accuracy = (b - a) / (N + 1);
		std::cout << "|   " << N << "\t" << "    |\t" << std::setprecision(3)
			<< fun_min(Сoordinate_point, N) << "+-" << Search_accuracy << "  \t| " << std::endl;
		N++;
	}
	std::cout << "---------------------------------" << std::endl;


	//golden ratio
	double ak = a, bk = b, Golden_ratio = 1.618;
	double xk_left = ak + (1 - 1 / Golden_ratio) * (bk - ak);
	double xk_right = ak + (bk - ak) / Golden_ratio;
	std::cout << std::endl << "  Golden ratio  " << std::endl;
	std::cout << "_________________________________________" << std::endl;
	std::cout << "| Start\t|" << " End\t|" << "Length\t|" << " \t|" << "  \t|" << std::endl;
	std::cout << "| int\t|" << " int\t|" << " int \t|" << " f(ak)\t|" << " f(bk)\t|" << std::endl;
	std::cout << "| (ak)\t|" << "(bk) \t|" << " (l) \t|" << "  \t|" << "  \t|" << std::endl;
	std::cout << "-----------------------------------------" << std::endl;
	std::cout << "| " << ak << "\t|" << bk << " \t|" << bk - ak
		<< "\t|" << xk_left << "\t|" << xk_right << "\t|" << std::endl;
	while ((bk - ak) > 0.1)
	{
		double y1 = Function(xk_left);
		double y2 = Function(xk_right);

		if (y1 < y2)
		{
			bk = xk_right;
			xk_right = xk_left;
			xk_left = ak + bk - xk_right;
		}
		else
		{
			ak = xk_left;
			xk_left = xk_right;
			xk_right = ak + bk - xk_left;

		}
		std::cout << std::setprecision(3) << "| " << ak << "\t|" << bk << " \t|" << bk - ak
			<< "\t|" << xk_left << "\t|" << xk_right << "\t|" << std::endl;
	}
	std::cout << "-----------------------------------------" << std::endl;
	std::cout << "Minimum of function at " << std::setprecision(3) << (ak + bk) / 2
		<< " +- " << (bk - ak) / 2 << std::endl;
}