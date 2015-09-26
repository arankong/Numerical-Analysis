#include <stdio.h>
#include <math.h>

double f0( double x, double l, double t )
{
	return sqrt(1.0+l*l*t*t*cos(t*x)*cos(t*x));
}

double Integral(double a, double b, double (*f)(double x, double y, double z), double eps, 
				double l, double t);

int main()
{
	double a=0.0, b, eps=0.005, l, t;

	while (scanf("%lf %lf %lf", &l, &b, &t) != EOF)
		printf("%.2f\n", Integral(a, b, f0, eps, l, t));

	return 0;
}

