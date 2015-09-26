#include <stdio.h>
#include <math.h>

#define MAX_m 200
#define MAX_n 6

double f1(double x)
{
	return sin(x);
}

double f2(double x)
{
	return exp(x);
}

double f3(double x)
{
	return (1.0/x);
}

double f4(double x)
{
	return (0.5*cos(x)+0.3333*sin(x+x));
}

int OPA(double (*f)(double t), int m, double x[], double w[], double c[], double *eps);

void print_results( int n, double c[], double eps)
{	
	int i;

	printf("%d\n", n);
	for (i=0; i<=n; i++)
		printf("%12.4e ", c[i]);
	printf("\n");
	printf("error = %9.2e\n", eps);
	printf("\n");
}

int main()
{
	int m, i, n;
	double x[MAX_m], w[MAX_m], c[MAX_n+1], eps;

	m = 90;
	for (i=0; i<m; i++) {
		x[i] = 3.1415926535897932 * (double)(i+1) / 180.0;
		w[i] = 1.0;
	}
	eps = 0.001;
	n = OPA(f1, m, x, w, c, &eps);
	print_results(n, c, eps);
	eps = 1;
	n = OPA(f1, m, x, w, c, &eps);
	print_results(n, c, eps);
	eps = 0.0000000001;
	n = OPA(f1, m, x, w, c, &eps);
	print_results(n, c, eps);
	eps = 0.00000000000001;
	n = OPA(f1, m, x, w, c, &eps);
	print_results(n, c, eps);
	m = 200;
	for (i=0; i<m; i++) {
		x[i] = 0.01*(double)i;
		w[i] = 1.0;
	}
	eps = 0.001;
	n = OPA(f2, m, x, w, c, &eps);
	print_results(n, c, eps);
	m = 200;
	for (i=0; i<m; i++) {
		x[i] = 1.0+0.01*(double)i;
		w[i] = 1.0;
	}
	eps = 0.0001;
	n = OPA(f3, m, x, w, c, &eps);
	print_results(n, c, eps);
	m = 100;
	for (i=0; i<m; i++) {
		x[i] = 0.01*(double)i;
		w[i] = 1.0;
	}
	eps = 0.01;
	n = OPA(f4, m, x, w, c, &eps);
	print_results(n, c, eps);

	return 0;
}

