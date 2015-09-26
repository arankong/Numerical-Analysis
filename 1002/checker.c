#include <stdio.h>
#include <math.h>

#define ZERO 0.000000001 /* X is considered to be 0 if |X|<ZERO */
#define MAXN 11 	/* Max Polynomial Degree + 1 */

double Polynomial_Root(int n, double c[], double a, double b, double EPS);

int main()
{
	int n;
	double c[MAXN], a, b;
    double EPS = 0.00005;
	int i;

	while (scanf("%d", &n)!= EOF){
		for (i=n; i>=0; i--) 
			scanf("%lf", &c[i]);
		scanf("%lf %lf", &a, &b);
		printf("%.4lf\n", Polynomial_Root(n, c, a, b, EPS));
	}

	return 0;
}

