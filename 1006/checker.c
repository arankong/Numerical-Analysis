#include <stdio.h>

#define MAX_N 20

void Cubic_Spline(int n, double x[], double f[], int Type, double s0, double sn, 
				  double a[], double b[], double c[], double d[]);

double S( double t, double Fmax, 
		  int n, double x[], double a[], double b[], double c[], double d[] );

int main()
{
	int n, Type, m, i;
	double x[MAX_N], f[MAX_N], a[MAX_N], b[MAX_N], c[MAX_N], d[MAX_N];
	double s0, sn, Fmax, t0, tm, h, t;

	while (scanf("%d", &n) != EOF) {
		for (i=0; i<=n; i++) 
			scanf("%lf", &x[i]);
		for (i=0; i<=n; i++) 
			scanf("%lf", &f[i]);
		scanf("%d %lf %lf %lf", &Type, &s0, &sn, &Fmax);
		/* Fmax is the default value that S(x) will assume if x is out of the range [x[0], x[n]] */

		Cubic_Spline(n, x, f, Type, s0, sn, a, b, c, d);
		for (i=1; i<=n; i++)
			printf("%12.8e %12.8e %12.8e %12.8e \n", a[i], b[i], c[i], d[i]);

		scanf("%lf %lf %d", &t0, &tm, &m);
		h = (tm-t0)/(double)m;
		for (i=0; i<=m; i++) {
			t = t0+h*(double)i;
			printf("f(%12.8e) = %12.8e\n", t, S(t, Fmax, n, x, a, b, c, d));
		}
		printf("\n");
	}

	return 0;
}

