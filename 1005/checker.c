#include <stdio.h>

#define MAX_SIZE 100

int EigenV(int n, double b[][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN);

int main()
{
	int n, MAXN, m, i, j, k;
	double a[MAX_SIZE][MAX_SIZE], v[MAX_SIZE];
	double lambda, TOL;

	while (scanf("%d", &n) != EOF) {
		for (i=0; i<n; i++) 
			for (j=0; j<n; j++) 
				scanf("%lf", &a[i][j]);
		scanf("%lf %d", &TOL, &MAXN);
		scanf("%d", &m);
		for (i=0; i<m; i++) {
			scanf("%lf", &lambda);
			for (j=0; j<n; j++)
				scanf("%lf", &v[j]);
			switch (EigenV(n, a, &lambda, v, TOL, MAXN)) {
				case -1: 
					printf("%12.8f is an eigenvalue.\n", lambda );
					break;
				case 0:
					printf("Maximum number of iterations exceeded.\n");
					break;
				case 1:
					printf("%12.8f\n", lambda );
					for (k=0; k<n; k++)
						printf("%12.8f ", v[k]);
					printf("\n");
					break;
			}
		}
		printf("\n");
	}

	return 0;
}

