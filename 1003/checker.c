#include <stdio.h>

#define Max_size 10000 /* max number of dishes */

void Price( int n, double p[] );

int main()
{
	int n, i;
	double p[Max_size];

	while (scanf("%d", &n)!=EOF) {
		for (i=0; i<n; i++) 
			scanf("%lf", &p[i]);
		Price(n, p);
		for (i=0; i<n; i++)
			printf("%.2f ", p[i]);
		printf("\n");
	}

	return 0;
}

