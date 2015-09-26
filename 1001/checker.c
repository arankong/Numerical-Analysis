#include <stdio.h>

void Series_Sum( double sum[] );

int main()
{
	int i;
	double x, sum[3001];
	
	Series_Sum( sum );

	x = 0.0;
	for (i=0; i<3001; i++)
		printf("%6.2lf %16.12lf\n", x + (double)i * 0.10, sum[i]);

	return 0;
}

