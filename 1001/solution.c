/* ZOJ Numerical Analysis Problem Set - 1001
 * Numerical Summation of a Series
 * Algorithm: to accelerate the convergence speed by increasing the degree of the polynomials in the denominator
 */
 
void Series_Sum( double sum[] ){
	int i;
	double k, x;
	x = 0.00;
	const double g3 = (double)1.0 /18;
	const double h4 = (double)1.0 /96;  	
	for(i = 0; i < 3001; i++ ){
		sum[i] = 0.000;
		for(k = 1; k < 10000; k++)
			sum[i] += 1.0 / (k * (k + 1) * (k + 2) * (k + 3) * (k + 4) * (k + x));
			sum[i] = (1 - x) * ((2 - x) * ((3 - x) * ((4 - x) * sum[i] + h4) + g3) + 0.25) + 1.0;
		x += 0.10;
	}
}