/* ZOJ Numerical Analysis Problem Set - 1003
 * There is No Free Lunch
 * Algorithm: use LU Decomposition(Crout) to solve special linear system of equations 
 */
 
void Price(int n, double p[]) {
	double a[Max_size], b[Max_size], c[Max_size];
	int i;
	double temp;	
	for(i = 0; i < n; i++) {
		a[i] = 0.50;
		b[i] = 2.0;
		c[i] = 0.50;
	}
	a[0] /= b[0];
	c[0] /= b[0];
	p[0] /= b[0];
	for(i = 1; i < n-1; i++) {
		temp = b[i] - a[i]*c[i-1];
		c[i] /= temp;
		p[i] = (p[i] - a[i]*p[i-1])/temp;
		a[i] = -a[i]*a[i-1]/temp;
	}
	a[n-2] = -a[n-2] - c[n-2];
	for(i = n-3; i >= 0; i--) {
		a[i] = -a[i] - c[i]*a[i+1];
		p[i] -= c[i]*p[i+1];
	}
	p[n-1] -= (c[n-1]*p[0] + a[n-1]*p[n-2]);
	p[n-1] /= (c[n-1]*a[0] + a[n-1]*a[n-2] + b[n-1]);
	for(i = n-2; i >= 0; i--)
		p[i] += a[i]*p[n-1];
}