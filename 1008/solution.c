/* ZOJ Numerical Analysis Problem Set - 1008
 * Shape Roof
 * Algorithm: use Romberg Method to calculate numerical integral for given funtions 
 */

double CalculateIntegral(double a, double b, double (*f)(double x, double y, double z), double eps, double l, double t) {
#define MAX_ROW  1000
	double R[2][MAX_ROW];
	int i, j;
	double n, k;
	int m = 1;
	double h = b - a;
	double temp;
	int min = (int)(log(h*10.0) / log(2.0));
	for (i = 0; i < MAX_ROW; i++) {
		R[0][i] = 0.0;
		R[1][i] = 0.0;
	}
	
	R[0][0] = 0.50 * h * ((*f)(a, l, t) + (*f)(b, l, t)); 
	for(i = 2; i < MAX_ROW; i++) {
		temp = 0.0;
		for(k = 1 ; k <= m; k++)
			temp += (*f)(a + h*(k-0.50), l, t);
		R[1][0] = 0.50 * (R[0][0] + h*temp);
		n = 4.0;
		for(j = 1; j < i; j++) {
			R[1][j] = R[1][j-1] + (R[1][j-1] - R[0][j-1]) / (n - 1.0);
			n *= 4.0;
		}
		if(((fabs(R[1][i-1]-R[0][i-2])<eps) && (i>min)))
			return R[1][i-1];
		h *= 0.50;
		m *= 2;
		for(j = 0; j < i; j++)
			R[0][j] = R[1][j];
	}
	return R[1][MAX_ROW-1];
}

double Integral(double a, double b, double (*f)(double x, double y, double z), double eps, double l, double t) {
#define PI 3.1415926535897932
	int n;
	double T, remain, result;
	T = 2.0 * PI / t;
	n = (int)(floor)(b * t * 0.50 / PI);
	remain = b - (double)n * T;
	if(n)
		result = CalculateIntegral(a, T, f0, eps/(double)n, l, t);
	result = (double)n * result + CalculateIntegral(0.0, remain, f0, eps, l, t );
	result /= 100.0;
	return result;
}

