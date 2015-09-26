/* ZOJ Numerical Analysis Problem Set - 1006
 * Cubic Spline
 * Algorithm: use LD Decomposition(Crout) to help implement Cubic Spline
 */
 
double abs(double x) {
	if(x < 0)
		return -x;
	else
		return x;
}

void Crout(int n, double a[], double b[], double c[], double x[]) {
	double temp;
	double eps = 1e-7;
	int i;
	if(abs(b[0]) < eps) 
		return;
	c[0] /= b[0];
	x[0] /= b[0];
	for(i = 1; i < n-1; i++) {
		if (abs(temp = b[i]-a[i]*c[i-1]) < eps)  
			return;
		c[i] /= temp;
		x[i] = (x[i] - a[i]*x[i-1])/temp;
	}
	if(abs(temp = b[n-1] - a[n-1]*c[n-2]) < eps)  
		return;
	x[n-1] = (x[n-1]- a[n-1]*x[n-2])/temp;
	for(i = n-2; i >= 0; i--)
		x[i] -= c[i]*x[i+1];
	return;
}

void Cubic_Spline(int n, double x[], double f[], int Type, double s0, double sn, double a[], double b[], double c[], double d[]) {
	double h1, h2, M[MAX_N];
	int i;
	for(i = 1; i < n; i++) {
		b[i] = 2.0;
		h1 = x[i] - x[i-1];
		h2 = x[i+1] - x[i];
		a[i] = h1 / (h1+h2);
		c[i] = 1.0 - a[i];
		M[i] = 6.0 / (h1+h2) * ((f[i+1]-f[i])/h2 - (f[i]-f[i-1])/h1);
	}
	b[0] = b[n] = 2.0;
	switch(Type) {
	case 1:
		a[n] = c[0] = 1.0;
		M[0] = 6.0 * ((f[1]-f[0]) / (x[1]-x[0])-s0) / (x[1]-x[0]);
		M[n] = 6.0 * (sn-(f[n]-f[n-1]) / (x[n]-x[n-1])) / (x[n]-x[n-1]);
		break;
	case 2:
		a[n] = c[0] = 0.0;
		M[0] = s0 + s0;
		M[n] = sn + sn;
		break;
	}
	Crout(n+1, a, b, c, M);
	for(i = 1; i <= n; i++) {
		h1 = x[i] - x[i-1];
		a[i] = f[i-1];
		b[i] = (f[i]-f[i-1])/h1 - (M[i-1]*2.0+M[i])*h1/6.0;
		c[i] = M[i-1] * 0.50;
		d[i] = (M[i] - M[i-1]) / (6.0*h1);
	}
}

double S( double t, double Fmax, 
		  int n, double x[], double a[], double b[], double c[], double d[] )
{
	int i;
	double result;
	if((t < x[0]) || (t > x[n]))
		result = Fmax;
	else {
		i = 1;
		while(t > x[i])
			i++;
		t -= x[i-1];
		result = t*d[i] + c[i];
		result = result*t + b[i];
		result = result*t + a[i];
	}
	return result;
}

