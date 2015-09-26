/* ZOJ Numerical Analysis Problem Set - 1002
 * Root of a Polynomial
 * Algorithm: use Newton's method to find the root of a polynomial
 */
 
double Fx(int n, double c[], double x){
    int i; 
	double p = c[n];
    for (i = n-1; i >= 0; i--)
        p = p*x + c[i];
    return p;
}

double dFx(int n, double c[], double x){
    int i;
	double p = n*c[n];
    for(i = n-1; i > 0; i--)
        p = p*x + i*c[i];
    return p;
}

double d2Fx(int n, double c[], double x){
    int i;
	double p = n*(n-1)*c[n];
    for(i = n-1; i > 1; i--)
        p = p*x + i*(i-1)*c[i];
    return p;
}

double Polynomial_Root(int n, double c[], double a, double b, double EPS){
    double t;
    if(a > b){
        t = a; 
        a = b; 
        b = t;
    }
    double eps = EPS * 1e-7;
    double minError = 999;
    int MAX_ITERATION = 1000;
	double x, root, f, df, d2f;
    int i, j;
    for(i = 0; i < 10; i++){
        x = a + (b-a)*i/10;
        j = 0;
        while(j <= MAX_ITERATION){
            f = Fx(n, c, x);
            df = dFx(n, c, x);
            d2f = d2Fx(n, c, x);
            j++;
            if(fabs(df*df-f*d2f) < eps) 
                break;
            x = x - (f*df)/(df*df - f*d2f);
            if(x<a || x>b) 
                break;
        }
        double error = fabs(Fx(n, c, x));
        if(a<=x && x<=b && error<minError) {
            root = x;
			minError = error;
        }
    }
    if(fabs(root) < eps)
		return fabs(root);
	else
		return root;
}