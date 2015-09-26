/* ZOJ Numerical Analysis Problem Set - 1007
 * Orthogonal Polynomials Approximation
 * Algorithm: use Gram-Schimidt orthogonalization method to calculate a group of orthogonal polynomials
 */
 
double CalculatePoly(double x, double c[], int k) {
    double result = 0;
    int i;
    for(i = k; i >= 0; i--) 
        result = result * x + c[i];
    return result;
}

void CalculatePhi1(int m, double x[], double w[], double phi[], double phi1[]) {
	int i;
    double p, q, B1, Phi0;
    p = q =  0.0;
    for(i = 0; i < m; i++) {
    	Phi0 = CalculatePoly(x[i], phi, 5);
        p += x[i] * w[i] * Phi0 * Phi0;
        q += w[i] * Phi0 * Phi0;
    }
    B1 = p / q;
    phi1[0] = - B1;
	phi1[1] = 1.0;
}

void GramSchmidt(int m, double x[], double w[], double phi1[], double phi2[], double phik[]) {
	int i;
    double p, q, r, s, Bi, Ci, Phi1, Phi2;
    p = q =  r = s = 0.0;
    for(i = 0; i < m; i++) {
    	Phi1 = CalculatePoly(x[i], phi1, 5);
    	Phi2 = CalculatePoly(x[i], phi2, 5);
        p += x[i] * w[i] * Phi1 * Phi1;
        q += w[i] * Phi1 * Phi1;
        r += x[i] * w[i] * Phi1 * Phi2;
        s += w[i] * Phi2 * Phi2;
    }
    Bi = p / q;
    Ci = r / s;
    phik[0] = - Bi*phi1[0] - Ci*phi2[0];
	for(i = 1; i <= MAX_n; i++)
    	phik[i] = phi1[i-1] - Bi*phi1[i] - Ci*phi2[i];
}

double CalculateAi(int m, double x[], double w[], double phi[], double(*f)(double t)) {
    int i;
    double p, q, Phi, Ai;
    p = q = 0.0;
    for(i = 0; i < m; i++) {
		Phi = CalculatePoly(x[i], phi, 5);
        p += w[i] * Phi * (*f)(x[i]);
        q += w[i] * Phi * Phi; 
    }
    Ai = p / q;
    return Ai;
}

int OPA( double (*f)(double t), int m, double x[], double w[],
double c[], double *eps ) {
    double phi[6][6];
    double Phi, value;
    double a[6];
    double Error = 0.0;
	int i, j, k;
	for(i = 0; i <= MAX_n; i++){
		a[i] = 0;
		c[i] = 0;
    	for(j = 0; j <= MAX_n; j++)
    		phi[i][j] = 0.0;
    }
    phi[0][0] = 1.0;
    a[0] = CalculateAi(m, x, w, phi[0], f);
    k = 0;
	for(i = 0; i < m; i++) {
		value = 0.0;
		for(j = 0; j <= k; j++) {
			Phi = CalculatePoly(x[i], phi[j], 5);
			value += a[j]*Phi;
		}
		Error += w[i] * ((*f)(x[i])-value) * ((*f)(x[i])-value);
	}
	if(fabs(Error) < *eps) {
		k = 1;
	}
	else {
		Error = 0.0;
    	CalculatePhi1(m, x, w, phi[0], phi[1]);
   		a[1] = CalculateAi(m, x, w, phi[1], f);
		k = 1;
		for(i = 0; i < m; i++) {
			value = 0.0;
			for(j = 0; j <= k; j++) {
				Phi = CalculatePoly(x[i], phi[j], 5);
				value += a[j]*Phi;
			}
			Error += w[i] * ((*f)(x[i])-value) * ((*f)(x[i])-value);
		}
		if(fabs(Error) < *eps) {
			k = 2;
		}
		else {
			k = 2;
			while((k<=MAX_n) && (fabs(Error) >= *eps)) {
				GramSchmidt(m, x, w, phi[k-1], phi[k-2], phi[k]);
				a[k] = CalculateAi(m, x, w, phi[k], f);
				Error = 0.0;
				for(i = 0; i < m; i++) {
					value = 0.0;
					for(j = 0; j <= k; j++) {
						Phi = CalculatePoly(x[i], phi[j], 5);
						value += a[j]*Phi;
					}
					Error += w[i] * ((*f)(x[i])-value) * ((*f)(x[i])-value);
				}
				k++;
			}
		}
	}
	for(i = 0; i <= MAX_n; i++)
		for(j = 0; j <= MAX_n; j++)
			c[i] += a[j]*phi[j][i];
	*eps = Error;
	return k-1;
}
