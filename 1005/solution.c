/* ZOJ Numerical Analysis Problem Set - 1005
 * Approximating Eigenvalues
 * Algorithm: use Doolittle Method to approximate the eigenvalue of a matrix
 */
 
double fabs(double x) {
	if(x < 0)
		return -x;
	else
		return x;
}

int Doolittle(int n, double a[][MAX_SIZE]) {
	int i, j, k;
	double sum;
	if(a[0][0] == 0)
		return 0;
	for(i = 1; i < n; i++)
		a[i][0] /= a[0][0];
	for(i = 1; i < n; i++) {
		if(a[i][i] == 0)
			return 0;
		for(j = i; j < n; j++) {
			sum = 0;
			for(k = 0; k < i; k++)
				sum += a[i][k]*a[k][j];
			a[i][j] -= sum;
		}
		for(j = i+1; j < n; j++) {
			sum = 0;
			for(k = 0; k < i; k++)
				sum += a[j][k]*a[k][i];
			a[j][i] = (a[j][i]-sum) / a[i][i];
		}
	}
	return 1;
}

int EigenV(int n, double b[][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN) {
	int i, j, k, index;
	double miu, ERR, max, sum;
	double a[MAX_SIZE][MAX_SIZE], u[MAX_SIZE];
	double lamda = *lambda;
	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			a[i][j] = b[i][j];
	for(i = 0; i < n; i++)
		a[i][i] -= lamda;
	if(!Doolittle(n, a))
		return -1;
	for(i = 0; i < n; i++)
		if(a[i][i] == 0)
			return -1;
	k = 1;
	max = 0;
	for(i = 0; i < n; i++)
		if(fabs(v[i]) > fabs(max)) {
			max = v[i];
			index = i;
		}
	for(i = 0; i < n; i++)
		v[i] /= max;
	while(k <= MAXN) {
		for(i = 0; i < n; i++)
			u[i] = v[i];
		for(i = 0; i < n; i++){
			sum = 0;
			for(j = 0; j < i; j++)
				sum += a[i][j]*v[j];
			v[i] -= sum;
		}
		for(i = n-1; i >= 0; i--) {
			sum = 0.0;
			for(j = i+1; j < n; j++)
				sum += a[i][j]*v[j];
			v[i] = (v[i] - sum)/a[i][i];
		}
		miu = v[index];
		max = 0;
		for(i = 0; i < n; i++)
			if(fabs(v[i]) > fabs(max)) {
				max = v[i];
				index = i;
			}
		for(i = 0; i < n; i++)
			v[i] /= max;
		ERR = 0;
		for(i = 0; i < n; i++)
			if(fabs(u[i]-v[i]) > ERR)
				ERR = fabs(u[i]-v[i]);
		if(ERR < TOL) {
			*lambda = 1.0/miu + lamda;
			return 1;
		}
		k++;
	}
	return 0;
}

