/* ZOJ Numerical Analysis Problem Set - 1004
 * Compare Methods of Jacobi with Gauss-Seidel
 * Algorithm: implement both Jacobi Method and Gauss-Seidel Method
 */
 
bool IsZero(double x) {
	if(fabs(x) < ZERO)
		return true;
	else
		return false;
}

bool SetAii(int n, double a[][MAX_SIZE], double b[]) {
	int i, j, index;
	double max, temp;
	for(i = 0; i < n; i++) {
	max = 0;
	for(j = i; j < n; j++)
		if(fabs(a[j][i]) > max) {
			max = a[j][i];
			index = j;
		}
		if(IsZero(max)) {
			for(j = i-1; j >= 0; j--)
				if(fabs(a[j][i]) > max) {
					max = a[j][i];
					index = j;
				}
			if(IsZero(max))
				return false;
			else {
				for(j = 0; j < n; j++)
					a[i][j] += a[index][j];
				b[i] += b[index];
			}
		}
		else {
			for(j = 0; j < n; j++) {
				temp = a[i][j];
				a[i][j] = a[index][j];
				a[index][j] = temp;
			}
			temp = b[i];
			b[i] = b[index];
			b[index] = temp;
		}
	}
	return true;
}

int Jacobi( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN ) {
	int i, j, k;
	double temp, delta, y[MAX_SIZE];
	if(!SetAii(n, a, b))
		return -1;
	k = 1;
	while(k <= MAXN) {
		delta = 0;
		for(i = 0; i < n; i++) {
			temp = 0;
			for(j = 0; j < n; j++)
				if(j != i)
					temp += a[i][j]*x[j];
			y[i] = (-temp + b[i]) / a[i][i];
			if(fabs(y[i]) > bound)
				return -2;
			if(fabs(y[i] - x[i]) > delta)
				delta = fabs(y[i] - x[i]);
		}
		for(i = 0; i < n; i++)
			x[i] = y[i];
		if(delta < TOL)
			return k;
		k++;
	}
	return 0;
}

int Gauss_Seidel( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN ) {
	int i, j, k, index;
	double temp, delta;
	k = 1;
	if(!SetAii(n, a, b))
		return -1;
	while(k <= MAXN) {
		delta = 0;
		for(i = 0; i < n; i++) {
			temp = 0;
			for(j = 0; j < n; j++)
				if(j != i)
					temp += a[i][j]*x[j];
			temp = (-temp + b[i]) / a[i][i];
			if(fabs(temp) > bound)
				return -2;
			if(fabs(temp - x[i]) > delta)
				delta = fabs(temp - x[i]);
			x[i] = temp;
		}
		if(delta < TOL)
			return k;
		k++;
	}
	return 0;	
} 
