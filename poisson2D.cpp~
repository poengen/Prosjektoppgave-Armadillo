//C++ code for solving poisson 2D

#include <iostream>
#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;

#define PI 3.14159

int main(int argc, char** argv){
	
// THE GRID
	int n;
	double h;
	n = 100;
	h = 1/double(n+1);

	vec x = zeros<vec>(n+2);
	for (int i = 0; i<n+2; i++){
		x(i) = i*h;
	}
	//cout << h << endl;
	
// LOAD VECTOR
	mat fm = zeros<mat>(n,n);
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			fm(i,j) = 2*PI*PI*sin(PI*x(i+1))*sin(PI*x(j+1));
		}
	}
	vec f = vectorise(fm);
	//cout << f << endl;

// MATRICES
	mat A = zeros<mat>(n*n,n*n);
	for (int i = 0; i < n*n; i++){ A(i,i) = 4; }
	for (int i = 0; i < n*n-n; i++) { A(i,n+i) = -1; A(n+i,i) = -1; }
	for (int i = 0; i < n*n-1; i++) { A(i,i+1) = -1; A(i+1,i) = -1; }
	for (int i = 0; i < n-1; i++) { A((i+1)*n-1,(i+1)*n) = 0; A((i+1)*n,(i+1)*n-1) = 0; }
	A = 1/(h*h)*A;
	//cout << A << endl;
	
// SOLUTION VECTOR
	mat um = zeros<vec>(n*n);
	vec uv = vectorise(um);

// SOLVE LINEAR SYSTEM
	uv = solve(A,f);
	mat uh = reshape(uv,n,n);
	mat u = zeros<mat>(n+2,n+2);
	u(span(1,n),span(1,n)) = uh;

	//cout << u << endl;

// WRITE SOLUTION TO FILE
	u.save("u.dat", raw_ascii);
	
}
