/*
 * Calculate the mode mixing matrix for 2D power spectrum measurements
 * Input: C_l power spectrum of a mask, e.g. output of anafast(mask)
 * Output: Mll' which describes \tilde C_l = anafast(mask*map) = M_ll' C_l
 * lmax, the dimension of Mll', is hard coded in main()
 * The input power spectrum should start with the l=0 element, i.e. mean(mask^2)
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#define WMAX 10000
#define LFT_LMAX 10000

double lft[LFT_LMAX];

/*
 * Fill the table of log factorials
 */
void fill_lft()
{
	int i;
	for(i=1;i<LFT_LMAX;i++) {
		lft[i] = log(i);
	}
	lft[0] = lft[1];
	for(i=1;i<LFT_LMAX;i++) {
		lft[i] = lft[i] + lft[i-1];
	}
}

/*
 * Calculate a Wigner 3j symbol using the explicit formula from the appendix of the MASTER paper
 * For numerical stability, the rational function of factorials is computed as the
 * exponetial of the sum of logarithms of factorials.
 * Log factorials must be precomputed in the global array lft[]
 */
double wigner3j(int l1,int l2,int l3)
{
	double sign,term1,term2,r;
	int L = l1 + l2 + l3;
	assert(L+1 < LFT_LMAX);
	if(abs(l1 - l2) > l3)
		return 0.;
	if(l3 > l1 + l2)
		return 0.;
	if(L%2 == 1)
		return 0.;
	if((L/2 % 2) == 0)
		sign = 1.;
	else
		sign = -1.;

	term1 = 0.5*(lft[L-2*l1] + lft[L-2*l2] + lft[L-2*l3] - lft[L+1]);
	term2 = lft[L/2] - lft[L/2 - l1] - lft[L/2 - l2] - lft[L/2 - l3];

	r = sign * exp(term1 + term2);
	return r;
}

/*
 * Fill m with the mode mixing matrix Mll', with size lmax * lmax
 * w is the input power spectrum which must have size >= lmax
 */
void modemix(double *m, double *w, int lmax)
{
	int l1,l2,l3;
	double w3j,s,t;
	for(l1=0;l1<lmax;l1++) {
		printf("%d ",l1);
		fflush(stdout);
		for(l2=0;l2<lmax;l2++) {
			s = (2*l2 + 1) / (4. * M_PI);
			t = 0.;
			int l3a = abs(l1-l2);
			int l3b = l1 + l2;
			for(l3=l3a;l3<=l3b;l3++) {
				w3j = wigner3j(l1,l2,l3);
				t += (2*l3+1) * w[l3] * w3j*w3j;
			}
			m[l1*lmax + l2] = s * t;
		}
	}
	printf("\n\n");
}

void output_matrix(FILE *fo, double *m, int lmax)
{
	int l1,l2;
	for(l1=0;l1<lmax;l1++) {
		for(l2=0;l2<lmax;l2++) {
			fprintf(fo,"%e ",m[l1*lmax+l2]);
		}
		fprintf(fo,"\n");
	}
	fprintf(fo,"\n");
}

int main(int argc, char **argv)
{
	assert(argc == 4);
	fill_lft();
	int i;
	double w[WMAX];
	char *infn = argv[1];
	char *outfn = argv[2];
	int lmax = atoi(argv[3]);

	assert(LFT_LMAX >= 3*lmax+1);

	FILE *f = fopen(infn,"r");
	for(i=0;i<WMAX;i++) {
		int res = fscanf(f,"%lf",&w[i]);
		if(res != 1)
			break;
	}
	fclose(f);
	printf("input power spectrum size = %d\n", i);

	assert(i >= lmax);
	printf("lmax = %d\n",lmax);

	double *m = (double *)calloc(lmax*lmax,sizeof(double));

	modemix(m,w,lmax);
	
	FILE *fo = fopen(outfn,"w");
	output_matrix(fo,m,lmax);
	fclose(fo);

	return 0;
}
