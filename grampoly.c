#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

// program to compute savitksy-golay gram polynomial coefficients based
// on (s=0) or its sth derivative evaluated at i, order k, over 2m+1
// points, centered with offset t

int s=0;	// 0=smoothing, 1=first_derivative ...
int m=3;	// width = 2*m+1
int n=4;	// 2=quadratic fit
int t=0;	// how to center the estimate 

// General Least-Squares Smoothing and Differentiation by the
// Convolution (Savitzky-Golay) Method, Peter A. Gorry

char *progname;
void license(void);

double gp(int i, int m, int k, int s) {

    if (k>0) {
	return (4.0*k-2.0)/(k*(2.0*m-k+1.0))*(i*gp(i,m,k-1,s) +
	  s*gp(i,m,k-1,s-1))-((k-1.0)*(2.0*m+k))/(k*(2.0*m-k+1.0))*gp(i,m,k-2,s);    
    } else {
	if (k==0 && s==0) {
	    return 1.0;
	} else {
	    return 0.0;
	}
    }
}

// generalized factorial
// (a)(a-1)...(a-b+1)

double genfact(int a, int b) {
    int j;
    double gf;
    gf=1.0;
    for (j=(a-b)+1; j<=a; j++) {
	gf*=j;
    }
    return(gf);
}

// calculate the weight of the ith dat point for the tth
// lease-square point of the sth derivative, over 2m+1 points
// order n
double weight(int i, int t, int m, int n, int s) {
    int k;
    double w = 0;
    for (k=0; k<=n; k++) {
	w += (2*k+1)*(genfact(2*m,k)/genfact(2*m+k+1,k+1))*gp(i,m,k,0)*gp(t,m,k,s);
    }
    return w;
}

void usage(void)
{
    fprintf(stderr, "usage: %s [options]\n", progname);
    fprintf(stderr, "    [-l] print license\n", s);
    fprintf(stderr, "    [-s <smoothing>] 0=smoothing, 1=first_derivative (default %d)\n", s);
    fprintf(stderr, "    [-m <points>] width = 2*m+1 (default %d)\n", m);
    fprintf(stderr, "    [-n <order>] 2=quadratic fit (default %d)\n", n);
    fprintf(stderr, "    [-t <offset>] how to offset the estimate in window (default %d)\n", t);
}

int main(int argc, char **argv) {
    extern int optind;          /* argv index of next option */
    extern int opterr;
    extern char *optarg;
    int optval;
    opterr = 0;                 /* disables getopt's error msg's */
    // int errflag = 0;
    int i;

    progname = argv[0];

    while ((optval = getopt(argc, argv, "ls:m:n:t:")) != EOF)
        switch (optval) {
        case 'l':
	    license(); 
	    exit(1);
            break;
        case 's':
            s = atoi(optarg);
            break;
        case 'm':
            m = atoi(optarg);
            break;
        case 'n':
            n = atoi(optarg);
            break;
        case 't':
            t = atoi(optarg);
            break;
        case '?':
	default:
	    usage();
	    exit(1);
	    break;
    }

    printf("# s=%d m=%d n=%d t=%d\n", s, m, n, t);

    for (i=0; i<2*m+1; i++) {
	printf("%d %g\n",i, weight(i-m, t, m, n, s));
    }
}


void license(void) {
    printf("    grampoly(1), computes Savitsky-Golay FIR filter coefficients\n");
    printf("    Copyright (C) 2021-2024  by Richard C. Walker.\n");
    printf("\n");
    printf("    Author's email: walker AT omnisterra DOT com\n");
    printf("    Homepage: http://www.omnisterra.com/walker/linux/piglet/intro.htm\n");
    printf("\n");
    printf(" Permission is hereby granted, free of charge, to any person obtaining a copy\n");
    printf(" of this software and associated documentation files (the \"Software\"), to deal\n");
    printf(" in the Software without restriction, including without limitation the rights\n");
    printf(" to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n");
    printf(" copies of the Software, and to permit persons to whom the Software is\n");
    printf(" furnished to do so, subject to the following conditions:\n");
    printf(" \n");
    printf(" The above copyright notice and this permission notice shall be included in all\n");
    printf(" copies or substantial portions of the Software.\n");
    printf(" \n");
    printf(" THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n");
    printf(" IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n");
    printf(" FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n");
    printf(" AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n");
    printf(" LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n");
    printf(" OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n");
    printf(" SOFTWARE.\n");
    printf("\n");
}
