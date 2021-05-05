#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>

double czas(struct rusage *ru0, struct rusage *ru1){
	
	double utime = 0, stime = 0, ttime = 0;	

  	/* Obliczenie czasow. Aby mikrosekundy traktowac jako czesci sekund musimy je wymnozyc przez 10^-6*/
	utime = (double) ru1->ru_utime.tv_sec 
		+ 1.e-6 * (double) ru1->ru_utime.tv_usec 
		- ru0->ru_utime.tv_sec 
		- 1.e-6 * (double) ru0->ru_utime.tv_usec;
  	stime = (double) ru1->ru_stime.tv_sec
    		+ 1.e-6 * (double) ru1->ru_stime.tv_usec 
		- ru0->ru_stime.tv_sec
    		- 1.e-6 * (double) ru0->ru_stime.tv_usec;
	ttime = stime + utime;

	/*printf("user time: %3f\n", utime);
	printf("system time: %3f\n", stime);*/
	// printf("total time: %3f\n", ttime);

    return ttime;
} 

double LU_decomp(int n, gsl_matrix* m_main, gsl_vector_view b){
    puts("Metoda LU_decomp");
    struct rusage t0, t1;
    
    gsl_vector *x = gsl_vector_alloc (n);
    gsl_matrix *m = gsl_matrix_alloc (n,n);
    gsl_matrix_memcpy(m, m_main);

    int s;
    gsl_permutation * p = gsl_permutation_alloc (n);

    getrusage(RUSAGE_SELF, &t0);
    // dekompozycja i rozwiazanie
    gsl_linalg_LU_decomp (m, p, &s);
    gsl_linalg_LU_solve (m, p, &b.vector, x);

    getrusage(RUSAGE_SELF, &t1);

    // wyprintowanie wektora x
    // printf ("x = \n");
    // gsl_vector_fprintf (stdout, x, "%g");

    gsl_matrix_memcpy(m, m_main);

    // sprawdzenie czy A*x = b
    gsl_vector *result = gsl_vector_alloc(n);

    gsl_blas_dgemv(CblasNoTrans,
                  1.0, m, x,
                  0.0, result);

    // for(int i =0; i<n; i++){
    //     printf("Sprawdzam element: %d -> ",i);
    //     printf("%f = ", gsl_vector_get(&b.vector,i));
    //     printf("%f\n", gsl_vector_get(result,i));
    // }

    gsl_permutation_free (p);
    gsl_vector_free (x);
    gsl_vector_free(result);

    return czas(&t0, &t1);
}

double LU_invert(int n, gsl_matrix* m_main, gsl_vector_view b){
    puts("Metoda LU_invert");
    struct rusage t0, t1;
    int s;
    gsl_matrix *m = gsl_matrix_alloc (n,n);
    gsl_matrix_memcpy(m, m_main);
    gsl_permutation * p = gsl_permutation_alloc (n);

    getrusage(RUSAGE_SELF, &t0);

    // dekompozycja i odwrocenie macierzy
    gsl_linalg_LU_decomp(m, p, &s);
    gsl_matrix *inv = gsl_matrix_alloc(n,n);
    gsl_linalg_LU_invert(m, p, inv);

    gsl_vector *x = gsl_vector_alloc (n);

    // mnozenie macierzy inv i wektora b -> x bo x = A^{-1}*b
    gsl_blas_dgemv(CblasNoTrans,
                  1.0, inv, &b.vector,
                  0.0, x);

    // wyprintowanie wektora x
    // printf ("x = \n");
    // gsl_vector_fprintf (stdout, x, "%g");

    getrusage(RUSAGE_SELF, &t1);

    gsl_permutation_free (p);

    // utworzenie macierzy jednostkowej
    gsl_matrix *unit = gsl_matrix_alloc(n,n);
    gsl_matrix_set_identity(unit); 

    // sprawdzenie dwoch rzeczy: czy AA-1=I i A-1A=I
    gsl_matrix *result1 = gsl_matrix_alloc(n,n);
    gsl_matrix *result2 = gsl_matrix_alloc(n,n);

    gsl_matrix_memcpy(m, m_main);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                  1.0, m, inv,
                  0.0, result1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                  1.0, inv, m,
                  0.0, result2);


    if (gsl_matrix_equal(result1, unit) && gsl_matrix_equal(result2, unit)){
        puts("LU_inv -> Success!\n");
    }

    gsl_matrix_free(inv);
    gsl_matrix_free(result1);
    gsl_matrix_free(result2);

    return czas(&t0, &t1);

}

double QR_decomp(int n, gsl_matrix* m_main, gsl_vector_view b){
    puts("Metoda QR_decomp");
    struct rusage t0, t1;

    gsl_matrix *m = gsl_matrix_alloc (n,n);
    gsl_matrix_memcpy(m, m_main);
    gsl_vector *tau = gsl_vector_alloc (n);
    gsl_vector *x1 = gsl_vector_alloc (n);

    getrusage(RUSAGE_SELF, &t0);

    // dekompozycja i rozwiazanie
    gsl_linalg_QR_decomp(m, tau);
    gsl_linalg_QR_solve (m, tau, &b.vector, x1); 

    getrusage(RUSAGE_SELF, &t1);

    // wyprintowanie wektora x
    // printf ("x = \n");
    // gsl_vector_fprintf (stdout, x1, "%g");

    gsl_matrix_memcpy(m, m_main);

    // sprawdzenie czy A*x = b
    gsl_vector *result = gsl_vector_alloc(n);
    gsl_blas_dgemv(CblasNoTrans,
                  1.0, m, x1,
                  0.0, result);

    // for(int i =0; i<n; i++){
    //     printf("Sprawdzam element: %d -> ",i);
    //     printf("%f = ", gsl_vector_get(&b.vector,i));
    //     printf("%f\n", gsl_vector_get(result,i));
    // }

    gsl_vector_free(tau);
    gsl_vector_free (x1);
    gsl_vector_free(result);

    return czas(&t0, &t1);
}


int main (int argc, char const *argv[])
{

    int n = atoi(argv[1]);

    double a_data[n*n];
    double b_data[n];

    for(int i=0; i<=n*n; i++){
        int r = rand()%100;
        a_data[i] = r;
    }

    for(int i=0; i<=n; i++){
        int r = rand()%100;
        b_data[i] = r;
    }

    gsl_matrix_view m = gsl_matrix_view_array (a_data, n, n);
    gsl_vector_view b = gsl_vector_view_array (b_data, n);

    puts("Macierz A:");
    gsl_matrix_fprintf (stdout, &m.matrix, "%g");
    puts("\n");
    puts("Wektor B:");
    gsl_vector_fprintf (stdout, &b.vector, "%g");
    puts("\n");
    

    printf("Czas LU_decomp dla rozmiaru = %d: %f \n", n, LU_decomp(n,&m.matrix,b));
    printf("Czas LU_invert dla rozmiaru = %d: %f \n", n, LU_invert(n,&m.matrix,b));
    printf("Czas QR_decomp dla rozmiaru = %d: %f \n", n, QR_decomp(n,&m.matrix,b));



    // gsl_vector_fprintf (stdout, &b.vector, "b1: %g");
    // int time_LU_decomp = LU_decomp(n,&m.matrix,b);
    // int time_LU_invert = LU_invert(n,&m.matrix,b);
    // int time_QR_DECOMP = QR_decomp(n,&m.matrix,b);


    
    return 0;

}