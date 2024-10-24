//Programa que resuelve Ecuación de Laplace (ecuaciones diferenciales parciales)
// Método de diferencias finitas

#include<stdio.h>
#include<math.h>

int main () {

double p[200][200], h=0.01, pp, dif, max, tol=1.0e-6;
int i,j,k, n=100;

for( i=1 ; i<=n-1 ; i++ ){

    for( j=1 ; j<=n-1 ; j++ ){
		
        p[i][j]=50;  //guess inicial

    }


}

for( i=0 ; i<=n ; i++ ){  //condiciones de frontera
	
	p[0][i]=0;
	p[n][i]=100;
	p[i][0]=0;
	p[i][n]=100;
	
}

    do{
		max=0;
        for( i=1 ; i<=n-1 ; i++ ){


            for( j=1 ; j<=n-1 ; j++ ){
				
				pp=p[i][j];
                p[i][j]=(1.0/4.0)*(p[i-1][j]+p[i+1][j]+p[i][j-1]+p[i][j+1]);
				
				dif=fabs(p[i][j]-pp);
				
				if(dif>max)max=dif;
				
            }

        }

        k++;
    }
    while(max>tol);

printf("El número de iteraciones es: %d\n", k);

    FILE*arch=fopen("rel.txt", "w");


    for( i=0 ; i<=n ; i++ ){

        for(j=0;j<=n;j++){
            fprintf(arch, "%lf %lf %lf\n", i*h , j*h , p[i][j]);     
            }
            fprintf(arch, "\n");
    }

    fclose(arch);



    return 0;
}
