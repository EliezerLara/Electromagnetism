//Programa que resuelve Ecuación de Laplace (ecuaciones diferenciales parciales)
// Método de diferencias finitas

#include<stdio.h>
#include<math.h>

int main () {

float p[100][100], h=0.01;
int i,j,k;

for( i=1 ; i<=98 ; i++ ){

    for( j=1 ; j<=98 ; j++ ){

        p[i][j]=50;

    }

}

    do{
        for( i=1 ; i<=98 ; i++ ){


            for( j=1 ; j<=98 ; j++ ){

                p[i][j]=(1.0/4.0)*(p[i-1][j]+p[i+1][j]+p[i][j-1]+p[i][j+1]);

            }

        }

        k++;
    }
    while();



    FILE*arch=fopen("rel.txt", "w");


    for( i=1 ; i<=98 ; i++ ){

        for(j=0;j<=65;j++){
            fprintf(arch, "%lf %lf %lf\n", i*h , j*h , p[i][j]);     
            }
    }

    fclose(arch);



    return 0;
}