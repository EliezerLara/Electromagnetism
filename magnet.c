//Programa que resuelve Ecuación de Poisson para caso magnetostático de un conductor cerca de una placa metálica
// Método de diferencias finitas

#include<stdio.h>
#include<math.h>

int main () {

double p[201][201], h=0.02 , pp , dif , max , tol=1.0e-6 , w=1.95;
double Q[201][201] , D=19500.0 , c[201][201], bx[201][201], by[201][201];
int i,j,k=0, n=200;


for( i=1 ; i<=n-1 ; i++ ){

    for( j=1 ; j<=n-1 ; j++ ){
		
        p[i][j]=1.0;  //guess inicial de la función p y de la densidad de carga
        Q[i][j]=0.0;
        c[i][j]=0.0;
    }


}

for(i=110; i<=130 ; i++){
    for( j=0 ; j<=100 ; j++){  //definimos la densidad de carga del conductor

        Q[i][j]=D;
    }
}

for(i=170;i<=172;i++){

    for(j=0;j<=170;j++){

        c[i][j]=100.0;          //definimos la permeabilidad de la placa de metal
        p[i][j]=0;
    }
}

for( i=0 ; i<=n ; i++ ){  //condiciones de frontera de dirichlet
	
	p[0][i]=0;  //borde izquierdo
	p[n][i]=0;  //borde derecho
	p[i][n]=0;  //borde superior
	
}

for(i=0; i<=n ; i++){  //condicion de frontera de Neumann en el borde inferior

    for(j=0;j<=n;j++){
    
        if(i<170){
            p[i][0]=(1.0/4.0)*( ((2.0/101.0)*p[i+1][j])+p[i][j+1]+(200.0/101.0)*p[i-1][j]+p[i][j-1] );
        }
        else if(i>172){
            p[i][0]=(1.0/4.0)*( ((200.0/101.0)*p[i+1][j])+p[i][j+1]+(2.0/101.0)*p[i-1][j]+p[i][j-1] );
        }
    }
}

    do{
		max=0;
        for( i=1 ; i<=n-1 ; i++ ){


            for( j=1 ; j<=n-1 ; j++ ){

            //Método de diferencias finitas

				pp=p[i][j];
                p[i][j]=(1.0/4.0)*(p[i-1][j]+p[i+1][j]+p[i][j-1]+p[i][j+1]+(c[i][j]*Q[i][j]));
			
            //SOR
            
                p[i][j]=(1-w)*pp + w*p[i][j];

            //componentes del campo magnético

                bx[i][j]=(0.5*h)*(p[i][j+1]+p[i-1][j]-p[i+1][j]-p[i][j-1]);
                by[i][j]=-(0.5*h)*(p[i][j+1]-p[i-1][j]+p[i+1][j]-p[i][j-1]);

            //continuamos para escoger el máximo de las diferencias
            
				dif=fabs(p[i][j]-pp);
				
				if(dif>max)max=dif;
				
            }

        }

        k++;
    }
    while(max>tol);

printf("El número de iteraciones es: %d\n", k);

    FILE*arch=fopen("magnet.txt", "w");


    for( i=0 ; i<=n ; i++ ){

        for(j=0;j<=n;j++){
            fprintf(arch, "%lf %lf %lf %lf %lf \n", i*h , j*h , p[i][j], bx[i][j],by[i][j]);     
            }
            fprintf(arch, "\n");
    }

    fclose(arch);



    return 0;
}
