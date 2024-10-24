//Program that obtains the nuclear magnetic moment for an odd proton when 1/2<j<13/2

#include<stdio.h>
#include<math.h>

#define muN 3.1525E-8

int main (){

    int l;
    double j, mu1, mu2, test;

    FILE*arch=fopen("magnetic.txt", "w");

    for(l=0;l<=6;l++){

        j=l+0.5;

        /* print the values of j to make sure i'm using the correct set of data
        
         printf("%.1f \n", j); */

        //j=l+0.5;

        //test=j*2.0;

        mu1 = ((j-0.5)+(3.3514/2.0));

        mu2 = (((j*(j+1.5))/(j+1.0))-((3.3514*j)/(2.0*(j+1))));

        fprintf(arch, "%.1lf %E %.1lf %E\n", j, mu1, j, mu2);


    }

    fclose(arch);

    return 0;

}

