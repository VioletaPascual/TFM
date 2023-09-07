#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#define VacunacionConstante
#define VacunacionVariable
//#define FuerzaInfeccion
#define Incidencia

//Cabeceras Funciones.
void evolucion(double *S, double *I, double *R, double *V, double lambda, double delta, double mu, double alpha);

double funcS(double S, double I, double R, double lambda, double delta, double alpha);
double funcI(double I, double S, double lambda, double mu);
double funcR(double R, double I, double mu, double delta, double alpha);
double funcV(double S, double R, double alpha);

double sumaVector(double *vector);
void igualarVectores(double *a, double *b);
void intercambioVectoresDouble(double *v1, double *v2);
void intercambioVectoresInt(int *v1, int *v2);
void ordenacionBurbuja(double *vector, int *indices);
void ordenaGrupos(double *vectorDesordenado, int *vectorIndices);


//Varables globales.
int Npasos=15000;//Número máximo de pasos de tiempo.
int Ngrupos=16; //Número de grupos de edad.

FILE *poblacion, *matrizContactos;
FILE *SIRV_sinVacuna, *SIRV_conVacuna, *fuerzaInfeccion, *incidencia, *IRR, *vacuna;
/*
FICHEROS CARGADOS (input, ya existentes):

    *poblacion: contiene el número de individuos de cada grupo de edad (16 números). Los valores aparecen en una sola línea separados por espacios.

    *matrizContactos: contiene la matriz cuadrada de contactos. Cada cuadro [i][j] corresponde al número de contactos que el grupo "i" reporta que tiene con el grupo "j".


FICHEROS CREADOS (output, ficheros de resultados):

    *SIRV: contiene el número total de susceptibles, infectados, recuperados y vacunados para cada paso de tiempo hasta el final de la simutacion t=Npasos.
    En SIRV_sinVacuna se registran los datos sin aplicar vacunación. En SIRV_conVacuna se guardan los valores aplicando la vacunacion a partir de Npasos/2.

    *fuerzaInfeccion: contiene el valor de la fuerza de infección para cada grupo en cada paso temporal.
    La fuerza de infección es la probabilidad de que un individuo infecte a otro.
    No es una probabilidad constante, es distinta para cada grupo de edad. lambda[i], "i" el grupo de edad.
    Depende del número de contactos del grupo "i" con el resto de grupos (en concreto depende de la frecuencia relativa de contactos, tanto por 1)
    y del estado del sistema (el número de infectados de todos los grupos en "t" se utilizan para calcular la fuerza de infección de "i" en "t+1").
    Para calcular lambda[i] en "t+1" es necesario conocer el número de infectados de cada grupo en "t", no solo los infectados del grupo "i".

    *incidencia: contiene el número de susceptibles que se infectan en un paso temporal, es decir, el numero de casos nuevos. Incidencia[i](t+1) = lambda[i](t) *S[i](t).

    *IRR (Incidence Rate Reduction): contiene la tasa de reducción de la incidencia, es decir, el porcentaje de reducción de casos al vacunar
    respecto a los casos que había sin vacunación. IRR = ((incidencia_sinVacuna - incidencia_conVacuna)/incidencia_sinVacuna)*100.

    *vacuna: contiene la probabilidad de vacunación para cada grupo en cada paso temporal. alpha[i].
*/

//Directorios de los ficheros comentados.
char directorio_1[200]="/Users/violetapascuallaborda/Desktop/TFG/piramidesPoblacionales/poblacionSriLanka.txt";
char directorio_2[200]="/Users/violetapascuallaborda/Desktop/TFG/matricesContactos/SriLankaMC_Homogeneizada_Normalizada.txt";
char directorio_3[200]="/Users/violetapascuallaborda/Desktop/TFG/SriLanka/SriLanka_SIRV_sinVacuna.txt";
char directorio_4[200]="/Users/violetapascuallaborda/Desktop/TFG/SriLanka/SriLanka_SIRV_conVacuna.txt";
char directorio_5[200]="/Users/violetapascuallaborda/Desktop/TFG/SriLanka/SriLanka_fuerzaInfeccion.txt";
char directorio_6[200]="/Users/violetapascuallaborda/Desktop/TFG/SriLanka/SriLanka_incidencia.txt";
char directorio_7[200]="/Users/violetapascuallaborda/Desktop/TFG/SriLanka/SriLanka_IRR.txt";
char directorio_8[200]="/Users/violetapascuallaborda/Desktop/TFG/SriLanka/SriLanka_vacunacion.txt";

//Cuerpo del código.
int main()
{
//Parámetros.
    double N[Ngrupos], S[Ngrupos], I[Ngrupos], R[Ngrupos], V[Ngrupos], I_aux[Ngrupos];//población, susceptibles, infectados, recuperados, vacunados por grupo de edad. I_aux es un vector necesario para no perder los infectados del paso de tiempo anterior.
    double k[Ngrupos][Ngrupos], lambda[Ngrupos], mu[Ngrupos], betha, gamma, delta;//frecuencia de contactos, fuerza de infección, probabilidad de recuperación.
    double vectorCamp[Ngrupos], incidencia_sinVacuna[Ngrupos], incidencia_conVacuna[Ngrupos];//vetorCampaña se utilizará para determinar el orden de la campaña de vacunación. Los otros dos vectores se utilizan para calcular el IRR.
    int vectorIndices[Ngrupos];
    int t, i, j, l, m;//Tiempo y distintos contadores.
    double a, b, c; //Constantes que describen el riesgo percibido usado para calcular la tasa vacunación.


#ifdef VacunacionConstante
    double alpha; //Probabilidad de vacunación.
    alpha=0.001;    //Porcentaje de vacunados. Se multiplica la población por alpha para dar el número de vacunados.
    //Este número se utiliza en caso de que la tasa de vacunación no dependa del estado del sistema y sea un porcentaje constante.
#endif // VacunacionConstante

#ifdef VacunacionVariable
    double alpha[Ngrupos]; //probabilidad de vacunación.
#endif // VacunacionVariable


//Rate constants.
    betha=0.75;  //Probabilidad de transmisión intrínseca por contacto.
    gamma=0.05;  //Constante usada para describir la probabilidad de que un infectado se recupere naturalmente.
    delta=0.1;   //Probabilidad de que un recuperado pierda la inmunidad y vuelva a ser susceptible.

    /*
    Se carga el número de personas de cada grupo (N[i]) con los datos de la pirámide poblacional.
    Se inicializan los infectados, susceptibles, recuperados y vacunados de cada grupo. Los recuperados y vacunados serán 0.
    Cargo también los datos de la matriz de contactos en k[i][j].
    k[i][j] representa la frecuencia relativa de contactos que un individuo del grupo "i" tiene con un individuo del grupo "j".
    */
    poblacion=fopen(directorio_1, "r");
    matrizContactos=fopen(directorio_2, "r");

    for(i=0; i<Ngrupos; i++)
        fscanf(poblacion, "%lf", &N[i]);

    for(i=0; i<Ngrupos; i++)
    {
        I[i]=0.001*sumaVector(N);
        S[i]=N[i]-I[i];
        R[i]=V[i]=0;

        mu[i]=gamma*exp(-0.025*i);

        for(j=0; j<Ngrupos; j++)
            fscanf(matrizContactos, "%lf", &k[i][j]);
    }

    fclose(poblacion);
    fclose(matrizContactos);

    a=sumaVector(N)/2.0;
    b=2.5;
    c=3;

    /*
    Se utiliza un vector auxiliar I_aux[i] para guardar los infectados del paso temporal "t".
    Para calcular I[i] en "t+1" se usa el número de infectados de todos los grupos en "t".
    Se necesita un vector auxiliar ya que el vector I[i] se actualiza para cada grupo.
    Si no existiera el vector I_aux[i], estaríamos calculando los nuevos infectados de un grupo en "t+1" con valores de infectados en "t+1" en lugar de calcularlos con valores en "t".
    */
    igualarVectores(I, I_aux);

//Se abren los ficheros de resultados.
    SIRV_sinVacuna=fopen(directorio_3, "wt");
    SIRV_conVacuna=fopen(directorio_4, "wt");
    fuerzaInfeccion=fopen(directorio_5, "wt");
    incidencia=fopen(directorio_6, "wt");
    IRR=fopen(directorio_7, "wt");
    vacuna=fopen(directorio_8, "wt");

//Se escribe en la primera línea los valores iniciales.
    fprintf(SIRV_sinVacuna, "%d %lf %lf %lf %lf\n", 0, sumaVector(S), sumaVector(I), sumaVector(R), sumaVector(V));
    fprintf(SIRV_conVacuna, "%d %lf %lf %lf %lf\n", 0, sumaVector(S), sumaVector(I), sumaVector(R), sumaVector(V));
    fprintf(fuerzaInfeccion, "%d ", 1);
    fprintf(incidencia, "%d ", 1);

    for(i=0; i<Ngrupos; i++)
    {
        fprintf(fuerzaInfeccion, "%lf ", lambda[i]);
        fprintf(incidencia, "%lf ", lambda[i]*S[i]);
    }

    fprintf(fuerzaInfeccion, "\n");
    fprintf(incidencia, "\n");

//Nucleo del programa.

//Evolución de la epidemia sin vacuna.
    for(t=1; t<Npasos; t++)
    {
        fprintf(fuerzaInfeccion, "%d ", t+1);
        fprintf(incidencia, "%d ", t+1);

        for(i=0; i<Ngrupos; i++)
        {
            lambda[i]=0; //Para no calcular el nuevo lambda[i] sobre el anterior lambda[i]. Sumar los nuevos incrementos al lambda anterior.

            for(j=0; j<Ngrupos; j++)
                lambda[i]+=k[i][j]*I_aux[j]/N[j];

            lambda[i]=betha*lambda[i];

            evolucion(&S[i], &I[i], &R[i], &V[i], lambda[i], mu[i], delta, 0);

            fprintf(fuerzaInfeccion, "%lf ", lambda[i]);
            fprintf(incidencia, "%lf ", lambda[i]*S[i]);

            if (t==Npasos-1)
            {
                incidencia_sinVacuna[i]=lambda[i]*S[i];
#ifdef FuerzaInfeccion
                vectorCamp[i]=lambda[i];
#endif // FuerzaInfeccion

#ifdef Incidencia
                vectorCamp[i]=incidencia_sinVacuna[i];
#endif // Incidencia
            }
        }

        if(((sumaVector(I))<0.01))
            break;

        igualarVectores(I, I_aux);

        fprintf(SIRV_sinVacuna, "%d %lf %lf %lf %lf\n", t, sumaVector(S), sumaVector(I), sumaVector(R), sumaVector(V));

        fprintf(fuerzaInfeccion, "\n");
        fprintf(incidencia, "\n");
    }

//Evolución de la epidemia con vacuna. Vacunando en orden de mayor fuerza de infección o de incidencia dependiendo del #define.
    for (l=0; l<Ngrupos; l++)
    {
        for(i=0; i<Ngrupos; i++)
        {
            I[i]=0.001*sumaVector(N);
            S[i]=N[i]-I[i];
            R[i]=V[i]=0;
        }

        igualarVectores(I, I_aux);

        for(t=1; t<Npasos; t++)
        {

#ifdef VacunacionVariable
            if(l==15)
                if (t>=(5000))
                    fprintf(vacuna, "%d ", t);
#endif // VacunacionVariable

            for(i=0; i<Ngrupos; i++)
            {
                lambda[i]=0; //Para no calcular el nuevo lambda[i] sobre el anterior lambda[i]. Sumar los nuevos incrementos al lambda anterior.

                for(j=0; j<Ngrupos; j++)
                    lambda[i]+=k[i][j]*I_aux[j]/N[j];

                lambda[i]=betha*lambda[i];

                if (t<(5000))
                    evolucion(&S[i], &I[i], &R[i], &V[i], lambda[i], mu[i], delta, 0);
                else
                {
                    int conditionMet = 0; // Variable para controlar si se cumple la condición del if

                    for (m=0; m<l+1; m++)
                        if (i==vectorIndices[m])
                        {
                            conditionMet = 1;
                            break;
                        }

                    if (conditionMet)
                    {
#ifdef VacunacionConstante
                        evolucion(&S[i], &I[i], &R[i], &V[i], lambda[i], mu[i], delta, alpha);
#endif // VacunacionConstante

#ifdef VacunacionVariable
                        alpha[i]=lambda[i]*(1-exp(-pow(sumaVector(I_aux)/a,b)))/(1+c*sumaVector(V)/sumaVector(N));
                        evolucion(&S[i], &I[i], &R[i], &V[i], lambda[i], mu[i], delta, alpha[i]);
                        if(l==15)
                            fprintf(vacuna, "%lf ", alpha[i]);
#endif // VacunacionVariable
                    }
                    else
                        evolucion(&S[i], &I[i], &R[i], &V[i], lambda[i], mu[i], delta, 0);


                }
                if (t==Npasos-1)
                    incidencia_conVacuna[i]=lambda[i]*S[i];
            }

            //if(((sumaVector(I))<0.01))
                //break;

            igualarVectores(I, I_aux);

            //Solo guarda los datos para todos los grupos vacunados. Se puede cambiar el if para elegir cuantos grupos vacunados.
            if(l==15)
            {
                fprintf(SIRV_conVacuna, "%d %lf %lf %lf %lf\n", t, sumaVector(S), sumaVector(I), sumaVector(R), sumaVector(V));
#ifdef VacunacionVariable
                if (t>=(5000))
                    fprintf(vacuna, "\n");
#endif // VacunaVariable
            }

        }

        fprintf(IRR, "%lf ", ((sumaVector(incidencia_sinVacuna)-sumaVector(incidencia_conVacuna))/sumaVector(incidencia_sinVacuna))*100);
    }

//Se cierran los ficheros.
    fclose(SIRV_sinVacuna);
    fclose(SIRV_conVacuna);
    fclose(fuerzaInfeccion);
    fclose(incidencia);
    fclose(IRR);
    fclose(vacuna);
    return 0;
}


//Runge kutta 4.
void evolucion(double *S, double *I, double *R, double *V, double lambda, double mu, double delta, double alpha)
{
    double k1, k2, k3, k4, k, l1, l2, l3, l4, l, m1, m2, m3, m4, m, n1, n2, n3, n4, n;
    double h;
    h=0.1;

    k1=funcS(*S, *I, *R, lambda, delta, alpha);
    l1=funcI(*I, *S, lambda, mu);
    m1=funcR(*R, *I, mu, delta, alpha);
    n1=funcV(*S, *R, alpha);

    k2=funcS((*S+k1*0.5), (*I+l1*0.5), (*R+m1*0.5), lambda, delta, alpha);
    l2=funcI((*I+l1*0.5), (*S+k1*0.5), lambda, mu);
    m2=funcR((*R+m1*0.5), (*I+l1*0.5), mu, delta, alpha);
    n2=funcV((*S+k1*0.5), (*R+m1*0.5), alpha);

    k3=funcS((*S+k2*0.5), (*I+l2*0.5), (*R+m2*0.5), lambda, delta, alpha);
    l3=funcI((*I+l2*0.5), (*S+k2*0.5), lambda, mu);
    m3=funcR((*R+m2*0.5), (*I+l2*0.5), mu, delta, alpha);
    n3=funcV((*S+k2*0.5), (*R+m2*0.5), alpha);

    k4=funcS((*S+k3), (*I+l3), (*R+m3), lambda, delta, alpha);
    l4=funcI((*I+l3), (*S+k3), lambda, mu);
    m4=funcR((*R+m3), (*I+l3), mu, delta, alpha);
    n4=funcV((*S+k3), (*R+m3), alpha);

    k=((k1+2*k2+2*k3+k4)/6.0);
    l=((l1+2*l2+2*l3+l4)/6.0);
    m=((m1+2*m2+2*m3+m4)/6.0);
    n=((n1+2*n2+2*n3+n4)/6.0);

    *S=*S+k*h;
    *I=*I+l*h;
    *R=*R+m*h;
    *V=*V+n*h;

    if(*S<0)
        *S=0;
    if(*I<0)
        *I=0;
    if(*R<0)
        *R=0;
    if(*V<0)
        *V=0;
}

double funcS(double S, double I, double R, double lambda, double delta, double alpha)
{
    return -lambda*S + delta*R - alpha*S;
}

double funcI(double I, double S, double lambda, double mu)
{
    return lambda*S - mu*I;
}

double funcR(double R, double I, double mu, double delta, double alpha)
{
    return mu*I - delta*R - alpha*R;
}

double funcV(double S, double R, double alpha)
{
    return alpha*(S+R);
}

double sumaVector(double *vector)
{
    double suma;
    int n;

    suma=0;
    for(n=0; n<Ngrupos; n++)
        suma+=vector[n];

    return suma;
}

void igualarVectores(double *a, double *b)
{
    int n;
    for(n=0; n<Ngrupos; n++)
        b[n]=a[n];
}

void intercambioVectoresDouble(double *v1, double *v2)
{
    double temporal = *v1;
    *v1 = *v2;
    *v2 = temporal;
}

void intercambioVectoresInt(int *v1, int *v2)
{
    int temporal = *v1;
    *v1 = *v2;
    *v2 = temporal;
}
void ordenacionBurbuja(double *vector, int *indices)
{
    int i, j;
    for (i = 0; i < Ngrupos - 1; i++)
        for (j = 0; j < Ngrupos - i - 1; j++)
            if (vector[j] < vector[j + 1])
            {
                intercambioVectoresDouble(&vector[j], &vector[j + 1]);
                intercambioVectoresInt(&indices[j], &indices[j + 1]);
            }
}

//La función ordena los índices de *vectorDesdenado en base a los valores de mayor a menor y los almacena en *vectorIndices.
void ordenaGrupos(double *vectorDesordenado, int *vectorIndices)
{
    // Inicializar el vector de índices
    for (int i=0; i<Ngrupos; i++)
        vectorIndices[i] = i;

    // Ordenar los índices en función de los valores correspondientes en el vector desordenado
    ordenacionBurbuja(vectorDesordenado, vectorIndices);
}
