#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
int var3;
int i;

struct variabili{
    int x;
    int y;
    int z;
};
int diff(int x, int y){
    return x-y;
}

DWORD WINAPI GeneticA(struct variabili *tr1)
{
    FILE *fp;
    fp=fopen("output.txt","w");
    int var1,var2;
    var1=50;
    var2=30;
    var3 = diff(var1,var2);
    tr1->z=tr1->x*tr1->y;
    while(1) {
        fprintf(fp,"ciaone1");
        rewind(fp);
        printf("eseguo genetic %d\n", i);
        i++;
        Sleep(1000);

    }
}

DWORD WINAPI TabuSearch(struct variabili *tr2)
{
    FILE *fp;
    fp=fopen("output.txt","w");
    tr2->z=tr2->x*tr2->y;
    while(1) {
        fprintf(fp,"ciaone2");
        printf("eseguo tabu %d\n", i);
        rewind(fp);
        i++;
        Sleep(1000);
    }
}
int main(int argc, char*argv[]) {
    int t=10;
    t=t*1000; //da secondi a millisecondi
    struct variabili tr1,tr2;
    tr1.x=5;
    tr1.y=3;
    tr2.x=8;
    tr2.y=2;
    HANDLE thread1 = CreateThread(NULL,//attributo sicurezza
                                  0,//dimensione stack default
                                  (LPTHREAD_START_ROUTINE)GeneticA,//funzione thread
                                  &tr1,//argomento thread
                                  0,//flag creazione
                                  NULL//ritorno thread
    );

    HANDLE thread2 = CreateThread(NULL,//attributo sicurezza
                                  0,//dimensione stack default
                                  (LPTHREAD_START_ROUTINE)TabuSearch,//funzione thread
                                  &tr2,//argomento thread
                                  0,//flag creazione
                                  NULL//ritorno thread
    );
    //lettura del file e caricamento della struttura dati->
    Sleep(t);
    printf("termino thread1\n");
    TerminateThread(thread1, 0);

    Sleep(5000);
    printf("termino thread2\n");
    TerminateThread(thread2, 0);
    printf("programma terminato con %d e %d",tr1.z,tr2.z);
    //alla return del programma chiamante i thread sembra che smettano di girare ma non ne sono sicuro, nel dubbio
    //si chiudono a mano
    // eventuali free e termine del programma vanno fatte qui poichè i thread vengono chiusi senza preavviso.
    //Variabili allocate dinamicamente dentro al thread non ho idea di come farci una free prima di uscire.
    //Suppongo che nel momento in cui vengano terminati i thread implicitamente le free vengano fatte
    return 0;
}