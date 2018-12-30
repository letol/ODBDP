#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>

#define true 1
#define false 0
#define numN 10

typedef struct instances{
    int Q, I, C, M;
    int **Eci;
    int *Fi, *Mi;
    int **Gcq;
} Instances;

typedef struct sol{
    int **Xcq;
    int *Zi;
    int mem;
    int gain;
    int feasible;
}Sol;

typedef struct vett{
    int *vettX;
    float *vettR;
    int mem;
    int gain;
    int *Zi;
    int feasible;
}Vett;

void initialization(FILE *fin, Sol *best, Vett *pop, Vett *neighbor, Instances *in);
void letturavet(int *v, FILE *fin, int r);
void letturamat(int **m, FILE *fin, int r, int c);


void TSinit(Vett *pop, Instances *in);
void createZi(Vett *v, Instances *in);
void calculateOF(Vett *v, Instances *in);
int searchCforQ(int q, int *confmem, Instances *in);
int check3(Vett *v, Instances *in);
void changeConfMem(int *confmem, Vett *x, Instances *in);

void Vett2Sol(Vett *v, Sol *s, Instances *in);

int main(int argc, char* argv[])
{
    int c, q, i, iteration=0;
    Instances in;
    Sol best;
    FILE *fin, *fout;
    time_t start=time(NULL);
    int timelimit=0;
    Vett temp, *nei;


    assert(argc == 4);

    assert(strcmp(argv[2], "-t") == 0);

    timelimit = atoi(argv[3]);

    fin = fopen(argv[1], "r");

    assert(fin != NULL);

    assert( (nei = malloc(numN*sizeof(Vett))) != NULL);

    initialization(fin, &best, &temp, nei, &in);

    TSinit(&temp, &in);

    for(q=0; q<in.Q; q++)
    {
        printf("%d ", temp.vettX[q]);
    }
    printf("\nvettore zi:\n ");
    for(i=0; i<in.I; i++)
    {
        printf("%d ", temp.Zi[i]);
    }
    printf("\ngain: %d\nmemory: %d\nfeasible:%d \n", temp.gain, temp.mem, temp.feasible);

    Vett2Sol(&temp, &best, &in);

    while ((time(NULL) - start)<=timelimit)
    {
        iteration++;
    }
    /*
        while(clock()<30000){}
    */

    printf("Best solution is:\n");
    for(c=0; c<in.C; c++)
    {
        for(q=0; q<in.Q; q++)
        {
            printf("%d ", best.Xcq[c][q]);
        }
        printf("\n");
    }
    printf("\n");
    for(i=0; i<in.I; i++)
    {
        printf("%d ", best.Zi[i]);
    }

    printf("\ngain: %d\nmemory: %d\nfeasible: %d", best.gain, best.mem, best.feasible);
    printf("\n\nIteration: %d", iteration);

    //TODO: Ricordiamoci delle free

    return 0;
}

void initialization(FILE *fin, Sol *best, Vett *pop, Vett *neighbor, Instances *in)
{
    char lecture[40];
    int c=0;

    assert( (fscanf(fin, "%s %d", lecture, &in->Q)) != EOF );
    assert( (fscanf(fin, "%s %d", lecture, &in->I)) != EOF );
    assert( (fscanf(fin, "%s %d", lecture, &in->C)) != EOF );
    assert( (fscanf(fin, "%s %d", lecture, &in->M)) != EOF );

    assert((in->Eci=malloc(in->C*sizeof(int*))) != NULL);

    for(c=0; c<in->C; c++)
    {
        assert( (in->Eci[c]=malloc(in->I*sizeof(int))) !=NULL  );
    }

    assert( (in->Fi = malloc(in->I*sizeof(int))) !=NULL );
    assert( (in->Mi = malloc(in->I*sizeof(int))) !=NULL );

    assert( (in->Gcq = malloc(in->C*sizeof(int*))) != NULL );

    for(c=0; c<in->C; c++)
    {
        assert( (in->Gcq[c]=malloc(in->Q*sizeof(int))) !=NULL  );
    }

    assert( (best->Xcq = malloc(in->C*sizeof(int*))) != NULL);

    for(c=0; c<in->C; c++)
    {
        assert( (best->Xcq[c]=malloc(in->Q*sizeof(int))) !=NULL  );
    }

    assert( (best->Zi = malloc(in->I*sizeof(int))) !=NULL );

    assert( (fscanf(fin, "%s", lecture)) != EOF );

    letturamat(in->Eci, fin, in->C, in->I);

    assert( (fscanf(fin, "%s", lecture)) != EOF );

    letturavet(in->Fi, fin, in->I);

    assert( (fscanf(fin, "%s", lecture)) != EOF );

    letturavet(in->Mi, fin, in->I);

    assert( (fscanf(fin, "%s", lecture)) != EOF );

    letturamat(in->Gcq, fin, in->C, in->Q);

    assert( (pop->vettX = malloc(in->Q*sizeof(int))) != NULL);
    assert( (pop->vettR = malloc(in->Q*sizeof(float))) != NULL);
    assert( (pop->Zi = malloc(in->I*sizeof(int))) != NULL);

    for(c=0; c<numN; c++)
    {
        assert( (neighbor[c].vettX = malloc(in->Q*sizeof(int))) != NULL);
        assert( (neighbor[c].vettR = malloc(in->Q*sizeof(float))) != NULL);
        assert( (neighbor[c].Zi = malloc(in->I*sizeof(int))) != NULL);
    }

    return;
}


void letturavet(int *v, FILE *fin, int r)
{
    int i;

    for(i=0; i<r; i++)
    {
        assert(fscanf(fin,"%d", &v[i]) != EOF);
    }

    return;
}

void letturamat(int **m, FILE *fin,int r, int c)
{
    int i,j;

    for(i=0; i<r; i++)
    {
        for(j=0; j<c; j++)
        {
            assert(fscanf(fin,"%d", &m[i][j]) != EOF);
        }
    }
    return;
}

void TSinit(Vett *pop, Instances *in)
//questa inizializzazione permette di avere una soluzione che utilizzi il maggior numero di memoria possibile controllando un certo
//numero di configurazioni casuali che verranno aggiunte ad una query o sostituite dove il gain risulta maggiore
{
    int i, c, q, count;
    int *confmem; //vettore del costo di una configurazione. Si modifica in base agli indici attualmente attivi tramite la funzione changeConfMem
    int Memory, numconf; //memory tiene il conto della memoria attualmente disponibile e permette di verificare che una configurazione non causi infeasibility
    pop->feasible=-1;
    pop->gain=-1;

    assert(confmem=malloc(in->C*sizeof(int)));

    while(pop->feasible==0 || pop->gain <=0)
    {
        count=0;
        Memory=in->M;
        numconf=0;

        for(c=0; c<in->C; c++)
        {
            confmem[c]=0;
            for(i=0; i<in->I; i++)
            {
                if(in->Eci[c][i]==1)
                {
                    confmem[c]+=in->Mi[i]; //init confmem
                }
            }
        }

        for(q=0; q<in->Q;q++)
        {
            pop->vettX[q]=-1;
        }

        createZi(pop, in);

        for(c=rand()%in->C; count<(in->C/5); count++, c=rand()%in->C)
        {
            if( (Memory-confmem[c]) > 0) //se non supero la memoria disponibile.
            {
                numconf++;
                for(q=0;q<in->Q;q++)
                {
                    if((pop->vettX[q]) == -1 && (in->Gcq[c][q]) > 0 ) //assegno la c a q se fornisce guadagno e la query non è servita
                    {
                        pop->vettX[q]=c;
                    }
                    else if( (pop->vettX[q]) > -1)
                    {
                        if(in->Gcq[c][q] > in->Gcq[pop->vettX[q]][q])   //se la q è servita da una c che fornisce un guadagno minore, sostituiscila
                        {
                            pop->vettX[q]=c;
                        }
                    }
                }

                Memory-=confmem[c];
                createZi(pop, in);
                changeConfMem(confmem, pop, in);
                calculateOF(pop, in);
                check3(pop, in);
            }
        }
    }
    printf("configuzioni usate %d\n", numconf);
    return;
}

void createZi(Vett *v, Instances *in)
{
    int i, q;

    for(i=0; i<in->I; i++)
    {
        v->Zi[i]=0;
    }

    for(q=0; q<in->Q; q++)
    {
        if(v->vettX[q]>=0)
        {
            for(i=0; i<in->I; i++)
            {
                if(in->Eci[v->vettX[q]][i]==1)
                {
                    v->Zi[i]=1;
                }
            }
        }
    }

    return;
}

void calculateOF(Vett *v, Instances *in)
{
    int q, i, gain=0, cost=0;
    v->mem=0;

    for (q=0; q<in->Q; q++) {
        if (v->vettX[q] >= 0) {
            gain += in->Gcq[v->vettX[q]][q];
        }
    }

    for (i=0; i<in->I; i++) {
        if (v->Zi[i] == 1) {
            cost += in->Fi[i];
            v->mem+=in->Mi[i];
        }
    }

    v->gain = (gain-cost);

    return;
}

void Vett2Sol(Vett *v, Sol *s, Instances *in)
// Translate type Vett into type Sol solution
{
    int c, q, i;
    for(c=0; c<in->C; c++){
        for(q=0; q<in->Q; q++) {
            if(c == v->vettX[q]) {
                s->Xcq[c][q] = 1;
            } else {
                s->Xcq[c][q] = 0;
            }
        }
    }

    for(i=0;i<in->I;i++)
    {
        s->Zi[i] = v->Zi[i];
    }

    s->gain = v->gain;
    s->mem = v->mem;
    s->feasible = v->feasible;

    return;
}

int searchCforQ(int q, int *confmem, Instances *in)
{
    float temp=0, best=0;
    int c, confB=0;

    for(c=0; c<in->C; c++)
    {
        if(confmem[c]==0 && in->Gcq[c][q]>0)
        {
            if(best > (temp=(1/in->Gcq[c][q])) )
            {
                confB=c;
                best=temp;
            }

        }
        else if(in->Gcq[c][q]>0)
        {
            if(best > (temp=(confmem[c]/in->Gcq[c][q])) )
            {
                confB=c;
                best=temp;
            }
        }
    }


    return confB;
}

int check3(Vett *v, Instances *in)
// relaxation of constraint 3
// returns 0 if limit M is not exceeded, otherwise it returns the excess amount, scaled by a factor.
{
    //TODO: aggiungere fattore di scalamento

    if(v->mem > in->M)
    {
        v->feasible=0;
        return v->mem - in->M;
    }
    v->feasible=1;
    return 0;
}

void changeConfMem(int *confmem, Vett *x, Instances *in)
{
    int i, j;

    for(j=0; j<in->C; j++)
    {
        confmem[j]=0;
        for(i=0; i<in->I; i++)
        {
            if(in->Eci[j][i]==1)
            {
                confmem[j]+=in->Mi[i];
            }
        }
    }

    for(i=0; i<in->I; i++)
    {
        for(j=0; j<in->C; j++)
        {
            if(in->Eci[j][i] == 1 && x->Zi[i]==1)
            {
                confmem[j]-=in->Mi[i];
            }
        }

    }
}