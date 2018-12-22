#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>

#define numP 10

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
    float fitness;
    int mem;
    int gain;
    int *Zi;
    int feasible;
}Vett;

void initialization(FILE *fin, Instances *in, Sol *best, Vett *pop);
void letturavet(int *v, FILE *fin, int r);
void letturamat(int **m, FILE *fin, int r, int c);
void calculateOF(Vett *v, Instances *in);
int relax3(Vett *v, Instances *in);
void GAinit(Vett *pop, Instances *in);
void createZi(Vett *v, Instances *in);
void searchMax(Vett *pop, Sol *best, Instances *in);
void crossover(Sol *p1, Sol *p2, Sol *son1, Sol *son2, Instances *in);
void calculateR(Vett *v, Instances *in);
void calculateFitness(Vett *v, Instances *in);
void Vett2Sol(Vett *v, Sol *s, Instances *in);
void averageGain(int *av,Instances *in);

int main(int argc, char* argv[])
{
    int j, c, q, i, iteration=0;
    Instances in;
    Sol best;
    FILE *fin, *fout;
    time_t start=time(NULL);
    int timelimit=0;
    Vett population[numP];

    assert(argc == 4);

    assert(strcmp(argv[2], "-t") == 0);

    timelimit = atoi(argv[3]);

    fin = fopen(argv[1], "r");

    assert(fin != NULL);

    initialization(fin, &in, &best, population);

    GAinit(population, &in);

    for(j=0; j<numP; j++)
    {
        for(q=0; q<in.Q; q++)
        {
            printf("%d ", population[j].vettX[q]);
        }
        printf("\nvettore zi%d\n ", j);
        for(i=0; i<in.I; i++)
        {
            printf("%d ", population[j].Zi[i]);
        }
        printf("\ngain: %d\nmemory: %d\nfeasible:%d \n", population[j].gain, population[j].mem,population[j].feasible);
    }


    best.gain=0;

    searchMax(population, &best, &in);

    while ((time(NULL) - start)<=timelimit)
    {

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

void initialization(FILE *fin, Instances *in, Sol *best, Vett *pop)
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

    for(c=0; c<numP; c++)
    {
        assert( (pop[c].vettX = malloc(in->Q*sizeof(int))) != NULL);
        assert( (pop[c].vettR = malloc(in->Q*sizeof(float))) != NULL);
        assert( (pop[c].Zi = malloc(in->Q*sizeof(int))) != NULL);
    }
}

void letturavet(int *v, FILE *fin, int r)
{
    int i;

    for(i=0; i<r; i++)
    {
        assert(fscanf(fin,"%d", &v[i]) != EOF);
    }
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
}

void calculateOF(Vett *v, Instances *in)
{
    int q, i, gain=0, cost=0;
    v->mem=0;

    for (q=0; q<in->Q; q++) {
        if (v->vettX[q] != -1) {
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
}

void GAinit(Vett *pop, Instances *in)
{
    int j=0, k=0,c=0, q=0;
    int todoN;
    int av=0;

    averageGain(&av, in);

    for(j=0; j<numP; j++)
    {
        pop[j].feasible=-1;
        pop[j].gain=-1;

        while(pop[j].feasible==0 || pop[j].gain <=0)
        {
            pop[j].feasible=-1;
            pop[j].gain=-1;
            //while(pop[j].feasible==0 || pop[j].gain <=0)
            //{

            todoN=0;
            for(q=0; q<in->Q; q++)
            {
                pop[j].vettX[q]=-1;
            }

            for(q=rand()%in->Q; todoN<in->Q;)
            {
                if(pop[j].vettX[q]==-1)
                {
                    c=rand()%in->C;

                    if(in->Gcq[c][q] >= av)
                    {
                        //temp->Xcq[c][q] = 1;
                        //q++;

                        for(k=0; k<in->Q; k++)
                        {
                            if(in->Gcq[c][k] >= av && pop[j].vettX[k]==-1)
                            {
                                pop[j].vettX[k] = c;
                                todoN++;
                            }
                        }
                        q=rand()%in->Q;
                    }
                }
                else
                {
                    q=rand()%in->Q;
                }

            }

            createZi(&pop[j], in);
            calculateOF(&pop[j], in);
            relax3(&pop[j], in);
        }
        //}
        calculateR(&pop[j], in);
        calculateFitness(&pop[j], in);
    }
}


void searchMax(Vett *pop, Sol *best, Instances *in)
{
    int i;
    for(i=0; i<numP; i++)
    {
        if(best->gain < pop[i].gain && pop[i].feasible==1)
        {
            Vett2Sol(&pop[i], best, in);
        }
    }
}

void crossover(Sol *p1, Sol *p2, Sol *son1, Sol *son2, Instances *in)
{
    int x,y;
    int c,q;
    x=rand()%in->Q;
    y=rand()%in->Q;

    for(q=0; q<in->Q; q++)
    {
        if(x>=y)
        {
            if(q>=x && q<=y)
            {
                for(c=0;c<in->C;c++)
                {
                    son1->Xcq[c][q]=p2->Xcq[c][q];
                    son2->Xcq[c][q]=p1->Xcq[c][q];
                }
            }
            else
            {
                for(c=0;c<in->C;c++)
                {
                    son1->Xcq[c][q]=p1->Xcq[c][q];
                    son2->Xcq[c][q]=p2->Xcq[c][q];
                }
            }
        }
        else
        {
            if(q>=y && q<=x)
            {
                for(c=0;c<in->C;c++)
                {
                    son1->Xcq[c][q]=p2->Xcq[c][q];
                    son2->Xcq[c][q]=p1->Xcq[c][q];
                }
            }
            else
            {
                for(c=0;c<in->C;c++)
                {
                    son1->Xcq[c][q]=p1->Xcq[c][q];
                    son2->Xcq[c][q]=p2->Xcq[c][q];
                }
            }
        }
    }
}

int relax3(Vett *v, Instances *in)
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

void calculateR(Vett *v, Instances *in)
// Calculates vettR from vettX
{
    int q, i;
    int sum;

    for(q=0; q<in->Q; q++) {
        for(i=0, sum=0; i<in->I; i++) {
            if(in->Eci[v->vettX[q]][i] == 1) {
                sum += in->Mi[i] + in->Fi[i];
            }
        }
        v->vettR[q] =(float) in->Gcq[v->vettX[q]][q] / sum;
    }
}

void calculateFitness(Vett *v, Instances *in)
{
    v->fitness = v->gain - relax3(v, in);
}

void Vett2Sol(Vett *v, Sol *s, Instances *in)
// Translate type Vett into type Sol solution
{
    int c, q;
    for(c=0; c<in->C; c++){
        for(q=0; q<in->Q; q++) {
            if(c == v->vettX[q]) {
                s->Xcq[c][q] = 1;
            } else {
                s->Xcq[c][q] = 0;
            }
        }
    }
    s->Zi = v->Zi;
    s->gain = v->gain;
    s->mem = v->mem;
    s->feasible = v->feasible;
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
        for(i=0; i<in->I; i++)
        {
            if(in->Eci[v->vettX[q]][i]==1)
            {
                v->Zi[i]=1;
            }
        }
    }
}

void averageGain(int *av,Instances *in)
{
    int c,q;
    for(c=0;c<in->C;c++)
    {
        for(q=0;q<in->Q;q++)
        {
            *av+=in->Gcq[c][q];
        }
    }
    *av= *av/(in->C*in->Q);

    return;

}
