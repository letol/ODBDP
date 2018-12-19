#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>

#define true 1
#define false 0
#define numP 8

typedef struct instances{
    int Q, I, C, M;
    int **Eci;
    int *Fi, *Mi;
    int **Gcq;
} Instances;

typedef struct sol{
    int **Xcq;
    //TODO: int *Yc; //Dobbiamo aggiungerlo?
    int *Zi;
    int mem;
    int gain;
    int feasible;
}Sol;

void initialization(FILE *fin, Instances *in, Sol *best, Sol *temp);
void letturavet(int *v, FILE *fin, int r);
void letturamat(int **m, FILE *fin, int r, int c);
void calculateOF(Sol *temp, Instances *in);
int check(Sol *temp, Instances *in);
//int check1(Sol *temp, Instances *in);
//int check2(Sol *temp, Instances *in);
//int check3(Sol *temp, Instances *in);
void GAinit(Sol *pop, Sol* newpop, Instances *in);
void createZi(Sol *temp, Instances *in);
void avoidC2(Sol *temp, Instances *in);
void searchMax(Sol *pop, Sol *best, Instances *in);
void cpySol(Sol *best, Sol *temp, Instances *in);
void crossoverC(Sol *p1, Sol *p2, Sol *son1, Sol *son2, Instances *in);
void crossoverR(Sol *p1, Sol *p2, Sol *son1, Sol *son2, Instances *in);
void changePop(Sol *pop1, Sol *pop2, Instances *in);

int main(int argc, char* argv[])
{
    int j, c, q, i, flag=0, iteration=0;
    Instances in;
    Sol best, temp;
    FILE *fin, *fout;
    time_t start=time(NULL);
    int timelimit=0;
    Sol population[numP], newpop[numP];

    assert(argc == 4);

    assert(strcmp(argv[2], "-t") == 0);

    timelimit = atoi(argv[3]);

    fin = fopen(argv[1], "r");

    assert(fin != NULL);

    initialization(fin, &in, &best, &temp);

    GAinit(population, newpop, &in);
/*
    for(j=0; j<numP; j++)
    {
        for(c=0; c<in.C; c++)
        {
            for(q=0; q<in.Q; q++)
            {
                printf("%d ", population[j].Xcq[c][q]);
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");

    for(j=0; j<numP; j++)
    {
        printf("vector Z%d: ", j);
        for(i=0; i<in.I; i++)
        {
            printf("%d ", population[j].Zi[i]);
        }
        printf("\ngain: %d\nmemory: %d\nfeasible: %d", population[j].gain, population[j].mem, population[j].feasible);
        printf("\n\n");
    }
    printf("\n");
*/
    best.gain=0;

    searchMax(population, &best, &in);

    while ((time(NULL) - start)<=timelimit){

        for(j=0; j<numP; j+=2)
        {
            flag=rand()%2;
            if(flag==0)
            {
                crossoverC(&population[j], &population[j+1], &newpop[j], &newpop[j+1], &in);

                createZi(&newpop[j], &in);
                createZi(&newpop[j+1], &in);

                calculateOF(&newpop[j], &in);
                calculateOF(&newpop[j+1], &in);

                newpop[j].feasible=check(&newpop[j], &in);
                newpop[j+1].feasible=check(&newpop[j+1], &in);
            }
            else
            {
                crossoverR(&population[j], &population[j+1], &newpop[j], &newpop[j+1], &in);

                createZi(&newpop[j], &in);
                createZi(&newpop[j+1], &in);

                calculateOF(&newpop[j], &in);
                calculateOF(&newpop[j+1], &in);

                newpop[j].feasible=check(&newpop[j], &in);
                newpop[j+1].feasible=check(&newpop[j+1], &in);
            }

        }

        changePop(population, newpop, &in);
        searchMax(population, &best, &in);
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

void initialization(FILE *fin, Instances *in, Sol *best, Sol *temp)
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

    assert( (temp->Xcq = malloc(in->C*sizeof(int*))) != NULL);

    for(c=0; c<in->C; c++)
    {
        assert( (temp->Xcq[c]=malloc(in->Q*sizeof(int))) !=NULL  );
    }

    assert( (temp->Zi = malloc(in->I*sizeof(int))) !=NULL );

    assert( (fscanf(fin, "%s", lecture)) != EOF );

    letturamat(in->Eci, fin, in->C, in->I);

    assert( (fscanf(fin, "%s", lecture)) != EOF );

    letturavet(in->Fi, fin, in->I);

    assert( (fscanf(fin, "%s", lecture)) != EOF );

    letturavet(in->Mi, fin, in->I);

    assert( (fscanf(fin, "%s", lecture)) != EOF );

    letturamat(in->Gcq, fin, in->C, in->Q);

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

void calculateOF(Sol *temp, Instances *in)
{
    int c, q, i, gain=0, cost=0;
    temp->mem=0;

    for (c=0; c<in->C; c++) {
        for (q=0; q<in->Q; q++) {
            if (temp->Xcq[c][q] == 1) {
                gain += in->Gcq[c][q];
            }
        }
    }

    for (i=0; i<in->I; i++) {
        if (temp->Zi[i] == 1) {
            cost += in->Fi[i];
            temp->mem+=in->Mi[i];
        }
    }

    temp->gain = (gain-cost);

    return;
}

int check(Sol *temp, Instances *in)
// Returns 1 if temp is feasible, 0 otherwise.
{
    int c, q;
    int con2=0;

    //Constraint (2)
    for(q=0; q<in->Q; q++)
    {
        con2=0;
        for(c=0; c<in->C; c++)
        {
            con2+=temp->Xcq[c][q];
        }
        if(con2>1)
        {
            return 0;
        }
    }

    //Constraint (3)
    if(temp->mem > in->M)
    {
        return 0;
    }
    return 1;
}

void GAinit(Sol* pop, Sol* newpop, Instances *in)
{
    int c=0, q=0, j=0;

    for(j=0; j<numP; j++)
    {
        assert( (pop[j].Xcq = malloc(in->C*sizeof(int*))) != NULL);

        for(c=0; c<in->C; c++)
        {
            assert( (pop[j].Xcq[c]=malloc(in->Q*sizeof(int))) !=NULL  );
        }

        assert( (pop[j].Zi = malloc(in->I*sizeof(int))) !=NULL );

        assert( (newpop[j].Xcq = malloc(in->C*sizeof(int*))) != NULL);

        for(c=0; c<in->C; c++)
        {
            assert( (newpop[j].Xcq[c]=malloc(in->Q*sizeof(int))) !=NULL  );
        }

        assert( (newpop[j].Zi = malloc(in->I*sizeof(int))) !=NULL );
    }
    for(j=0; j<numP; j++)
    {
        for(c=0; c<in->C; c++)
        {
            for(q=0; q<in->Q; q++)
            {
                if(in->Gcq[c][q]!=0)
                {
                    pop[j].Xcq[c][q]=rand()%2;
                }
                else
                {
                    pop[j].Xcq[c][q]=0;
                }
            }
        }
/*
        for(q=0; q<in->Q; q++)
        {
            c=rand()%(in->C);
            for(c=rand()%(in->C), x=0; pop[j].Xcq[c][q]==0 && x<50; c=rand()%(in->C), x++);

            if (pop[j].Xcq[c][q] == 1)
            {
                for(k=0; k<in->C; k++)
                {
                    if(k!=c)
                    {
                        pop[j].Xcq[k][q]=0;
                    }
                }
            }
        }
*/
        avoidC2(&pop[j], in);
        createZi(&pop[j], in);
        calculateOF(&pop[j], in);
        pop[j].feasible=check(&pop[j], in);
    }

    return;
}


void createZi(Sol* temp, Instances *in)
{
    int i, c, q;


    for(i=0; i<in->I; i++)
    {
        temp->Zi[i]=0;
    }

    for(c=0; c<in->C; c++)
    {
        for(q=0; q<in->Q; q++)
        {
            for(i=0; i<in->I; i++)
            {
                if(temp->Xcq[c][q]==1 && in->Eci[c][i]==1)
                {
                    temp->Zi[i]=1;
                }
            }
        }
    }
    return;
}

void avoidC2(Sol *temp, Instances *in)
{
    int c,q,x,k;

    for(q=0; q<in->Q; q++)
    {
        c=rand()%(in->C);
        for(c=rand()%(in->C), x=0; temp->Xcq[c][q]==0 && x<50; c=rand()%(in->C), x++);

        if (temp->Xcq[c][q] == 1)
        {
            for(k=0; k<in->C; k++)
            {
                if(k!=c)
                {
                    temp->Xcq[k][q]=0;
                }
            }
        }
    }
}

void searchMax(Sol *pop, Sol *best, Instances *in)
{
    int i;
    for(i=0; i<numP; i++)
    {
        if(best->gain <= pop[i].gain && pop[i].feasible==1)
        {
            cpySol(best, &pop[i], in);
        }
    }
}

void cpySol(Sol *best, Sol *temp, Instances *in)
{
    int i, j;
    for(i=0; i<in->C; i++)
    {
        for(j=0; j<in->Q; j++)
        {
            best->Xcq[i][j]=temp->Xcq[i][j];
        }
    }

    for(i=0; i<in->I; i++)
    {
        best->Zi[i]=temp->Zi[i];
    }
    best->gain=temp->gain;
    best->mem=temp->mem;
    best->feasible=temp->feasible;
}

void crossoverC(Sol *p1, Sol *p2, Sol *son1, Sol *son2, Instances *in)
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

void crossoverR(Sol *p1, Sol *p2, Sol *son1, Sol *son2, Instances *in)
{
    int x,y;
    int c,q;
    x=rand()%in->Q;
    y=rand()%in->Q;

    for(c=0; c<in->C; c++)
    {
        if(x>=y)
        {
            if(c>=x && c<=y)
            {
                for(q=0;q<in->Q;q++)
                {
                    son1->Xcq[c][q]=p2->Xcq[c][q];
                    son2->Xcq[c][q]=p1->Xcq[c][q];
                }
            }
            else
            {
                for(q=0;q<in->Q;q++)
                {
                    son1->Xcq[c][q]=p1->Xcq[c][q];
                    son2->Xcq[c][q]=p2->Xcq[c][q];
                }
            }
        }
        else
        {
            if(c>=x && c<=y)
            {
                for(q=0;q<in->Q;q++)
                {
                    son1->Xcq[c][q]=p2->Xcq[c][q];
                    son2->Xcq[c][q]=p1->Xcq[c][q];
                }
            }
            else
            {
                for(q=0;q<in->Q;q++)
                {
                    son1->Xcq[c][q]=p1->Xcq[c][q];
                    son2->Xcq[c][q]=p2->Xcq[c][q];
                }
            }
        }
    }
}

void changePop(Sol *pop1, Sol *pop2, Instances *in)
{
    int i;
    for(i=0;i<numP;i++)
    {
        cpySol(&pop1[i], &pop2[i], in);
    }
}
/*
int check(Sol *temp, Instances *in)
// Returns 1 if temp is feasible, 0 otherwise.
{
    if (check1(temp,in) > 0 ||
        check2(temp,in) > 0 ||
        check3(temp,in) >0) {
        return 1;
    } else {
        return 0;
    }
}

int check1(Sol *temp, Instances *in)
// check constraint 1
// returns the number of indexes that should be present and are not, scaled by a factor.
{
    int c, q, i;
    int p1=0;
    //TODO: aggiungere fattore di scalamento

    for (c=0; c<in->C; c++) {                                       //TODO: eventualmente da riscrivere usando Yc
        for (q=0; q<in->Q; q++) {
            if (temp->Xcq[c][q] == 1) {
                for (i=0; i<in->I; i++) {
                    if (in->Eci[c][i] == 1 && temp->Zi == 0) {
                        p1++;
                    }
                }
            }
        }

    }

    return p1;
}

int check2(Sol *temp, Instances *in)
// check constraint 2
// returns 0 if every query has at most one configuration associated, otherwise it returns a value grater then 0, scaled by a factor.
{
    int q, c;
    int p2, p2max=0;
    //TODO: aggiungere fattore di scalamento

    for(q=0; q<in->Q; q++)
    {
        p2=0;
        for(c=0; c<in->C; c++)
        {
            p2+=temp->Xcq[c][q];
        }
        if(p2>p2max)
        {
            p2max=p2;
        }
    }

    return p2max - 1;
}

int check3(Sol *temp, Instances *in)
// check constraint 3
// returns 0 if limit M is not exceeded, otherwise it returns the excess amount, scaled by a factor.
{
    //TODO: aggiungere fattore di scalamento

    if(temp->mem > in->M)
    {
        return temp->mem - in->M;
    }

    return 0;
}
*/
