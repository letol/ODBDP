#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>

#define true 1
#define false 0

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
int check1(Sol *temp, Instances *in);
int check2(Sol *temp, Instances *in);
int check3(Sol *temp, Instances *in);
void createZi(Sol* temp, Instances *in);

int main(int argc, char* argv[])
{
    Instances in;
    Sol best, temp;
    FILE *fin, *fout;
    time_t start=time(NULL);
    int timelimit=0;

    assert(argc == 4);

    assert(strcmp(argv[2], "-t") == 0);

    timelimit = atoi(argv[3]);

    fin = fopen(argv[1], "r");

    assert(fin != NULL);

    initialization(fin, &in, &best, &temp);


    while ((time(NULL) - start)<timelimit){

    }
    /*
        while(clock()<30000){}
    */

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

//int check(Sol *temp, Instances *in)
//// Returns 1 if temp is feasible, 0 otherwise.
//{
//    int c, q;
//    int con2=0;
//
//    //Constraint (2)
//    for(q=0; q<in->Q; q++)
//    {
//        con2=0;
//        for(c=0; c<in->C; c++)
//        {
//            con2+=temp->Xcq[c][q];
//        }
//        if(con2>1)
//        {
//            return 0;
//        }
//    }
//
//    //Constraint (3)
//    if(temp->mem > in->M)
//    {
//        return 0;
//    }
//    return 1;
//}

int check(Sol *temp, Instances *in)
// Returns 1 if temp is feasible, 0 otherwise.
{
    if (check1(temp,in) > 0 ||
        check2(temp,in) > 0 ||
        check3(temp,in) > 0) {
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

void createZi(Sol* temp, Instances *in)
{
    int i, c, q;

    for(c=0; c<in->C; c++)
    {
        for(q=0; q<in->Q; q++)
        {
            for(i=0; i<in->I; i++)
            {
                temp->Zi[i]=0;
            }
        }
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
