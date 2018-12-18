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

}Sol;

void initialization(FILE *fin, Instances *in, Sol *best, Sol *temp);
void letturavet(int *v, FILE *fin, int r);
void letturamat(int **m, FILE *fin, int r, int c);
void calculateOF(Sol *temp, Instances *in);
int check(Sol *temp, Instances *in);
int check1(Sol *temp, Instances *in);
int check2(Sol *temp, Instances *in); //TODO
int check3(Sol *temp, Instances *in); //TODO
Sol solGen(Sol *temp, Instances *in); //TODO

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
    int i, j, gain=0, cost=0;

    for (i=0; i<in->C; i++) {
        for (j=0; j<in->Q; j++) {
            if (temp->Xcq[i][j] == 1) {
                gain += in->Gcq[i][j];
            }
        }
    }

    for (i=0; i<in->I; i++) {
        if (temp->Zi[i] == 1) {
            cost += in->Fi[i];
        }
    }

    temp->gain = (gain-cost);
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

int check1(Sol *temp, Instances *in)
// check constraint 1
{
    int c, q, i;
    int p1=0;

    for (c=0; c<in->C; c++) {
        for (q=0; q<in->Q; q++) {
            if (temp->Xcq[c][q] == 1) {
                for (i=0; i<in->I; i++) {
                    if (in->Eci[c][i] == 1 && temp->Zi == 0) {
                        pi++;
                    }
                }
            }
        }

    }

    return p1;
}

int check2(Sol *temp, Instances *in)
// check constraint 2
{

}