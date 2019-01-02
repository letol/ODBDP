#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>

#define true 1
#define false 0
#define numN 5

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

typedef struct vett{
    int *vettX;
    float *vettR;
    int mem;
    int gain;
    int *Zi;
    int feasible;
    int *Yc;
    int *Fc;
    int *confmem;
    int *gainc;
    int availableMem;
}Vett;

void initialization(FILE *fin, Sol *best, Vett *pop, Instances *in);
void letturavet(int *v, FILE *fin, int r);
void letturamat(int **m, FILE *fin, int r, int c);


void TSinit(Vett *pop, Instances *in);
void createZi(Vett *v, Instances *in);
void calculateOF(Vett *v, Instances *in);
int check3(Vett *v, Instances *in);
void changeConfMem(Vett *x, Instances *in);

void Vett2Sol(Vett *v, Sol *s, Instances *in);
void changeFc(Vett *x, Instances *in);
void cpyVettX(Vett *v1, Vett *v2, Instances *in);
void createYc(Vett *temp, Instances *in);
void createVettXfromYc(Vett *x, Instances *in);
void changeGainC(Vett *x, Instances *in);
void configureC(Vett *x, Instances *in);
int evaluateGain(int c, Vett *temp, Instances *in);

void createZifromYc(Vett *v, Instances *in);

int main(int argc, char* argv[])
{
    int c, q, i, iteration=0;
    int count, actualMemory, gainC;
    Instances in;
    Sol best;
    FILE *fin, *fout;
    time_t start=time(NULL);
    int timelimit=0;
    Vett temp, *nei;
    int worst, worstVal, bestq, bestc;
    int bestVal;
    int *changed;
    int worstmem=0;
    int fc=0;

    assert(argc == 4);

    assert(strcmp(argv[2], "-t") == 0);

    timelimit = atoi(argv[3]);

    fin = fopen(argv[1], "r");

    assert(fin != NULL);

    initialization(fin, &best, &temp, &in);

    TSinit(&temp, &in);

    assert( (nei = malloc(numN*sizeof(Vett))) != NULL);

    for(c=0; c<numN; c++)
    {
        assert( (nei[c].vettX = malloc(in.Q*sizeof(int))) != NULL);
        assert( (nei[c].vettR = malloc(in.Q*sizeof(float))) != NULL);
        assert( (nei[c].Zi = malloc(in.I*sizeof(int))) != NULL);
    }


    assert( (changed=malloc(in.C*sizeof(int)))!= NULL );
/*
    printf("\navailablemem: %d\n", temp.availableMem);



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

    printf("\nvettore Yc:\n ");

    for(c=0; c<in.C; c++)
    {
        printf("%d ", temp.Yc[c]);
    }

    printf("\nvettore confmem:\n ");

    for(c=0; c<in.C; c++)
    {
        printf("%d ", temp.confmem[c]);
    }

    printf("\nvettore Fc:\n ");

    for(c=0; c<in.C; c++)
    {
        printf("%d ", temp.Fc[c]);
    }

    printf("\nvettore gainc:\n ");

    for(c=0; c<in.C; c++)
    {
        printf("%d ", temp.gainc[c]);
    }
*/
    Vett2Sol(&temp, &best, &in);

    for(q=0; q<in.C; q++)
    {
        changed[q]=0;
    }

    count=0;



    while ((time(NULL) - start)<=timelimit)
    {
        worst=0;
        worstVal=10000;
        bestq=0;
        bestc=0;
        bestVal=0;
        worstmem=0;
        count++;

        for(q=0;q<in.C;q++)
        {
            if(changed[q]==count)
            {
                changed[q]=0;
            }
        }

        if(temp.feasible==1)
        {
            for(c=0; c<in.C; c++)
            {
                if(temp.Yc[c]==1 && changed[c]==0)
                {
                    for(fc=0,i=0; i<in.I; i++)
                    {
                        if(in.Eci[c][i]==1)
                        {
                            fc+=in.Fi[i];
                        }
                    }
                    if(worstVal > (temp.gainc[c]-fc))
                    {
                        worst=c;
                        worstVal=temp.gainc[c]-fc;
                    }
                }
                else if(temp.Yc[c]==0 && bestVal < (bestq=evaluateGain(c,&temp,&in)-temp.Fc[c]) && changed[c]==0)
                {
                    bestVal=bestq;
                    bestc=c;
                }
            }
            //    printf(" fea %d %d\n", worst, bestc);
            temp.Yc[worst]=0;
            temp.Yc[bestc]=1;
            changed[worst]=count;
            changed[bestc]=count;
        }
        else
        {
            /*
            for(i=0; i<in.I; i++)
            {
                if(bestVal < in.Mi[i] && temp.Zi[i]==1)
                {
                    bestc=i;
                    bestVal = in.Mi[i];
                }
            }
            temp.Zi[bestc]=0;
            for(c=0; c<in.C;c++)
            {
                if(in.Eci[c][bestc]==1)
                {
                    temp.Yc[c]=0;
                }
            }
            */

            for(c=0; c<in.C; c++)
            {
                if(temp.Yc[c]==1 && changed[c]==0)
                {
                    for(bestVal=0, i=0; i<in.I; i++)
                    {
                        if(in.Eci[c][i]==1)
                        {
                            bestVal+=in.Mi[i];
                        }
                    }
                    if(worstmem < bestVal)
                    {
                        worstmem = bestVal;
                        worst=c;
                    }
                }
                else if(temp.Yc[c]==0 && changed[c]==0)
                {
                    if(worstVal>temp.confmem[c])
                    {
                        bestc=c;
                        worstVal=temp.confmem[c];
                    }
                }
            }
            temp.Yc[worst]=0;
            temp.Yc[bestc]=1;
            changed[worst]=count;
            changed[bestc]=count;
            //  printf(" no fea %d %d\n", worst, bestc);
        }

        createVettXfromYc(&temp, &in);
        createYc(&temp, &in);
        configureC(&temp, &in);
        createZi(&temp, &in);
        calculateOF(&temp, &in);
        check3(&temp, &in);

        if( (iteration%600) == 599)
        {
            TSinit(&temp, &in);
            iteration=0;
        }

        if(temp.gain > best.gain && temp.feasible==1)
        {
            Vett2Sol(&temp, &best, &in);
            iteration=0;
            count=0;
            for(q=0;q<in.C;q++)
            {
                changed[q]=0;
            }

        }
        else
        {
            iteration++;
        }

        if(count==7)
        {
            count=0;
        }

        //printf("gain: %d feasible: %d\n", temp.gain, temp.feasible);


    }
    /*
        while(clock()<30000){}
    */


    printf("\nBest solution is:\n");
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

void initialization(FILE *fin, Sol *best, Vett *pop, Instances *in)
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
    assert( (pop->Yc=malloc(in->C*sizeof(int))) != NULL);
    assert( (pop->Fc=malloc(in->C*sizeof(int))) != NULL );
    assert( (pop->confmem=malloc(in->C*sizeof(int))) != NULL );
    assert( (pop->gainc=malloc(in->C*sizeof(int))) != NULL );

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
    int Memory, numconf; //memory tiene il conto della memoria attualmente disponibile e permette di verificare che una configurazione non causi infeasibility
    pop->feasible=-1;
    pop->gain=-1;


    while(pop->feasible==0 || pop->gain <=0)
    {
        count=0;
        Memory=in->M;
        numconf=0;

        for(c=0; c<in->C; c++)
        {
            pop->confmem[c]=0;
            for(i=0; i<in->I; i++)
            {
                if(in->Eci[c][i]==1)
                {
                    pop->confmem[c]+=in->Mi[i]; //init confmem
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
            if( (Memory - pop->confmem[c]) > 0) //se non supero la memoria disponibile.
            {
                numconf++;
                for(q=0;q<in->Q;q++)
                {
                    if((pop->vettX[q]) == -1 && (in->Gcq[c][q]) > 0 ) //assegno la c a q se fornisce guadagno e la query non � servita
                    {
                        pop->vettX[q]=c;
                    }
                    else if( (pop->vettX[q]) > -1)
                    {
                        if(in->Gcq[c][q] > in->Gcq[pop->vettX[q]][q])   //se la q � servita da una c che fornisce un guadagno minore, sostituiscila
                        {
                            pop->vettX[q]=c;
                        }
                    }
                }

                Memory-=pop->confmem[c];
                createZi(pop, in);
                changeConfMem(pop, in);
                calculateOF(pop, in);
                check3(pop, in);
                createYc(pop, in);
                changeFc(pop, in);
                changeGainC(pop, in);
            }
        }
    }

    pop->availableMem=in->M - pop->mem;

    for(c=0, numconf=0; c<in->C; c++)
    {
        if(pop->confmem[c] <= pop->availableMem)
        {
            pop->Yc[c]=1;
            pop->availableMem-=pop->confmem[c];
            pop->confmem[c]=0;
            numconf++;
        }
    }

    createVettXfromYc(pop, in);
    createYc(pop, in);
    createZifromYc(pop, in);
    configureC(pop, in);
    calculateOF(pop, in);
    check3(pop, in);
    printf("configurazioni usate %d\n", numconf);

    pop->availableMem=in->M - pop->mem;
    printf("confmem %d\n", pop->availableMem);

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
void createYc(Vett *temp, Instances *in)
{
    int q, c;

    for(c=0; c<in->C; c++)
    {
        temp->Yc[c]=0;
    }

    for(q=0; q<in->Q; q++)
    {
        if(temp->vettX[q] > -1)
        {
            temp->Yc[temp->vettX[q]]=1;
        }
    }
}

void changeConfMem(Vett *x, Instances *in)
{
    int i, j;

    for(j=0; j<in->C; j++)
    {
        x->confmem[j]=0;
        for(i=0; i<in->I; i++)
        {
            if(in->Eci[j][i]==1)
            {
                x->confmem[j]+=in->Mi[i];
            }
        }
    }

    for(i=0; i<in->I; i++)
    {
        for(j=0; j<in->C; j++)
        {
            if(in->Eci[j][i] == 1 && x->Zi[i]==1)
            {
                x->confmem[j]-=in->Mi[i];
            }
        }

    }
}

void changeFc(Vett *x, Instances *in)
{
    int i, c;

    for(c=0; c<in->C; c++)
    {
        x->Fc[c]=0;
        for(i=0; i<in->I; i++)
        {
            if(in->Eci[c][i]==1)
            {
                x->Fc[c]+=in->Fi[i];
            }
        }
    }

    for(i=0; i<in->I; i++)
    {
        for(c=0; c<in->C; c++)
        {
            if(in->Eci[c][i] == 1 && x->Zi[i]==1)
            {
                x->Fc[c]-=in->Fi[i];
            }
        }

    }
}



void changeGainC(Vett *x, Instances *in)
{
    int q,c;

    for(c=0;c<in->C; c++)
    {
        x->gainc[c]=0;
    }
    for(q=0; q<in->Q; q++)
    {
        if(x->vettX[q]>-1)
        {
            x->gainc[x->vettX[q]]+=in->Gcq[x->vettX[q]][q];
        }
    }

}

void configureC(Vett *x, Instances *in)
{
    changeFc(x,in);
    changeConfMem(x,in);
    changeGainC(x,in);
}

void createVettXfromYc(Vett *x, Instances *in)
{
    int c,q;

    for(q=0; q<in->Q; q++)
    {
        x->vettX[q]=-1;
    }

    for(c=0;c<in->C; c++)
    {
        for(q=0; q<in->Q; q++)
        {
            if(x->vettX[q]!=-1)
            {
                if(x->Yc[c]==1 && in->Gcq[c][q] > in->Gcq[x->vettX[q]][q])
                {
                    x->vettX[q]=c;
                }
            }
            else if( x->vettX[q]==-1 && x->Yc[c]==1 && in->Gcq[c][q]>0)
            {
                x->vettX[q]=c;
            }
        }
    }
    return;
}

void cpyVettX(Vett *v1, Vett *v2, Instances *in)
{
    int q;

    for(q=0; q<in->Q; q++)
    {
        v1->vettX[q]=v2->vettX[q];
    }
    for(q=0; q<in->C; q++)
    {
        v1->Yc[q]=v2->Yc[q];
    }
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

int evaluateGain(int c, Vett *temp, Instances *in)
{
    int q;
    int gain=0, mom=0;

    for(q=0;q<in->Q; q++)
    {
        if(temp->vettX[q]==-1 && in->Gcq[c][q]>0)
        {
            gain+=in->Gcq[c][q];
        }
        else if( temp->vettX[q] > -1)
        {
            if( (mom=in->Gcq[c][q] - in->Gcq[temp->vettX[q]][q]) > 0)
            {
                gain+=mom;
            }
        }
    }
    return gain;
}

void createZifromYc(Vett *v, Instances *in)
{
    int i, c;

    for(i=0; i<in->I; i++)
    {
        v->Zi[i]=0;
    }

    for(c=0; c<in->C; c++)
    {
        if(v->Yc[c]==1)
        {
            for(i=0; i<in->I; i++)
            {
                if(in->Eci[c][i]==1)
                {
                    v->Zi[i]=1;
                }
            }
        }

    }

    return;
}
