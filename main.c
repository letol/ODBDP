#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <math.h>
#include <windows.h>

#define debug 1

typedef struct instances{
    int Q, I, C, M;
    int **Eci;
    int *Fi, *Mi;
    int **Gcq;
    char *name;
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
    float *vettC;
    float *vettR;
    float fitness;
    int mem;
    int *confmem;
    int availableMem;
    int gain;
    int *Zi;
    int feasible;
    int *Yc;
    int *Fc;
    int *eventualMemc;
    int *actualMemc;
    int *eventualGainc;
    int *actualGainc;
}Vett;

int numP=50;
float probOffset=.01;
int TSiter=50;
int numN=5;

Sol best;

void mainInitialization(FILE *fin, Instances *in, Sol *best);
void letturavet(int *v, FILE *fin, int r);
void letturamat(int **m, FILE *fin, int r, int c);
void bestWrite(Sol *best, Instances *in);
char *GetFileName(const char *path);


void TSinit(Vett *pop, Instances *in);
void TSdataInit(Vett *neighbor, Vett *temp, Sol *solu, Sol *TSbest, Instances *in);
void createZi(Vett *v, Instances *in);
int check3(Vett *v, Instances *in);

void eventualMemCalc(Vett *x, Instances *in);
void changeFc(Vett *x, Instances *in);
void changeGainC(Vett *x, Instances *in);
void actualMemCalc(Vett *temp, Instances *in);
void configureC(Vett *x, Instances *in);

void calculateOF(Vett *v, Instances *in);

void createYc(Vett *temp, Instances *in);
void createVettXfromYc(Vett *x, Instances *in);
void createZifromYc(Vett *v, Instances *in);

void eventualGainCalc(Vett *temp, Instances *in);
int evaluateGain(int c, Vett *temp, Instances *in);

void Vett2Sol(Vett *v, Sol *s, Instances *in);
void createSolution(Vett *x, Instances *in);
void cpyVett(Vett *v1, Vett *v2, Instances *in);


void GAdataInit(Instances *in, Vett *pop, Vett *temp, Vett *selectedParents, Vett *children, Vett *soluz,
                Vett *neighbor, int **chosen, int **neibest, int **neiworst);
void GAinit(Vett *pop, Instances *in);
void searchMax(Vett *pop, int dim, Sol *best, Instances *in, time_t time);
void calculateFitness(Vett *v, Instances *in);
void mutation(Vett *v, Instances *in);
void selectParents(Vett *pop, Vett *selPar, Instances *in);
void crossover(Vett *p1, Vett *p2, Vett *c1, Vett *c2, Instances *in);
void sobstitute(Vett *pop, Vett *children, Instances *in);
void calculateC(Vett *v, Instances *in);
void calculateR(Vett *v, Instances *in);
void cpySol(Vett *dst, Vett *src, Instances *in);
void changeConfMem(Vett *x, Instances *in);
void TS(Vett *localBest, Vett *temp, Vett *soluz, Vett *neighbor, int *chosen, int *neibest, int *neiworst,
        Instances *in);


DWORD WINAPI Hybrid(Instances *in)
{
    int j, i, Xcount=0;

    Vett population[numP], temp;
    Vett selectedParents[numP/2], children[numP/2];

    Vett soluz, neighbor;
    int *chosen, *neibest, *neiworst;

    time_t start=time(NULL);

#if debug == 0
    srand((unsigned int)time(NULL));
#endif

    GAdataInit(in, population, &temp, selectedParents, children,
               &soluz, &neighbor, &chosen, &neibest, &neiworst);

    GAinit(population, in);

#if debug >= 2
    printf("HYB: INITIAL POPULATION\n\n ");
    for(j=0; j<numP; j++)
    {
        printf("\nHYB: ");
        for(q=0; q<in.Q; q++)
        {
            printf("%d ", population[j].vettX[q]);
        }
        printf("\nHYB: vettore zi%d\nHYB: ", j);
        for(i=0; i<in.I; i++)
        {
            printf("%d ", population[j].Zi[i]);
        }
        printf("\nHYB: vettore c%d\nHYB: ", j);
        for(i=0; i<in.I; i++)
        {
            printf("%f ", population[j].vettC[i]);
        }
        printf("\nHYB: vettore r%d\nHYB: ", j);
        for(i=0; i<in.I; i++)
        {
            printf("%f ", population[j].vettR[i]);
        }
        printf("\nHYB: gain: %d\nHYB: memory: %d\nHYB: feasible:%d \n", population[j].gain, population[j].mem,population[j].feasible);
    }
#endif

    searchMax(population, numP, &best, in, time(NULL)-start);

    for(i=0; i<numP; i++) {
        calculateFitness(&population[i], in);
        mutation(&population[i], in);
    }

    while (1) {
        selectParents(population, selectedParents, in);
#if debug >= 2
        printf("\nHYB: SELECTED PARENTS\n");
        for (j = 0; j < numP / 2; j++) {
            if (selectedParents[j].feasible != -2) {
                printf("\nHYB: P%d\nHYB: ", j);
                for (q = 0; q < in->Q; q++) {
                    printf("%d ", selectedParents[j].vettX[q]);
                }
                printf("\nHYB: vettore zi%d\nHYB: ", j);
                for (i = 0; i < in->I; i++) {
                    printf("%d ", selectedParents[j].Zi[i]);
                }
                printf("\nHYB: vettore c%d\nHYB: ", j);
                for (i = 0; i < in->Q; i++) {
                    printf("%f ", selectedParents[j].vettC[i]);
                }
                printf("\nHYB: vettore r%d\nHYB: ", j);
                for (i = 0; i < in->Q; i++) {
                    printf("%f ", selectedParents[j].vettR[i]);
                }
                printf("\nHYB: gain: %d\nHYB: memory: %d\nHYB: feasible:%d \n", selectedParents[j].gain, selectedParents[j].mem,
                       selectedParents[j].feasible);
            }
        }
#endif

#if debug >= 2
        printf("\nHYB: CROSSOVER & MUTATION\n");
#endif
        for (i = 0; i < numP / 2; i++) {
            children[i].feasible = -2;
        }
        i = 0;
        while (i < numP / 2) {
            if (selectedParents[i].feasible != -2) {
                j = i + 1;
                while (j < numP / 2) {
                    if (selectedParents[j].feasible != -2) {
                        crossover(&selectedParents[i], &selectedParents[j], &children[i], &children[j], in);
                        Xcount++;

                        createYc(&children[i], in);
                        createZi(&children[i], in);
                        calculateOF(&children[i], in);
                        check3(&children[i], in);
                        configureC(&children[i], in);
                        eventualGainCalc(&children[i], in);

                        TS(&children[i], &temp, &soluz, &neighbor, chosen, neibest, neiworst, in);

                        calculateC(&children[i], in);
                        calculateR(&children[i], in);
                        calculateFitness(&children[i], in);


                        createYc(&children[j], in);
                        createZi(&children[j], in);
                        calculateOF(&children[j], in);
                        check3(&children[j], in);
                        configureC(&children[j], in);
                        eventualGainCalc(&children[j], in);

                        TS(&children[j], &temp, &soluz, &neighbor, chosen, neibest, neiworst, in);

                        calculateC(&children[j], in);
                        calculateR(&children[j], in);
                        calculateFitness(&children[j], in);
                    }
                    j++;
                }
            }
            i++;
        }

#if debug >= 2
        printf("\nHYB: PRODUCED CHILDREN\n");
        for (j = 0; j < numP / 2; j++) {
            if (children[j].feasible != -2) {
                printf("\nHYB: C%d\nHYB: ", j);
                for (q = 0; q < in->Q; q++) {
                    printf("%d ", children[j].vettX[q]);
                }
                printf("\nHYB: vettore zi%d\nHYB: ", j);
                for (i = 0; i < in->I; i++) {
                    printf("%d ", children[j].Zi[i]);
                }
                printf("\nHYB: vettore c%d\nHYB: ", j);
                for (i = 0; i < in->Q; i++) {
                    printf("%f ", children[j].vettC[i]);
                }
                printf("\nHYB: vettore r%d\nHYB: ", j);
                for (i = 0; i < in->Q; i++) {
                    printf("%f ", children[j].vettR[i]);
                }
                printf("\nHYB: gain: %d\nHYB: memory: %d\nHYB: feasible:%d \n", children[j].gain, children[j].mem,
                       children[j].feasible);
            }
        }
#endif

        searchMax(children, numP / 2, &best, in, time(NULL)-start);

        sobstitute(population, children, in);

        if (Xcount == 0) {
#if debug >= 2
            printf("\nHYB: NO CROSSOVER -> ONLY MUTATION\n");
#endif
            for (i = 0; i < numP; i++) {
                mutation(&population[i], in);
            }
#if debug >= 2
            for(j=0; j<numP; j++)
            {
                printf("\nHYB: ");
                for(q=0; q<in.Q; q++)
                {
                    printf("%d ", population[j].vettX[q]);
                }
                printf("\nHYB: vettore zi%d\nHYB: ", j);
                for(i=0; i<in.I; i++)
                {
                    printf("%d ", population[j].Zi[i]);
                }
                printf("\nHYB: vettore c%d\nHYB: ", j);
                for(i=0; i<in.Q; i++)
                {
                    printf("%f ", population[j].vettC[i]);
                }
                printf("\nHYB: vettore r%d\nHYB: ", j);
                for(i=0; i<in.Q; i++)
                {
                    printf("%f ", population[j].vettR[i]);
                }
                printf("\nHYB: gain: %d\nHYB: memory: %d\nHYB:feasible:%d \n", population[j].gain, population[j].mem,population[j].feasible);
            }
#endif
        }

        Xcount = 0;

        searchMax(population, numP, &best, in, 0);

    }

    return 0;
}

DWORD WINAPI TabuSearch(Instances *in)
{
    int c, q, i, iteration=0,j;
    int count;
    Vett temp, neighbor;
    Sol soluz, TSbest;
    int bestc, worst;
    int bestVal, worstVal;
    int nof=0, fea=0;
    int flag=0;
    int *chosen, *neibest, *neiworst;
    int tabu[2][8];

    time_t start=time(NULL);

#if debug == 0
    srand((unsigned int)time(NULL));
#endif

    TSdataInit(&neighbor, &temp, &soluz, &TSbest, in);
    assert( (neibest = malloc(numN*sizeof(int))) != NULL);
    assert( (neiworst = malloc(numN*sizeof(int))) != NULL);
    assert( (chosen = malloc(in->C*sizeof(int))) != NULL);

    TSinit(&temp, in);

    iteration=0;

    count=0;

    Vett2Sol(&temp, &TSbest, in);

    soluz.gain=0;

    while (1)
    {
        count++;

        for(i=0; i<numN; i++)
        {
            neibest[i]=-1;
            neiworst[i]=-1;
        }

        if(temp.feasible==1)
        {

            for(c=0; c<in->C; c++)
            {
                chosen[c]=0;
            }
            for(i=0; i<numN; i++)
            {
                bestc=-1;
                bestVal=-1000000000;
                worst=-1;
                worstVal=1000000000;

                for(c=0; c<in->C; c++)
                {
                    if(temp.Yc[c]==0 && chosen[c]==0 && bestVal < temp.eventualGainc[c])
                    {
                        bestVal = temp.eventualGainc[c];
                        bestc=c;
                    }
                    else if(temp.Yc[c]==1 && chosen[c]==0 && worstVal > temp.actualGainc[c])
                    {
                        worstVal=temp.actualGainc[c];
                        worst=c;
                    }
                }

                neibest[i]=bestc;
                neiworst[i]=worst;
                chosen[bestc]=1;
                chosen[worst]=1;
            }

            for(c=0;c<in->C; c++)
            {
                neighbor.Yc[c]=temp.Yc[c];
            }

            bestVal=-1000000;

            for(i=0; i<numN; i++)
            {
                if(neibest[i]!=-1)
                {
                    for(j=0; j<numN; j++)
                    {
                        flag=1;
                        if(neiworst[j]!=-1)
                        {
                            for(c=1;c<8;c++)
                            {
                                if( (tabu[0][c]==neibest[i] && tabu[1][c]==neiworst[j]) || (tabu[0][c]==neiworst[j] && tabu[1][c]==neibest[i]) )
                                {
                                    flag=0;
                                }
                            }
                            if(flag==1)
                            {

                                neighbor.Yc[neibest[i]]=1;
                                neighbor.Yc[neiworst[j]]=0;
                                createVettXfromYc(&neighbor, in);
                                createZi(&neighbor, in);
                                calculateOF(&neighbor, in);
                                if(neighbor.gain > bestVal)
                                {
                                    bestVal=neighbor.gain;
                                    bestc=neibest[i];
                                    worst=neiworst[j];
                                }
                                neighbor.Yc[neibest[i]]=0;
                                neighbor.Yc[neiworst[j]]=1;
                            }
                        }
                    }
                }
            }

            temp.Yc[bestc]=1;
            temp.Yc[worst]=0;

            tabu[0][count]=bestc;
            tabu[1][count]=worst;
            nof=0;
            fea++;

        }
        else
        {
            while(temp.feasible==0)
            {
                for(c=0; c<in->C; c++)
                {
                    chosen[c]=0;
                }

                for(i=0; i<numN; i++)
                {
                    bestVal=10000000;
                    bestc=-1;
                    worstVal=-1;
                    worst=-1;

                    for(c=0; c<in->C; c++)
                    {
                        if(temp.Yc[c]==0 && chosen[c]==0 && bestVal > temp.eventualMemc[c] && temp.eventualMemc!=0)
                        {
                            bestVal = temp.eventualMemc[c];
                            bestc=c;
                        }
                        else if(temp.Yc[c]==1 && chosen[c]==0 && worstVal < temp.actualMemc[c])
                        {
                            worstVal=temp.actualMemc[c];
                            worst=c;
                        }
                    }

                    neibest[i]=bestc;
                    neiworst[i]=worst;
                    chosen[bestc]=1;
                    chosen[worst]=1;
                }


                for(c=0;c<in->C; c++)
                {
                    neighbor.Yc[c]=temp.Yc[c];
                }

                bestVal=-10000000;

                for(i=0; i<numN; i++)
                {
                    if(neibest[i]!=-1)
                    {
                        for(j=0; j<numN; j++)
                        {
                            flag=1;
                            if(neiworst[j]!=-1)
                            {
                                for(c=1;c<=7;c++)
                                {
                                    if( (tabu[0][c]==neibest[i] && tabu[1][c]==neiworst[j]) || (tabu[0][c]==neiworst[j] && tabu[1][c]==neibest[i]) )
                                    {
                                        flag=0;
                                    }
                                }
                                if(flag==1)
                                {
                                    neighbor.Yc[neibest[i]]=1;
                                    neighbor.Yc[neiworst[j]]=0;
                                    createVettXfromYc(&neighbor, in);
                                    createZi(&neighbor, in);
                                    calculateOF(&neighbor, in);
                                    if(neighbor.gain > bestVal)
                                    {
                                        bestVal=neighbor.gain;
                                        bestc=neibest[i];
                                        worst=neiworst[j];
                                    }
                                    neighbor.Yc[neibest[i]]=0;
                                    neighbor.Yc[neiworst[j]]=1;
                                }
                            }
                        }
                    }
                }

                temp.Yc[bestc]=1;
                temp.Yc[worst]=0;

                tabu[0][count]=bestc;
                tabu[1][count]=worst;
                nof++;
                fea=0;

                createZifromYc(&temp, in);
                check3(&temp, in);

                if(nof>=20)
                {
                    nof=0;
                    while(temp.feasible==0)
                    {
                        worstVal=-1;
                        for(c=rand()%in->C; temp.Fc[c]==1; c=rand()%in->C);

                        if(temp.Yc[c]==1)
                        {
                            temp.Yc[c]=0;
                        }
                        nof=0;
                        createZifromYc(&temp, in);
                        check3(&temp, in);
                    }

                }
            }
        }

        createSolution(&temp,in);

        if(fea>=5)
        {

            for(q=rand()%in->Q, flag=0; temp.vettX[q]==-1 && flag<in->Q; q=rand()%in->Q, flag++);

            for(c=rand()%in->C, flag=0; flag<in->C; c=rand()%in->C, flag++)
            {
                if(in->Gcq[c][q]>0 && temp.Yc[c]==0)
                {
                    temp.Yc[c]=1;
                    flag=in->C;
                }
            }
            createSolution(&temp,in);
            fea=0;
        }


        if( (iteration%50) == 49)
        {
            TSinit(&temp, in);
            iteration=0;
            count=0;
        }

        if(temp.gain > TSbest.gain && temp.feasible==1)
        {
            Vett2Sol(&temp, &TSbest, in);
            cpyVett(&neighbor, &temp, in);

            bestc=0;

            while(bestc>-1)
            {
                bestc=-1;
                bestVal=0;
                for(c=0;c<in->C; c++)
                {
                    if(neighbor.eventualGainc[c] > bestVal && neighbor.eventualMemc[c] < (in->M - neighbor.mem))
                    {
                        bestVal=neighbor.eventualGainc[c];
                        bestc=c;
                    }
                    if(bestc>-1)
                    {
                        neighbor.Yc[bestc]=1;
                        createSolution(&neighbor,in);
                    }
                }
            }
            if(neighbor.gain > soluz.gain)
            {
                Vett2Sol(&neighbor, &soluz, in);
            }
            if(soluz.gain > TSbest.gain)
            {
                Vett2Sol(&neighbor, &TSbest, in);
            }
            if(TSbest.gain > best.gain)
            {
                Vett2Sol(&neighbor, &best, in);
                bestWrite(&TSbest, in);
#if debug >= 1
                printf("\nTS: gain: %d\nTS: memory: %d\nTS: feasible: %d\nTS: time: %ld\n", soluz.gain, soluz.mem, soluz.feasible, time(NULL)-start);
#endif
            }

            iteration=0;
        }
        else
        {
            iteration++;
        }

        if(count==7)
        {
            count=0;
        }
    }

    return 0;
}

int main(int argc, char*argv[]) {
    int t=atoi(argv[3]);
    DWORD timelimit=(DWORD)t*1000;
    Instances in;
    FILE *fin;

    fin = fopen(argv[1], "r");
    assert(fin != NULL);
    in.name = GetFileName(argv[1]);
    mainInitialization(fin, &in, &best);

    fclose(fin);

    best.gain=0;
    HANDLE thread1 = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)Hybrid, &in, 0, NULL);

    //HANDLE thread2 = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)TabuSearch, &in, 0, NULL );

    Sleep(timelimit);
    TerminateThread(thread1, 0);
    //TerminateThread(thread2, 0);

#if debug >= 1
    printf("\nBest solution is:\n");
    for(int c=0; c<in.C; c++)
    {
        for(int q=0; q<in.Q; q++)
        {
            printf("%d ", best.Xcq[c][q]);
        }
        printf("\n");
    }
    printf("\n");
    for(int i=0; i<in.I; i++)
    {
        printf("%d ", best.Zi[i]);
    }

    printf("\ngain: %d\nmemory: %d\nfeasible: %d", best.gain, best.mem, best.feasible);
#endif

    return 0;
}

void mainInitialization(FILE *fin, Instances *in, Sol *best) {
    char lecture[40];
    int c = 0;

    assert((fscanf(fin, "%s %d", lecture, &in->Q)) != EOF);
    assert((fscanf(fin, "%s %d", lecture, &in->I)) != EOF);
    assert((fscanf(fin, "%s %d", lecture, &in->C)) != EOF);
    assert((fscanf(fin, "%s %d", lecture, &in->M)) != EOF);

    assert((in->Eci = malloc(in->C * sizeof(int *))) != NULL);

    for (c = 0; c < in->C; c++) {
        assert((in->Eci[c] = malloc(in->I * sizeof(int))) != NULL);
    }

    assert((in->Fi = malloc(in->I * sizeof(int))) != NULL);
    assert((in->Mi = malloc(in->I * sizeof(int))) != NULL);

    assert((in->Gcq = malloc(in->C * sizeof(int *))) != NULL);

    for (c = 0; c < in->C; c++) {
        assert((in->Gcq[c] = malloc(in->Q * sizeof(int))) != NULL);
    }

    assert((best->Xcq = malloc(in->C * sizeof(int *))) != NULL);

    for (c = 0; c < in->C; c++) {
        assert((best->Xcq[c] = malloc(in->Q * sizeof(int))) != NULL);
    }

    assert((best->Zi = malloc(in->I * sizeof(int))) != NULL);

    assert((fscanf(fin, "%s", lecture)) != EOF);

    letturamat(in->Eci, fin, in->C, in->I);

    assert((fscanf(fin, "%s", lecture)) != EOF);

    letturavet(in->Fi, fin, in->I);

    assert((fscanf(fin, "%s", lecture)) != EOF);

    letturavet(in->Mi, fin, in->I);

    assert((fscanf(fin, "%s", lecture)) != EOF);

    letturamat(in->Gcq, fin, in->C, in->Q);
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

void TSinit(Vett *pop, Instances *in)
{
    int i, c, count;
    int Memory=0;
    pop->feasible=-1;
    pop->gain=-1;


    while(pop->feasible==0 || pop->gain <=0)
    {
        count=0;
        Memory=in->M;

        for(c=0; c<in->C; c++)
        {
            pop->eventualMemc[c]=0;

            for(i=0; i<in->I; i++)
            {
                if(in->Eci[c][i]==1)
                {
                    pop->eventualMemc[c]+=in->Mi[i];
                }
            }
        }

        for(c=0; c<in->C;c++)
        {
            pop->Yc[c]=0;
        }

        for(c=rand()%in->C; count<(in->C/2); count++, c=rand()%in->C)
        {
            if( (Memory - pop->eventualMemc[c]) > 0)
            {
                pop->Yc[c]=1;

                Memory-=pop->eventualMemc[c];
                createZifromYc(pop, in);
                eventualMemCalc(pop, in);
            }
        }
        createSolution(pop, in);
    }

    for(c=0; c<in->C; c++)
    {
        if(pop->eventualMemc[c] <= Memory)
        {
            pop->Yc[c]=1;
            Memory-=pop->eventualMemc[c];
            createZifromYc(pop, in);
            eventualMemCalc(pop, in);
        }
    }

    createSolution(pop, in);

}

void TSdataInit(Vett *neighbor, Vett *temp, Sol *solu, Sol *TSbest, Instances *in) {
    int c;

    assert( (temp->vettX = malloc(in->Q*sizeof(int))) != NULL);
    assert( (temp->Zi = malloc(in->I*sizeof(int))) != NULL);
    assert( (temp->Yc=malloc(in->C*sizeof(int))) != NULL);
    assert( (temp->Fc=malloc(in->C*sizeof(int))) != NULL );
    assert( (temp->eventualMemc=malloc(in->C*sizeof(int))) != NULL );
    assert( (temp->actualMemc=malloc(in->C*sizeof(int))) != NULL );
    assert( (temp->actualGainc=malloc(in->C*sizeof(int))) != NULL );
    assert( (temp->eventualGainc=malloc(in->C*sizeof(int))) != NULL );

    assert( (neighbor->vettX = malloc(in->Q*sizeof(int))) != NULL);
    assert( (neighbor->Zi = malloc(in->I*sizeof(int))) != NULL);
    assert( (neighbor->Yc=malloc(in->C*sizeof(int))) != NULL);
    assert( (neighbor->Fc=malloc(in->C*sizeof(int))) != NULL );
    assert( (neighbor->eventualMemc=malloc(in->C*sizeof(int))) != NULL );
    assert( (neighbor->actualMemc=malloc(in->C*sizeof(int))) != NULL );
    assert( (neighbor->actualGainc=malloc(in->C*sizeof(int))) != NULL );
    assert( (neighbor->eventualGainc=malloc(in->C*sizeof(int))) != NULL );

    assert( (solu->Xcq = malloc(in->C*sizeof(int*))) != NULL);

    for(c=0; c<in->C; c++)
    {
        assert( (solu->Xcq[c]=malloc(in->Q*sizeof(int))) !=NULL  );
    }

    assert( (solu->Zi = malloc(in->I*sizeof(int))) !=NULL );


    assert( (TSbest->Xcq = malloc(in->C*sizeof(int*))) != NULL);

    for(c=0; c<in->C; c++)
    {
        assert( (TSbest->Xcq[c]=malloc(in->Q*sizeof(int))) !=NULL  );
    }

    assert( (TSbest->Zi = malloc(in->I*sizeof(int))) !=NULL );

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
}

void eventualMemCalc(Vett *x, Instances *in) {

    int i, c;

    for (c = 0; c < in->C; c++)
    {
        x->eventualMemc[c] = 0;
        for (i = 0; i < in->I; i++)
        {
            if (in->Eci[c][i] == 1 && x->Zi[i] == 0)
            {
                x->eventualMemc[c] += in->Mi[i];
            }
        }
    }
}

void changeFc(Vett *x, Instances *in)
{
    int i, c;

    for (c = 0; c < in->C; c++)
    {
        x->Fc[c] = 0;
        for (i = 0; i < in->I; i++)
        {
            if (in->Eci[c][i] == 1 && x->Zi[i] == 0)
            {
                x->Fc[c] += in->Fi[i];
            }
        }
    }
}

void changeGainC(Vett *x, Instances *in)
{
    int q,c,i;
    int fc=0;

    for(c=0;c<in->C; c++)
    {
        x->actualGainc[c]=0;
    }


    for(q=0; q<in->Q; q++)
    {
        if(x->vettX[q]>-1)
        {
            x->actualGainc[x->vettX[q]]+=in->Gcq[x->vettX[q]][q];
        }
    }

    for(c=0; c<in->C; c++)
    {
        if(x->actualGainc[c]>0)
        {
            for(fc=0, i=0; i<in->I; i++)
            {
                if(in->Eci[c][i]==1)
                {
                    fc+=in->Fi[i];
                }
            }
        }
    }
}

void actualMemCalc(Vett *temp, Instances *in)
{
    int i,c;

    for(c=0; c<in->C; c++)
    {
        temp->actualMemc[c]=0;
        for(i=0; i<in->I; i++)
        {
            if(in->Eci[c][i]==1 && temp->Zi[i]==1)
            {
                temp->actualMemc[c]+=in->Mi[i];
            }
        }
    }
}
void configureC(Vett *x, Instances *in)
{

    eventualMemCalc(x,in);

    changeFc(x,in);

    changeGainC(x,in);

    actualMemCalc(x, in);

}
int check3(Vett *v, Instances *in)
{
    int i;

    v->mem=0;

    for(i=0; i<in->I; i++)
    {
        if(v->Zi[i])
        {
            v->mem+=in->Mi[i];
        }
    }

    if(v->mem > in->M)
    {
        v->feasible=0;
    }
    else
    {
        v->feasible=1;
    }


    return (in->M - v->mem);
}
void calculateOF(Vett *v, Instances *in)
{
    int q, i, gain=0, cost=0;
    v->mem=0;

    for (q=0; q<in->Q; q++)
    {
        if (v->vettX[q] >= 0)
        {
            gain += in->Gcq[v->vettX[q]][q];
        }
    }

    for (i=0; i<in->I; i++)
    {
        if (v->Zi[i] == 1)
        {
            cost += in->Fi[i];
            v->mem += in->Mi[i];
        }
    }

    v->gain = (gain-cost);

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
}

void eventualGainCalc(Vett *temp, Instances *in)
{
    int c;

    for(c=0; c<in->C; c++)
    {
        temp->eventualGainc[c]=evaluateGain(c, temp, in);
    }
}

int evaluateGain(int c, Vett *temp, Instances *in)
{
    int q;
    int gain=0;

    for(q=0;q<in->Q; q++)
    {
        if(temp->vettX[q]==-1 && in->Gcq[c][q]>0)
        {
            gain+=in->Gcq[c][q];
        }
        else if(temp->vettX[q] > -1)
        {
            if(in->Gcq[c][q] > in->Gcq[temp->vettX[q]][q])
            {
                gain+=in->Gcq[c][q] - in->Gcq[temp->vettX[q]][q];
            }
        }
    }
    return (gain- temp->Fc[c]);
}



void Vett2Sol(Vett *v, Sol *s, Instances *in)
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
}

void createSolution(Vett *x, Instances *in)
{
    createVettXfromYc(x, in);
    createYc(x, in);
    createZi(x, in);
    calculateOF(x, in);
    check3(x, in);
    configureC(x, in);
    eventualGainCalc(x, in);
}

void cpyVett(Vett *v1, Vett *v2, Instances *in)
{
    int c;

    for(c=0; c<in->C; c++)
    {
        v1->Yc[c]=v2->Yc[c];
    }
    createSolution(v1,in);
}

void bestWrite(Sol *best, Instances *in)
{
    int c,q;
    FILE *fout;
    char file[100];

    sprintf(file, "%s_OMAMZ_group04.sol", in->name);
    fout = fopen(file, "w+");

    for (c = 0; c < in->C; c++)
    {
        for (q = 0; q < in->Q; q++)
        {
            fprintf(fout,"%d", best->Xcq[c][q]);
            if(q< (in->Q - 1))
                fprintf(fout," ");
        }
        if(c< (in->C - 1))
            fprintf(fout,"\n");
    }

    fclose(fout);
}
void GAdataInit(Instances *in, Vett *pop, Vett *temp, Vett *selectedParents, Vett *children, Vett *soluz,
                Vett *neighbor, int **chosen, int **neibest, int **neiworst)
{
    int c;
    for(c=0; c<numP; c++)
    {
        assert( (pop[c].vettX = malloc(in->Q*sizeof(int))) != NULL);
        assert( (pop[c].vettC = malloc(in->Q*sizeof(float))) != NULL);
        assert( (pop[c].vettR = malloc(in->Q*sizeof(float))) != NULL);
        assert( (pop[c].Zi = malloc(in->I*sizeof(int))) != NULL);
        assert( (pop[c].Yc = malloc(in->C*sizeof(int))) != NULL);
        assert( (pop[c].Fc = malloc(in->C*sizeof(int))) != NULL);
        assert( (pop[c].eventualMemc = malloc(in->C*sizeof(int))) != NULL);
        assert( (pop[c].actualMemc = malloc(in->C*sizeof(int))) != NULL);
        assert( (pop[c].eventualGainc = malloc(in->C*sizeof(int))) != NULL);
        assert( (pop[c].actualGainc = malloc(in->C*sizeof(int))) != NULL);
        assert( (pop[c].confmem = malloc(in->C*sizeof(int))) != NULL);
    }

    for(c=0; c<numP/2; c++)
    {
        assert((children[c].vettX = malloc(in->Q * sizeof(int))) != NULL);
        assert((children[c].vettC = malloc(in->Q * sizeof(float))) != NULL);
        assert((children[c].vettR = malloc(in->Q * sizeof(float))) != NULL);
        assert((children[c].Zi = malloc(in->I * sizeof(int))) != NULL);
        assert((children[c].Yc = malloc(in->C*sizeof(int))) != NULL);
        assert((children[c].Fc = malloc(in->C*sizeof(int))) != NULL);
        assert((children[c].eventualMemc = malloc(in->C*sizeof(int))) != NULL);
        assert((children[c].actualMemc = malloc(in->C*sizeof(int))) != NULL);
        assert((children[c].eventualGainc = malloc(in->C*sizeof(int))) != NULL);
        assert((children[c].actualGainc = malloc(in->C*sizeof(int))) != NULL);
        assert((children[c].confmem = malloc(in->C*sizeof(int))) != NULL);

        assert((selectedParents[c].vettX = malloc(in->Q * sizeof(int))) != NULL);
        assert((selectedParents[c].vettC = malloc(in->Q * sizeof(float))) != NULL);
        assert((selectedParents[c].vettR = malloc(in->Q * sizeof(float))) != NULL);
        assert((selectedParents[c].Zi = malloc(in->I * sizeof(int))) != NULL);
        assert((selectedParents[c].Yc = malloc(in->C * sizeof(int))) != NULL);
        assert((selectedParents[c].Fc = malloc(in->C * sizeof(int))) != NULL);
        assert((selectedParents[c].eventualMemc = malloc(in->C * sizeof(int))) != NULL);
        assert((selectedParents[c].actualMemc = malloc(in->C * sizeof(int))) != NULL);
        assert((selectedParents[c].eventualGainc = malloc(in->C * sizeof(int))) != NULL);
        assert((selectedParents[c].actualGainc = malloc(in->C * sizeof(int))) != NULL);
        assert((selectedParents[c].confmem = malloc(in->C * sizeof(int))) != NULL);
    }

    assert( (temp->vettX = malloc(in->Q*sizeof(int))) != NULL);
    assert( (temp->vettC = malloc(in->Q*sizeof(float))) != NULL);
    assert( (temp->vettR = malloc(in->Q*sizeof(float))) != NULL);
    assert( (temp->Zi = malloc(in->I*sizeof(int))) != NULL);
    assert( (temp->Yc=malloc(in->C*sizeof(int))) != NULL);
    assert( (temp->Fc=malloc(in->C*sizeof(int))) != NULL );
    assert( (temp->eventualMemc=malloc(in->C*sizeof(int))) != NULL );
    assert( (temp->actualMemc=malloc(in->C*sizeof(int))) != NULL );
    assert( (temp->eventualGainc=malloc(in->C*sizeof(int))) != NULL );
    assert( (temp->actualGainc=malloc(in->C*sizeof(int))) != NULL );
    assert( (temp->confmem=malloc(in->C*sizeof(int))) != NULL );

    assert( (soluz->vettX = malloc(in->Q*sizeof(int))) != NULL);
    assert( (soluz->vettC = malloc(in->Q*sizeof(float))) != NULL);
    assert( (soluz->vettR = malloc(in->Q*sizeof(float))) != NULL);
    assert( (soluz->Zi = malloc(in->I*sizeof(int))) != NULL);
    assert( (soluz->Yc=malloc(in->C*sizeof(int))) != NULL);
    assert( (soluz->Fc=malloc(in->C*sizeof(int))) != NULL );
    assert( (soluz->eventualMemc=malloc(in->C*sizeof(int))) != NULL );
    assert( (soluz->actualMemc=malloc(in->C*sizeof(int))) != NULL );
    assert( (soluz->eventualGainc=malloc(in->C*sizeof(int))) != NULL );
    assert( (soluz->actualGainc=malloc(in->C*sizeof(int))) != NULL );
    assert( (soluz->confmem=malloc(in->C*sizeof(int))) != NULL );

    assert( (neighbor->vettX = malloc(in->Q*sizeof(int))) != NULL);
    assert( (neighbor->vettC = malloc(in->Q*sizeof(float))) != NULL);
    assert( (neighbor->vettR = malloc(in->Q*sizeof(float))) != NULL);
    assert( (neighbor->Zi = malloc(in->I*sizeof(int))) != NULL);
    assert( (neighbor->Yc=malloc(in->C*sizeof(int))) != NULL);
    assert( (neighbor->Fc=malloc(in->C*sizeof(int))) != NULL );
    assert( (neighbor->eventualMemc=malloc(in->C*sizeof(int))) != NULL );
    assert( (neighbor->actualMemc=malloc(in->C*sizeof(int))) != NULL );
    assert( (neighbor->eventualGainc=malloc(in->C*sizeof(int))) != NULL );
    assert( (neighbor->actualGainc=malloc(in->C*sizeof(int))) != NULL );
    assert( (neighbor->confmem=malloc(in->C*sizeof(int))) != NULL );

    assert( (*neibest = malloc(numN*sizeof(int))) != NULL);
    assert( (*neiworst = malloc(numN*sizeof(int))) != NULL);
    assert( (*chosen = malloc(in->C*sizeof(int))) != NULL);

}

void GAinit(Vett *pop, Instances *in)
{
    int j=0, c=0, q=0;
    int count, Memory, numconf;
    double iteration;

    for(j=0; j<numP; j++) {
        pop[j].feasible = -1;
        pop[j].gain = -1;
        iteration = 0;

        while ((pop[j].feasible == 0 || pop[j].gain <=  -exp(iteration)) && iteration < 6) //bisogna compilare con l'opzione -lm
        {

            iteration+=0.5;

            count = 0;
            Memory = in->M;
            numconf = 0;

            for (q = 0; q < in->Q; q++) {
                pop[j].vettX[q] = -1;
            }

            createZi(&pop[j], in);

            for (c = rand() % in->C; count < (in->C / 5); count++, c = rand() % in->C)
            {
                changeConfMem(&pop[j], in);
                if ((Memory - pop[j].confmem[c]) > 0) {
                    numconf++;
                    for (q = 0; q < in->Q; q++) {
                        if ((pop[j].vettX[q]) == -1 && (in->Gcq[c][q]) > 0) {
                            pop[j].vettX[q] = c;
                        } else if ((pop[j].vettX[q]) > -1) {
                            if (in->Gcq[c][q] > in->Gcq[pop[j].vettX[q]][q]) {
                                pop[j].vettX[q] = c;
                            }
                        }
                    }

                    Memory -= pop[j].confmem[c];
                    createZi(&pop[j], in);

                    calculateOF(&pop[j], in);
                    calculateC(&pop[j], in);
                    calculateR(&pop[j], in);
                    check3(&pop[j], in);
                    createYc(&pop[j], in);
                    changeFc(&pop[j], in);
                }
            }

        }

        pop[j].availableMem = in->M - pop[j].mem;

        for (c = 0, numconf = 0; c < in->C; c++) {
            if (pop[j].confmem[c] <= pop[j].availableMem) {
                pop[j].Yc[c] = 1;
                pop[j].availableMem -= pop[j].confmem[c];
                pop[j].confmem[c] = 0;
                numconf++;
            }
        }

        createVettXfromYc(&pop[j], in);
        createYc(&pop[j], in);
        createZifromYc(&pop[j], in);
        configureC(&pop[j], in);
        calculateOF(&pop[j], in);
        check3(&pop[j], in);
#if debug >= 2
        printf("HYB: configuzioni usate %d\n", numconf);
#endif
        pop[j].availableMem = in->M - pop[j].mem;
#if debug >= 2
        printf("HYB: confmem %d\n", pop[j].availableMem);
#endif
    }
}

void searchMax(Vett *pop, int dim, Sol *best, Instances *in, time_t time)
{
    int i;
    for(i=0; i<dim; i++)
    {
        if(best->gain < pop[i].gain && pop[i].feasible==1)
        {
            Vett2Sol(&pop[i], best, in);
            bestWrite(best,in);
#if debug >= 1
            printf("\nHYB: gain: %d\nHYB: memory: %d\nHYB: feasible: %d\nHYB: time: %ld\n", best->gain, best->mem, best->feasible, time);
#endif
        }
    }
}

void calculateFitness(Vett *v, Instances *in)
{
    v->fitness = v->gain - check3(v, in);
}

void mutation(Vett *v, Instances *in)
{
    int i, q, c = 0;
    float compatibility, maxCompatibility;
    int divisore, cMax;

    int qMin;
    float rMin;

    int flag = 0;
    float soglia = 0.17;

    int tabu[in->Q];
    for (q = 0; q < in->Q; ++q) {
        tabu[q] = 0;
    }
    q = 0;
    while (1)
    {
        rMin = INT_MAX;
        qMin = -1;
        c = 0;

        for (q = 0; q < in->Q; q++) {
            if (((v->vettC[q] != 0) ^ flag) && v->vettC[q] < rMin  && tabu[q] != 1)
            {
                qMin = q;
                rMin = v->vettC[q];
            }
        }

        if (qMin == -1)
        {
            if (flag == 1)
                break;
            flag = 1;
        } else {
            tabu[qMin] = 1;
            cMax = -1;
            maxCompatibility = 0;
            while (c < in->C)
            {
                if (in->Gcq[c][qMin] > 0)
                {
                    divisore = 0;
                    compatibility = 0;
                    for (i = 0; i < in->I; i++) {
                        if (v->Zi[i] == 1) {
                            divisore++;
                            if (in->Eci[c][i] == v->Zi[i]) {
                                compatibility++;
                            }
                        }
                    }
                    compatibility /= divisore;
                    if (compatibility > maxCompatibility) {
                        maxCompatibility = compatibility;
                        cMax = c;
                    }
                }
                c++;
            }
            if (flag == 1)
            {
                soglia = 0.21;
            }
            if (maxCompatibility >= soglia) {

                v->vettX[qMin] = cMax;
            }
            if (maxCompatibility < soglia) {
                v->vettX[qMin] = -1;

            }
            createZi(v, in);
            calculateOF(v, in);
            calculateC(v, in);
            calculateR(v, in);
            check3(v, in);
        }
    }

    int qIn;
    for (q = 0; q < in->Q; ++q) {

        if (v->vettX[q] != -1) {
            c = v->vettX[q];

            for (qIn = 0; qIn < in->Q; ++qIn) {
                if (v->vettX[(q + qIn) % in->Q] != -1)
                {
                    if (in->Gcq[c][(q + qIn) % in->Q] > in->Gcq[v->vettX[(q + qIn) % in->Q]][(q + qIn)% in->Q]) {
                        v->vettX[(q + qIn) % in->Q] = c;
                    }
                }
                else {
                    if (in->Gcq[c][(q + qIn) % in->Q] >0)
                        v->vettX[(q + qIn) % in->Q] = c;
                }
            }
        }
    }
    createZi(v, in);
    calculateOF(v, in);
    calculateC(v, in);
    calculateR(v, in);
    check3(v, in);

}

void selectParents(Vett *pop, Vett *selPar, Instances *in)
{
    int i;
    double sumFitness=0, prob, val;

    for (i=0; i<numP/2; i++) {
        selPar[i].feasible = -2;
    }

    for (i=0; i<numP; i++) {
        sumFitness += pop[i].fitness;
    }

    for (i=0; i<numP/2; i++) {
        prob = (pop[i].fitness/sumFitness)+probOffset;
        val =(double) rand()/RAND_MAX;
        if (val < prob) {
            cpySol(&selPar[i], &pop[i], in);
        }
    }
}

void crossover(Vett *p1, Vett *p2, Vett *c1, Vett *c2, Instances *in)
{
    int q;
    int x = (rand() % (in->Q - 1)) + 1;
    int i = x-1, j = x;
    float sumR1 = p1->vettR[i], sumR2 = p1->vettR[j];

    for(q = 0; q < in->Q; q++) {
        if (sumR1 != 0 && sumR2 != 0 &&
            (sumR1 > 2 * sumR2 || sumR2 > 2 * sumR1)) {
            break;
        }
        if (j < in->Q-1 && p1->vettR[j+1] > p2->vettR[j+1]) {
            j++;
            sumR1 += p1->vettR[j];
            sumR2 += p2->vettR[j];
        } else if (i > 0){
            i--;
            sumR1 += p1->vettR[i];
            sumR2 += p2->vettR[i];
        }
    }

    for (q = 0; q < i; q++) {
        c1->vettX[q] = p1->vettX[q];
        c2->vettX[q] = p2->vettX[q];
    }

    for (q = i; q < j; q++) {
        c1->vettX[q] = p2->vettX[q];
        c2->vettX[q] = p1->vettX[q];
    }

    for (q = j; q < in->Q; q++) {
        c1->vettX[q] = p1->vettX[q];
        c2->vettX[q] = p2->vettX[q];
    }
}

void sobstitute(Vett *pop, Vett *children, Instances *in)
{
    int i;
    double sumFitness=0, prob, val;

    for (i=0; i<numP; i++) {
        sumFitness += 1/children[i].fitness;
    }

    for (i=0; i<numP/2; i++) {
        prob = (children[i].fitness/sumFitness)+probOffset;
        val =(double) rand()/RAND_MAX;
        if (val < prob && children[i].feasible != -2) {
            cpySol(&pop[i], &children[i], in);
        }
    }
}

void calculateC(Vett *v, Instances *in)
{
    int q, i;
    float compatibility;
    int divisore;

    for (q = 0; q < in->Q; q++) {
        if (v->vettX[q] != -1) {
            compatibility = 0;
            divisore = 0;
            for (i = 0; i < in->I; i++) {
                if (v->Zi[i] == 1) {
                    divisore++;
                    if (in->Eci[v->vettX[q]][i] == v->Zi[i]) {
                        compatibility ++;
                    }
                }
            }
            v->vettC[q] = (float) compatibility / divisore;
        } else {
            v->vettC[q] = 0;
        }
    }
}

void calculateR(Vett *v, Instances *in)
{
    int q, i;
    int sum;

    for(q=0; q<in->Q; q++) {
        if(v->vettX[q] != -1) {
            for (i = 0, sum = 0; i < in->I; i++) {
                if (in->Eci[v->vettX[q]][i] == 1) {
                    sum += in->Mi[i] + in->Fi[i];
                }
            }
            v->vettR[q] = (float) in->Gcq[v->vettX[q]][q] / sum;
        } else {
            v->vettR[q] = 0;
        }
    }
}

void cpySol(Vett *dst, Vett *src, Instances *in)
{
    int i;

    dst->fitness = src->fitness;
    dst->feasible = src->feasible;
    dst->gain = src->gain;
    dst->mem = src->mem;
    dst->availableMem = src->availableMem;

    for (i=0; i<in->Q; i++) {
        dst->vettX[i] = src->vettX[i];
        dst->vettC[i] = src->vettC[i];
        dst->vettR[i] = src->vettR[i];
    }

    for (i=0; i<in->I; i++) {
        dst->Zi[i] = src->Zi[i];
    }

    for (i=0; i<in->C; i++) {
        dst->Yc[i] = src->Yc[i];
        dst->Fc[i] = src->Fc[i];
        dst->eventualMemc[i] = src->eventualMemc[i];
        dst->actualMemc[i] = src->actualMemc[i];
        dst->eventualGainc[i] = src->eventualGainc[i];
        dst->actualGainc[i] = src->actualGainc[i];
        dst->confmem[i] = src->confmem[i];
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
            if(in->Eci[j][i]==1 && x->Zi[i]==0)
            {
                x->confmem[j]+=in->Mi[i];
            }
        }
    }
}

void TS(Vett *localBest, Vett *temp, Vett *soluz, Vett *neighbor, int *chosen, int *neibest, int *neiworst,
        Instances *in)
{
    int c, q, i, iteration, j;
    int count;
    int bestc, worst;
    int bestVal, worstVal;
    int nof=0, fea=0;
    int flag=0;
    int tabu[2][8];

    cpySol(temp, localBest, in);
    cpySol(soluz, localBest, in);

#if debug >= 2
    if (temp->feasible != -2) {
        printf("\nHYB: Starting child\nHYB: ");
        for (q = 0; q < in->Q; q++) {
            printf("%d ", temp->vettX[q]);
        }
        printf("\nHYB: vettore Yc\nHYB: ");
        for (i = 0; i < in->C; i++) {
            printf("%d ", temp->Yc[i]);
        }
        printf("\nHYB: gain: %d\nHYB: memory: %d\nHYB: feasible:%d \n", temp->gain, temp->mem, temp->feasible);
    }
#endif

    iteration=0;

    count=0;

    while (iteration <= TSiter)
    {
        count++;

        for(i=0; i<numN; i++)
        {
            neibest[i]=-1;
            neiworst[i]=-1;
        }

        if(temp->feasible==1)
        {

            for(c=0; c<in->C; c++)
            {
                chosen[c]=0;
            }
            for(i=0; i<numN; i++)
            {
                bestc=-1;
                bestVal=-1000000000;
                worst=-1;
                worstVal=1000000000;

                for(c=0; c<in->C; c++)
                {
                    if(temp->Yc[c]==0 && chosen[c]==0 && bestVal < temp->eventualGainc[c])
                    {
                        bestVal = temp->eventualGainc[c];
                        bestc=c;
                    }
                    else if(temp->Yc[c]==1 && chosen[c]==0 && worstVal > temp->actualGainc[c])
                    {
                        worstVal=temp->actualGainc[c];
                        worst=c;
                    }
                }

                neibest[i]=bestc;
                neiworst[i]=worst;
                chosen[bestc]=1;
                chosen[worst]=1;
            }

            for(c=0;c<in->C; c++)
            {
                neighbor->Yc[c]=temp->Yc[c];
            }

            bestVal=-1000000;

            for(i=0; i<numN; i++)
            {
                if(neibest[i]!=-1)
                {
                    for(j=0; j<numN; j++)
                    {
                        flag=1;
                        if(neiworst[j]!=-1)
                        {
                            for(c=1;c<8;c++)
                            {
                                if( (tabu[0][c]==neibest[i] && tabu[1][c]==neiworst[j]) || (tabu[0][c]==neiworst[j] && tabu[1][c]==neibest[i]) )
                                {
                                    flag=0;
                                }
                            }
                            if(flag==1)
                            {

                                neighbor->Yc[neibest[i]]=1;
                                neighbor->Yc[neiworst[j]]=0;
                                createVettXfromYc(neighbor, in);
                                createZi(neighbor, in);
                                calculateOF(neighbor, in);
                                if(neighbor->gain > bestVal)
                                {
                                    bestVal=neighbor->gain;
                                    bestc=neibest[i];
                                    worst=neiworst[j];
                                }
                                neighbor->Yc[neibest[i]]=0;
                                neighbor->Yc[neiworst[j]]=1;
                            }
                        }
                    }
                }
            }

            temp->Yc[bestc]=1;
            temp->Yc[worst]=0;

            tabu[0][count]=bestc;
            tabu[1][count]=worst;
            nof=0;
            fea++;

        }
        else
        {
            while(temp->feasible==0)
            {
                for(c=0; c<in->C; c++)
                {
                    chosen[c]=0;
                }

                for(i=0; i<numN; i++)
                {
                    bestVal=10000000;
                    bestc=-1;
                    worstVal=-1;
                    worst=-1;

                    for(c=0; c<in->C; c++)
                    {
                        if(temp->Yc[c]==0 && chosen[c]==0 && bestVal > temp->eventualMemc[c] && temp->eventualMemc!=0)
                        {
                            bestVal = temp->eventualMemc[c];
                            bestc=c;
                        }
                        else if(temp->Yc[c]==1 && chosen[c]==0 && worstVal < temp->actualMemc[c])
                        {
                            worstVal=temp->actualMemc[c];
                            worst=c;
                        }
                    }

                    neibest[i]=bestc;
                    neiworst[i]=worst;
                    chosen[bestc]=1;
                    chosen[worst]=1;
                }


                for(c=0;c<in->C; c++)
                {
                    neighbor->Yc[c]=temp->Yc[c];
                }

                bestVal=-10000000;

                for(i=0; i<numN; i++)
                {
                    if(neibest[i]!=-1)
                    {
                        for(j=0; j<numN; j++)
                        {
                            flag=1;
                            if(neiworst[j]!=-1)
                            {
                                for(c=1;c<=7;c++)
                                {
                                    if( (tabu[0][c]==neibest[i] && tabu[1][c]==neiworst[j]) || (tabu[0][c]==neiworst[j] && tabu[1][c]==neibest[i]) )
                                    {
                                        flag=0;
                                    }
                                }
                                if(flag==1)
                                {
                                    neighbor->Yc[neibest[i]]=1;
                                    neighbor->Yc[neiworst[j]]=0;
                                    createVettXfromYc(neighbor, in);
                                    createZi(neighbor, in);
                                    calculateOF(neighbor, in);
                                    if(neighbor->gain > bestVal)
                                    {
                                        bestVal=neighbor->gain;
                                        bestc=neibest[i];
                                        worst=neiworst[j];
                                    }
                                    neighbor->Yc[neibest[i]]=0;
                                    neighbor->Yc[neiworst[j]]=1;
                                }
                            }
                        }
                    }
                }

                temp->Yc[bestc]=1;
                temp->Yc[worst]=0;

                tabu[0][count]=bestc;
                tabu[1][count]=worst;
                nof++;
                fea=0;

                createZifromYc(temp, in);
                check3(temp, in);

                if(nof>=20)
                {
                    nof=0;
                    while(temp->feasible==0)
                    {
                        worstVal=-1;
                        for(c=rand()%in->C;  temp->Fc[c]==1; c=rand()%in->C);
                        if(temp->Yc[c]==1)
                        {
                            temp->Yc[c]=0;
                        }
                        nof=0;
                        createZifromYc(temp, in);
                        check3(temp, in);
                    }

                }
            }
        }

        createSolution(temp,in);

        if(fea>=5)
        {

            for(q=rand()%in->Q, flag=0; temp->vettX[q]==-1 && flag<in->Q; q=rand()%in->Q, flag++);

            for(c=rand()%in->C, flag=0; flag<in->C; c=rand()%in->C, flag++)
            {
                if(in->Gcq[c][q]>0 && temp->Yc[c]==0)
                {
                    temp->Yc[c]=1;
                    flag=in->C;
                }
            }

            createSolution(temp,in);
            fea=0;
        }

        if(temp->gain > localBest->gain && temp->feasible==1)
        {
            cpySol(neighbor, temp, in);

            bestVal = 0;
            bestc = 0;

            while(bestc>-1)
            {
                bestc = -1;
                bestVal = 0;
                for (c = 0; c < in->C; c++) {
                    if (neighbor->eventualGainc[c] > bestVal && neighbor->eventualMemc[c] < (in->M - neighbor->mem)) {
                        bestVal = neighbor->eventualGainc[c];
                        bestc = c;
                    }
                    if (bestc > -1) {
                        neighbor->Yc[bestc] = 1;
                        createSolution(neighbor, in);
                    }
                }
            }
            if (neighbor->gain > soluz->gain) {
                cpySol(soluz, neighbor, in);
            }
            if (soluz->gain > localBest->gain) {
                cpySol(localBest, soluz, in);
#if debug >= 2
                printf("HYB: TS\nHYB: gain: %d\nHYB: memory: %d\nHYB: feasible: %d\n", soluz->gain, soluz->mem, soluz->feasible);
#endif
            }
        }

        iteration++;

        if(count==7)
        {
            count=0;
        }

    }
}

char *GetFileName(const char *path)
{
    char *filename = strrchr(path, '\\');
    if (filename == NULL)
        filename = (char *) path;
    else
        filename++;
    return filename;
}