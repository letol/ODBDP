#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <limits.h>

#define numP 250
#define debug 1
#define probOffset .1


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
    float *vettC;
    float *vettR;
    float fitness;
    int mem;
    int gain;
    int *Zi;
    int feasible;
}Vett;

void initialization(FILE *fin, Instances *in, Sol *best, Vett *pop, Vett *selectedParents, Vett *children);
void letturavet(int *v, FILE *fin, int r);
void letturamat(int **m, FILE *fin, int r, int c);
void calculateOF(Vett *v, Instances *in);
int relax3(Vett *v, Instances *in);
void GAinit(Vett *pop, Instances *in);
void createZi(Vett *v, Instances *in);
void searchMax(Vett *pop, Sol *best, Instances *in);
void selectParents(Vett *pop, Vett *selPar, Instances *in);
void crossover(Vett *p1, Vett *p2, Vett *c1, Vett *c2, Instances *in);
void sobstitute(Vett *pop, Vett *children, Instances *in);
void calculateC(Vett *v, Instances *in);
void calculateR(Vett *v, Instances *in);
void calculateFitness(Vett *v, Instances *in);
void Vett2Sol(Vett *v, Sol *s, Instances *in);
void invertion(Vett *v, Instances *in);
void mutation(Vett *v, Instances *in);
void cpySol(Vett *dst, Vett *src, Instances *in);
void changeConfMem(int *confmem, Vett *x, Instances *in);
//void averageGain(int *av,Instances *in);

int main(int argc, char* argv[])
{
    int j, c, q, i, iteration=0, Xcount=0;
    Instances in;
    Sol best;
    FILE *fin, *fout;
    time_t start=time(NULL);
    int timelimit=0;
    Vett population[numP];
    Vett selectedParents[numP/2], children[numP/2];

#if debug == 0
    srand((unsigned int)time(NULL)); //per aumentare casualità quando non siamo in debug
#endif

    assert(argc == 4);

    assert(strcmp(argv[2], "-t") == 0);

    timelimit = atoi(argv[3]);

    fin = fopen(argv[1], "r");

    assert(fin != NULL);

    initialization(fin, &in, &best, population, selectedParents, children);

    GAinit(population, &in);

#if debug >= 2
    printf("INITIAL POPULATION\n\n");
    for(j=0; j<numP; j++)
    {
        printf("\n");
        for(q=0; q<in.Q; q++)
        {
            printf("%d ", population[j].vettX[q]);
        }
        printf("\nvettore zi%d\n ", j);
        for(i=0; i<in.I; i++)
        {
            printf("%d ", population[j].Zi[i]);
        }
        printf("\nvettore c%d\n ", j);
        for(i=0; i<in.I; i++)
        {
            printf("%f ", population[j].vettC[i]);
        }
        printf("\nvettore r%d\n ", j);
        for(i=0; i<in.I; i++)
        {
            printf("%f ", population[j].vettR[i]);
        }
        printf("\ngain: %d\nmemory: %d\nfeasible:%d \n", population[j].gain, population[j].mem,population[j].feasible);
    }
#endif

    best.gain=0;

    searchMax(population, &best, &in);

    for(i=0; i<numP; i++) {
        calculateFitness(&population[i], &in);
        mutation(&population[i], &in);
    }

    while ((time(NULL) - start)<=timelimit)
    {
        selectParents(population, selectedParents, &in);
#if debug >= 1
        printf("\nSELECTED PARENTS\n");
        for(j=0; j<numP/2; j++) {
            if (selectedParents[j].feasible != -2) {
                printf("\nP%d\n", j);
                for (q = 0; q < in.Q; q++) {
                    printf("%d ", selectedParents[j].vettX[q]);
                }
                printf("\nvettore zi%d\n ", j);
                for (i = 0; i < in.I; i++) {
                    printf("%d ", selectedParents[j].Zi[i]);
                }
                printf("\nvettore c%d\n ", j);
                for (i = 0; i < in.Q; i++) {
                    printf("%f ", selectedParents[j].vettC[i]);
                }
                printf("\nvettore r%d\n ", j);
                for (i = 0; i < in.Q; i++) {
                    printf("%f ", selectedParents[j].vettR[i]);
                }
                printf("\ngain: %d\nmemory: %d\nfeasible:%d \n", selectedParents[j].gain, selectedParents[j].mem, selectedParents[j].feasible);
            }
        }
#endif

#if debug >= 1
        printf("\nCROSSOVER & MUTATION\n");
#endif
        for (i=0; i<numP/2; i++) {
            children[i].feasible = -2;
        }
        i=0;
        while (i<numP/2) {
            if (selectedParents[i].feasible != -2) {
                j=i+1;
                while (j<numP/2) {
                    if (selectedParents[j].feasible != -2) {
                        crossover(&selectedParents[i], &selectedParents[j], &children[i], &children[j], &in);
                        Xcount++;

                        calculateC(&children[i], &in);
                        calculateR(&children[i], &in);
                        calculateC(&children[j], &in);
                        calculateR(&children[j], &in);

                        mutation(&children[i], &in);
                        mutation(&children[j], &in);

                        calculateC(&children[i], &in);
                        calculateR(&children[i], &in);
                        calculateC(&children[j], &in);
                        calculateR(&children[j], &in);

                        createZi(&children[i], &in);
                        calculateOF(&children[i], &in);
                        calculateFitness(&children[i], &in);

                        createZi(&children[j], &in);
                        calculateOF(&children[j], &in);
                        calculateFitness(&children[j], &in);
                    }
                    j++;
                }
            }
            i++;
        }

#if debug >= 1
        printf("\nPRODUCED CHILDREN\n");
        for(j=0; j<numP/2; j++) {
            if (children[j].feasible != -2) {
                printf("\nC%d\n", j);
                for (q = 0; q < in.Q; q++) {
                    printf("%d ", children[j].vettX[q]);
                }
                printf("\nvettore zi%d\n ", j);
                for (i = 0; i < in.I; i++) {
                    printf("%d ", children[j].Zi[i]);
                }
                printf("\nvettore c%d\n ", j);
                for (i = 0; i < in.Q; i++) {
                    printf("%f ", children[j].vettC[i]);
                }
                printf("\nvettore r%d\n ", j);
                for (i = 0; i < in.Q; i++) {
                    printf("%f ", children[j].vettR[i]);
                }
                printf("\ngain: %d\nmemory: %d\nfeasible:%d \n", children[j].gain, children[j].mem, children[j].feasible);
            }
        }
#endif

        sobstitute(population, children, &in);


        //If no crossover happened, apply mutation over population
        if(Xcount == 0) {
#if debug >= 1
            printf("\nNO CROSSOVER -> ONLY MUTATION\n");
#endif
            for(i=0; i<numP; i++) {
                mutation(&population[i], &in);
            }
#if debug >= 2
            for(j=0; j<numP; j++)
            {
                printf("\n");
                for(q=0; q<in.Q; q++)
                {
                    printf("%d ", population[j].vettX[q]);
                }
                printf("\nvettore zi%d\n ", j);
                for(i=0; i<in.I; i++)
                {
                    printf("%d ", population[j].Zi[i]);
                }
                printf("\nvettore c%d\n ", j);
                for(i=0; i<in.Q; i++)
                {
                    printf("%f ", population[j].vettC[i]);
                }
                printf("\nvettore r%d\n ", j);
                for(i=0; i<in.Q; i++)
                {
                    printf("%f ", population[j].vettR[i]);
                }
                printf("\ngain: %d\nmemory: %d\nfeasible:%d \n", population[j].gain, population[j].mem,population[j].feasible);
            }
#endif
        }
        Xcount = 0;

        searchMax(population, &best, &in);

#if debug >= 1
        printf("\nCurrent Best solution stats:\n");

        printf("\ngain: %d\nmemory: %d\nfeasible: %d", best.gain, best.mem, best.feasible);
        printf("\n\nIteration: %d", ++iteration);
#endif
    }

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

void initialization(FILE *fin, Instances *in, Sol *best, Vett *pop, Vett *selectedParents, Vett *children)
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
        assert( (pop[c].vettC = malloc(in->Q*sizeof(float))) != NULL);
        assert( (pop[c].vettR = malloc(in->Q*sizeof(float))) != NULL);
        assert( (pop[c].Zi = malloc(in->I*sizeof(int))) != NULL);
    }

    for(c=0; c<numP/2; c++)
    {
        assert((children[c].vettX = malloc(in->Q * sizeof(int))) != NULL);
        assert((children[c].vettC = malloc(in->Q * sizeof(float))) != NULL);
        assert((children[c].vettR = malloc(in->Q * sizeof(float))) != NULL);
        assert((children[c].Zi = malloc(in->I * sizeof(int))) != NULL);

        assert((selectedParents[c].vettX = malloc(in->Q * sizeof(int))) != NULL);
        assert((selectedParents[c].vettC = malloc(in->Q * sizeof(float))) != NULL);
        assert((selectedParents[c].vettR = malloc(in->Q * sizeof(float))) != NULL);
        assert((selectedParents[c].Zi = malloc(in->I * sizeof(int))) != NULL);
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
// Calculates Objective Function value of solution v
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
// Genetic Algorithm initialization
{
    int j=0, i=0,c=0, q=0;
    int count, Memory, numconf;
    int *confmem;
    //int av=0;
    assert(confmem=malloc(in->C*sizeof(int)));
    //averageGain(&av, in);

    for(j=0; j<numP; j++)
    {
        pop[j].feasible=-1;
        pop[j].gain=-1;

        while(pop[j].feasible==0 || pop[j].gain <=0)
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
                        confmem[c]+=in->Mi[i];
                    }
                }
            }
            for(q=0; q<in->Q;q++)
            {
                pop[j].vettX[q]=-1;
            }

            createZi(&pop[j], in);

            for(c=rand()%in->C; count<(in->C/5); count++, c=rand()%in->C)
            {
                if( (Memory-confmem[c]) > 0)
                {
                    numconf++;
                    for(q=0;q<in->Q;q++)
                    {
                        if((pop[j].vettX[q]) == -1 && (in->Gcq[c][q]) > 0 )
                        {
                            pop[j].vettX[q]=c;
                        }
                        else if( (pop[j].vettX[q]) > -1)
                        {
                            if(in->Gcq[c][q] > in->Gcq[pop[j].vettX[q]][q])
                            {
                                pop[j].vettX[q]=c;
                            }
                        }
                    }

                    Memory-=confmem[c];
                    createZi(&pop[j], in);
                    changeConfMem(confmem, &pop[j], in);
                    calculateOF(&pop[j], in);
                    calculateC(&pop[j], in);
                    calculateR(&pop[j], in);
                    relax3(&pop[j], in);
                }
            }
        }
#if debug >= 1
        printf("configuzioni usate %d\n", numconf);
#endif
    }
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
    return;
}

void searchMax(Vett *pop, Sol *best, Instances *in)
// Scan population for best solution
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

void crossover(Vett *p1, Vett *p2, Vett *c1, Vett *c2, Instances *in)
//scorro i due vettori contemporaneamente appena una dei due prevarica l'altro in fatto di gain offerto scelgo quello come punto di taglio
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
        if (j < in->Q-1 && p1->vettR[j+1] > p2->vettR[j-1]) {
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

void calculateC(Vett *v, Instances *in)
// Calculates vettC from vettX - compatibility based: guardo quanti indici attivi ha in comune la configurazione con Zi corrente
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
// Calculates vettR from vettX
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

void calculateFitness(Vett *v, Instances *in)
{
    v->fitness = v->gain - relax3(v, in);
}

void Vett2Sol(Vett *v, Sol *s, Instances *in)
// Translate type Vett into type Sol solution
{
    int i, c, q;
    for(c=0; c<in->C; c++){
        for(q=0; q<in->Q; q++) {
            if(c == v->vettX[q]) {
                s->Xcq[c][q] = 1;
            } else {
                s->Xcq[c][q] = 0;
            }
        }
    }
    for(i=0;i-in->I;i++) {
        s->Zi = v->Zi;
    }
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
/*
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
 */

void selectParents(Vett *pop, Vett *selPar, Instances *in)
// Select parents with fitness based probability
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

void sobstitute(Vett *pop, Vett *children, Instances *in)
// Substitute parents with children fitness based probability
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

void invertion(Vett *v, Instances *in)
// Invertion between genes with minimum r
{
    int q;
    int q1=-1, q2=-1;
    float r1=INT_MAX, r2=INT_MAX;
    int xTmp;
    float rTmp;

    for(q=0; q<in->Q; q++) {
        if (v->vettC[q] != 0) {
            if (v->vettC[q] < r2) {
                if (v->vettC[q] < r1) {
                    r2 = r1;
                    q2 = q1;
                    r1 = v->vettC[q];
                    q1 = q;
                } else {
                    r2 = v->vettC[q];
                    q2 = q;
                }
            }
        }
    }

    xTmp = v->vettX[q1];
    rTmp = v->vettC[q1];
    v->vettX[q1] = v->vettX[q2];
    v->vettC[q1] = v->vettC[q2];
    v->vettX[q2] = xTmp;
    v->vettC[q2] = rTmp;
}

void mutation(Vett *v, Instances *in)
//criticità: una conf non serve tutte le query che potrebbe, bisognerebbe dare una passata alla fine per vedere se c'è qualche query non servita con una sua conf attiva
//il valore soglia è empirico
//cambiare la soglia quando flag = 1 seconda la memoria disponibile, in generale credo proprio che la soglia vada in base alla mem
{
    int i, q, c = 0;
    float compatibility, maxCompatibility;
    int divisore,cMax;

    int qMin;
    float rMin;

    int flag = 0;
    float soglia = 0.15; //valore empirico


    int tabu[in->Q];
    for (int var = 0; var < in->Q; ++var) {
        tabu[var] = 0;
    }
    while (1) //ne esco dopo aver considerato una e una sola volta tutte le query
    {
        rMin = INT_MAX;
        qMin = -1;
        c = 0;
        compatibility = 0;

        for (q = 0; q < in->Q; q++) {
            if (((v->vettC[q] != 0) ^ flag) //voglio occuparmi prima delle query con già una config e poi di quelle senza
                && v->vettC[q] < rMin  //ricerca minimo
                //&& v->vettC[q] < soglia //se voglio escludere di cambiare le q con una compatibilità > soglia
                && tabu[q] != 1) //voglio fare un solo giro di vettX
            {
                qMin = q;
                rMin = v->vettC[q];
            }
        }

        if (qMin == -1) // mi permette di occuparmi prima delle query già servite da una config e poi delle altre
        {
            if (flag == 1)
                break;
            flag = 1;
        } else {
            tabu[qMin] = 1;
            cMax = -1;
            maxCompatibility = 0;
            while (c < in->C) //mi guardo tutte le configurazioni che con questa query mi danno gain positivo alla ricerca di quella con compatibilità max
            {
                if (in->Gcq[c][qMin] > 0) //calcolo conf con compatibilità max
                {
                    divisore = 0;
                    compatibility = 0;
                    for (i = 0; i < in->I; i++) {
                        if (v->Zi[i] == 1) {
                            divisore++;
                            if (in->Eci[c][i] == v->Zi[i]) {
                                compatibility ++;
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
            //printf("max %f\n", maxCompatibility);
            if (flag == 1) //soglia per le query che non erano ancora servite, più alta per non aggiurne troppe e andare fuori con la memoria (sarebbe bello decidere in base a quanta memoria mi rimane a disp)
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
            relax3(v, in);
        }
    }
}

void cpySol(Vett *dst, Vett *src, Instances *in)
// Copy solution from src to dst
{
    int i;

    dst->fitness = src->fitness;
    dst->feasible = src->feasible;
    dst->gain = src->gain;
    dst->mem = src->mem;

    for (i=0; i<in->Q; i++) {
        dst->vettX[i] = src->vettX[i];
        dst->vettC[i] = src->vettC[i];
        dst->vettR[i] = src->vettR[i];
    }

    for (i=0; i<in->I; i++) {
        dst->Zi[i] = src->Zi[i];
    }
}