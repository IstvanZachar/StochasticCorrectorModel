// Stochastic corrector model with fusion
// István Zachar
// 2017-2018

// Version history
//  5.0 - Switch to turn on fitness-dependent fusion
//  6.0 - Fitness dependent fusion was changed to be inversely proportionate to fitness (the smaller fitness is, the more fusion is likely)
//  7.0 - Skipped.
//  8.0 - Added a new statistic: MeanFitnessNondead
//  9.0 - Added fitness dependent fusion's limit as external parameter;
//      - Fixed: added missing test for checking if there are at least two cells matching the condition of fitness-dependent fusion
//      - Fixed: if simulation terminates before maxt is reached, duplicate last point and terminate.
// 10.0 - New fitness function: flux reduction involves junk in the numerator besides parasite.
//      - Fixed: Corrected a marginal bug in fitness calculation converting the `prod` term to integers in every iteration step
//      - Fixed: Corrected calculation of `distT`.
//      - Fixed: Saving 0.0 for time-averaged statistics if distT==0 (usually at first step).


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h> // `access`

#include "randomGenerator.c"


#define VERSION	(10)
#define CONFIG	("config.cfg")
#define OUTPUT	("data.dat")
#define MAXE    (20) // maximal number of enzyme types
#define MAXC    (10000) // maximal number of cells
#define SKIP    (fscanf(file, "%*[^\n]\n")) // skip rest of line to next line


// Globals

int N, PS, sumD = 0, sumV = 0;
double WMAX, sumW = 0, minW = 1, maxW = 0, act[2];

typedef struct CELL {
  double w;
	int    v, d, good[MAXE], med[MAXE], junk, par;
} CELL;
CELL cells[MAXC];



// Functions

int randomSplit(int n) { // random variable from a binomial distribution (between 0 and `n`)
	int i, k = 0;
	for(i = 0; i < n; i++) if(randd() < 0.5) k++;
	return(k);
}

int clip(int x, int min, int max) { // clip `x` between `(min, max)`
	if (x > max) return(max); else if (x < min) return(min); else return(x);
}

int selectCell(double max) { // select a cell with probability proportional to its fitness
	int i = 0;
	double sum = 0, r = randd()*max;
	while((sum < r) && (i < PS)) {
		sum += cells[i].w;
		i++;	
	}
	return(i-1);
}

int selectCellInv(double max) { // select a cell with probability inversely proportional to its fitness
	int i = 0;
	double invmax, sum = 0, r;
	invmax = (PS*WMAX)-max;
	r = randd()*invmax;
	while((sum < r) && (i < PS)) {
		sum += (WMAX - cells[i].w);
		i++;	
	}
	return(i-1);
}

int newRandomChoice(double *v, double max) { // randomly select an index of vector `v` with probability proportional to v[i]
	int i = 0;
	double r = randd()*max, sum = v[0];
	while(sum < r) {
		i++;
		sum += v[i];
	}
	return(i);
}

void printCell(int c) {
	int i;
	printf("C%d:\tG[", c);
	for(i = 0; i < N; i++) printf("%d ", cells[c].good[i]);
	printf("]\tM[");
	for(i = 0; i < N; i++) printf("%d ", cells[c].med[i]);
	printf("]\tJ[%d]\tP[%d]\t", cells[c].junk, cells[c].par);  
	printf("D%d\tV%d\tW%lf\n", cells[c].d, cells[c].v, cells[c].w);
	fflush(stdout);
}

void updateCell(int c) {
	int i, j, e = 0, n = 0, deadQ = 0, oldD = cells[c].d, oldV = cells[c].v;
	double a = 0, oldW = cells[c].w, p = 1;
	//printf("BEFORE: %lf\t", sumW); printCell(c); // Note: after molecular replication, cell stats are supposed to be NOT RIGHT, as they need to be recalculated!
	for(i = 0; i < N; i++) {
		e = cells[c].good[i] + cells[c].med[i];
		if (e > 0) {
			a = (act[0]*cells[c].good[i]) + (act[1]*cells[c].med[i]);
			if (a > 0) {
				p *= a; // calculating product ((A_gi*E_gi) + (A_mi*E_mi)) iteratively
				n += e; // calculating total enzyme count
			} else { // actual combined enzymatic activity is zero
				deadQ = 1; break; // no need to check enzymes further
			}
		} else { // missing both states of essential enzyme
			deadQ = 1; break; // no need to check enzymes further
		}
	}
	if (deadQ) { // cell is dead
		cells[c].d = 1;
		cells[c].v = 0;
		cells[c].w = 0;
  } else {     // cell is not dead
		cells[c].d = 0;
		cells[c].v = n + cells[c].junk + cells[c].par; // molecule count ("volume")
		// cells[c].w = (p/WMAX)*(1 - ((cells[c].junk + cells[c].par)/(double)cells[c].v)); // old fitness
		cells[c].w = p/WMAX; // new fitness
	}
	sumD = sumD + (cells[c].d - oldD);
	sumW = sumW + (cells[c].w - oldW);
	sumV = sumV + (cells[c].v - oldV);
	//printf("AFTER:  %lf\t", sumW); printCell(c);
}

void splitCell(int p, int o) {
	int i;
	for(i = 0; i < N; i++) {
		cells[o].good[i] = randomSplit(cells[p].good[i]);
		cells[p].good[i] = cells[p].good[i] - cells[o].good[i];
		cells[o].med[i]  = randomSplit(cells[p].med[i]);
		cells[p].med[i]  = cells[p].med[i] - cells[o].med[i];
	}
	cells[o].junk = randomSplit(cells[p].junk);
	cells[p].junk = cells[p].junk - cells[o].junk;
	cells[o].par  = randomSplit(cells[p].par);
	cells[p].par  = cells[p].par - cells[o].par;
	updateCell(p);
	updateCell(o);
}

void minMaxFitness(int t) { // finds minimum and maximum values of fitness, modifying globals 'minW' and 'maxW'
	int i;
	minW = 1;
	maxW = 0;
	for(i = 0; i < PS; i++) {
		if(cells[i].w < minW) {
			minW = cells[i].w; 
		} else {
			if(cells[i].w > maxW) maxW = cells[i].w;
		}
	}
	// printf("T%d   m = %lf M = %lf    [", t, minW, maxW);
	// for(i = 0; i < N; i++) printf("%lf%s", 100000*cells[i].w, (i==(N-1)?"]":" "));
	// printf("\n");
	
}


int main(int argc, char** argv) { // Main
	
	unsigned long int randomseed = 1;
	int     tmax, smax, sd, init[4], fs, saveAt, i, j, t = 0, c = -1, fdf = 0, oldT = 0;
	double  fusion = 0, limit = 1, rep[4], mutg[4], mutm[4], mutj[4], mutp[4], div;
	double  *propensities = NULL;
	char    fileName[1024];
	FILE    *file = NULL;
	int     eventC = 0, event[(2*MAXE)+2] = {0}, mutC[4][4] = {{0}};
	int     splitCount = 0, fusionCount = 0;
	
	{ // Read parameters
		if (access(CONFIG, F_OK) != -1) {
			file = fopen(CONFIG, "r");
			fscanf(file, "%lu", &randomseed); SKIP; // `signed long`: `seed` only accepts positive `long`, does not accept 0
			fscanf(file, "%d",  &tmax);       SKIP; // maximal number of time steps
			fscanf(file, "%d",  &smax);       SKIP; // maximal number of stored steps
			fscanf(file, "%d",  &PS);         SKIP; // population size
			fscanf(file, "%d",  &sd);         SKIP; // split density
			fscanf(file, "%d",  &N);          SKIP; // number of enzyme types
			fscanf(file, "%lf", &fusion);     SKIP; // fusion probability
			fscanf(file, "%d",  &fdf); 				SKIP; // fitness dependent fusion
			fscanf(file, "%lf", &limit); 			SKIP; // fitness dependent fusion's fitness limit
			// Read required amount of entries and then skip to next line.
			for(i = 0; i < 2; i++) fscanf(file, "%lf", &act[i]);  SKIP; // enzyme activity
			for(i = 0; i < 4; i++) fscanf(file, "%lf", &rep[i]);  SKIP; // enzyme replication affinity
			for(i = 0; i < 4; i++) fscanf(file, "%d",  &init[i]); SKIP; // max values for initial random range of (0, init[i])
			for(i = 0; i < 4; i++) fscanf(file, "%lf", &mutg[i]); SKIP; // mutation rates from good enzyme
			for(i = 0; i < 4; i++) fscanf(file, "%lf", &mutm[i]); SKIP; // mutation rates from medium enzyme
			for(i = 0; i < 4; i++) fscanf(file, "%lf", &mutj[i]); SKIP; // mutation rates from junk
			for(i = 0; i < 4; i++) fscanf(file, "%lf", &mutp[i]); SKIP; // mutation rates from parasite
			fclose(file);
		} else { printf("ERROR. Parameter file %s does not exist. Aborting.\n", CONFIG); exit(0); }

	}
	
	
	{ // Assertions
		double gSum = 0, mSum = 0, jSum = 0, pSum = 0;
		if (PS < 2) { printf("ERROR. Population size %d is too small. Aborting.\n", PS); exit(0); }
		if (randomseed < 1) { printf("ERROR. Random seed %lu must be a positive integer. Aborting.\n", randomseed); exit(0); }
		for(i = 0; i < 4; i++) gSum += mutg[i];
		if (gSum != 1) { printf("ERROR. Mutation rates of good enzyme do not sum up to 1. Aborting.\n"); exit(0); }
		for(i = 0; i < 4; i++) mSum += mutm[i];
		if (mSum != 1) { printf("ERROR. Mutation rates of medium enzyme do not sum up to 1. Aborting.\n"); exit(0); }
		for(i = 0; i < 4; i++) jSum += mutj[i];
		if (jSum != 1) { printf("ERROR. Mutation rates of junk enzyme do not sum up to 1. Aborting.\n"); exit(0); }
		for(i = 0; i < 4; i++) pSum += mutp[i];
		if (pSum != 1) { printf("ERROR. Mutation rates of parasite do not sum up to 1. Aborting.\n"); exit(0); }
		if ((init[0]+init[1]) == 0) { printf("ERROR. The sum of the initial ranges for good and medium enzymes must be at least 1, otherwise there will be no enzyme in the system.\n"); exit(0); }
	}
	
	
	{ // Initialization
		seed(randomseed);
		fs = (log10((double)tmax)+1); // field size of time
		saveAt = clip(tmax/smax, 1, 5000); // save at every `saveAt` timestep
		//WMAX = exp(N*log((double)sd/(double)N))/(double)sd; // old fitness calculation
		WMAX = pow((sd*act[0])/(double)N, N);
		div = (double)(PS*sd); // divisor, used for generating population averages
		propensities = (double*)calloc((N*2)+2, sizeof(double));
		file = fopen(OUTPUT, "w");
	}
	
	
	if (0) { // Verbosing
		printf("%lu %d %d %d %d %f %d %d %lf\n\n", randomseed, tmax, saveAt, PS, sd, fusion, fdf, N, WMAX);
		printf("random seed: %lu\n", randomseed);
		printf("max t:       %d\n",  tmax);
		printf("save at:     %d\n",  saveAt);
		printf("pop. size:   %d\n",  PS);
	  printf("enzyme no.:  %d\n",  N);
		printf("split dens.: %d\n",  sd);
		printf("fusion rate: %lf\n", fusion);
		printf("fusion type: %d\n",  fdf);
		printf("fusion lim.: %lf\n", limit);
		printf("max fitness: %lf\n", WMAX);
		printf("activity: "); for(i = 0; i < 2; i++) printf("%lf\t", act[i]);  printf("\n");
		printf("rep.rate: "); for(i = 0; i < 4; i++) printf("%lf\t", rep[i]);  printf("\n");
		printf("initials: "); for(i = 0; i < 4; i++) printf("%d\t",  init[i]); printf("\n");
		printf("mut.good: "); for(i = 0; i < 4; i++) printf("%lf\t", mutg[i]); printf("\n");
		printf("mut.med : "); for(i = 0; i < 4; i++) printf("%lf\t", mutm[i]); printf("\n");
		printf("mut.junk: "); for(i = 0; i < 4; i++) printf("%lf\t", mutj[i]); printf("\n");
		printf("mut.par : "); for(i = 0; i < 4; i++) printf("%lf\t", mutp[i]); printf("\n");
		printf("\n");
		fflush(stdout);
	}


	{ // Initial population
		for(i = 0; i < PS; i++) {
			cells[i].v = 0;
			cells[i].d = 0;
			cells[i].w = 0;
			while ((cells[i].v < 1) || (cells[i].v >= sd)) { // don't generate empty or full cells
				for(j = 0; j < N; j++) {
					cells[i].good[j] = 0;
					cells[i].med[j]  = 0;
					 // Note that actual range of molecular amounts will be (0, init[i])
					 // as `randl` returns between (0, n-1)!
					if(init[0]) cells[i].good[j] = clip((int)randl(init[0]+1), 0, sd);
					if(init[1]) cells[i].med[j]  = clip((int)randl(init[1]+1), 0, sd);
				}
				cells[i].junk = 0;
				cells[i].par  = 0;
				if(init[2]) cells[i].junk = clip((int)randl(init[2]+1), 0, sd);
				if(init[3]) cells[i].par  = clip((int)randl(init[3]+1), 0, sd);
				updateCell(i);
			}
			//printCell(i);
		}
	}

		
	while (t <= tmax) {
				
		if (((t % saveAt) == 0) || (t == tmax) || (sumD == PS))  { // Save intermediate results
			int goodt[MAXE] = {0}, medt[MAXE] = {0}, junkt = 0, part = 0, distT = t - oldT;
			double minW = 1, maxW = 0;
			
			oldT = t;
			
			//printf("%d\n", t);
			//printf("%d \r", t); // Monitor progress
			
			for(i = 0; i < PS; i++) if (cells[i].d == 0) {
				for(j = 0; j < N; j++) {
					goodt[j] += cells[i].good[j];
					medt[j]  += cells[i].med[j];
				}
				junkt += cells[i].junk;
				part  += cells[i].par;
				if(cells[i].w < minW) {
					minW = cells[i].w; 
				} else {
						if(cells[i].w > maxW) maxW = cells[i].w;
				}
			}
			fprintf(file, "%d\t", t);
			for(j = 0; j < N; j++) fprintf(file, "%lf\t", goodt[j]/div);
			for(j = 0; j < N; j++) fprintf(file, "%lf\t", medt[j]/div);
			fprintf(file, "%lf\t%lf\t", junkt/div, part/div);
			fprintf(file, "%lf\t%lf\t%lf\t", minW, maxW, sumW/(double)PS);
			fprintf(file, "%lf\t%lf\t%lf\t", sumV/div, sumD/(double)PS, sumW/(double)(PS-sumD));
			if(distT) fprintf(file, "%lf\t%lf\t", (double)splitCount/(double)distT, (double)fusionCount/(double)distT);
			else      fprintf(file, "0.0\t0.0\t");
			fprintf(file, "\n");
			fflush(file);
			
			/*
			if(0) {
				printf("%d ", t);
				for(j = 0; j < N; j++) printf("%lf ", goodt[j]/(double)PS);
				for(j = 0; j < N; j++) printf("%lf ", medt[j]/(double)PS);
				printf("%lf %lf ", junkt/(double)PS, part/(double)PS);
				printf("%lf %lf %lf ", minW, maxW, sumW/(double)PS);
				printf("%lf %lf\n", sumV/(double)PS , sumD/(double)PS);
				printf("%lf %d %d\n", sumW, sumV, sumD);
				for(i = 0; i < PS; i++) { printf("\t"); printCell(i); };
				fflush(stdout);
			}
			*/
			
			splitCount  = 0;
			fusionCount = 0;
		}
		
		if(sumD == PS) {
			if(t == tmax) break;
			t = tmax; //break; // if all cells are dead, terminate
		} else {
		
			c = selectCell(sumW);
						
			{ // Select reaction and multiply & mutate molecule
				int r = 0, from = -9, to = -9;
				double sumP = 0;
				
				for(i = 0; i < ((2*N)+2); i++) propensities[i] = 0;
				
				/*
				if (0) {
					printf("T%d\tD%d\tV%d\tW%lf\t", t, sumD, sumV, sumW);
					printf("C%d\td%d\tv%d\tw%lf\t", c, cells[c].d, cells[c].v, cells[c].w);
					printf("M[");
					for(i = 0; i < N; i++) printf("%d ", cells[c].good[i]);
					for(i = 0; i < N; i++) printf("%d ", cells[c].med[i]);
					printf("%d %d]\t", cells[c].junk, cells[c].par);
					printf("P["); for(i = 0; i < ((2*N)+2); i++) printf("%lf ", propensities[i]); printf("]");
					printf("\tP%lf\t%d: %d -> %d\n", sumP, r, from, to);
					fflush(stdout);
				}
				*/
				
				propensities[2*N]     = (rep[2]*cells[c].junk);
				propensities[(2*N)+1] = (rep[3]*cells[c].par);
				sumP = (propensities[2*N] + propensities[(2*N)+1]);
				for(i = 0; i < N; i++) {
					propensities[i]   = (rep[0]*cells[c].good[i]);
					propensities[N+i] = (rep[1]*cells[c].med[i]);
					sumP += (propensities[i] + propensities[N+i]);
				}
							
				r = newRandomChoice(propensities, sumP); // select reaction index `r`
				if(r < N) from = 0; else if (r < (2*N)) from = 1; else if (r == (2*N)) from = 2; else if(r == (2*N)+1) from = 3;
				switch (from) {
					case 0: to = newRandomChoice(mutg, 1); break;
					case 1: to = newRandomChoice(mutm, 1); break;
					case 2: to = newRandomChoice(mutj, 1); break;
					case 3: to = newRandomChoice(mutp, 1);
				}
				switch (to) {
					case 0: cells[c].good[r % N]++; break;
					case 1: cells[c].med[r % N]++;  break;
					case 2: cells[c].junk++;        break;
					case 3: cells[c].par++;
				}
				updateCell(c);
				
				event[r]++;
				mutC[from][to]++;
				eventC++;
				
				/*
				if (0) {
					printf("T%d\tD%d\tV%d\tW%lf\t", t, sumD, sumV, sumW);
					printf("C%d\td%d\tv%d\tw%lf\t", c, cells[c].d, cells[c].v, cells[c].w);
					printf("M[");
					for(i = 0; i < N; i++) printf("%d ", cells[c].good[i]);
					for(i = 0; i < N; i++) printf("%d ", cells[c].med[i]);
					printf("%d %d]\t", cells[c].junk, cells[c].par);
					printf("P["); for(i = 0; i < ((2*N)+2); i++) printf("%lf ", propensities[i]); printf("]");
					printf("\tP%lf\t%d: %d -> %d\n", sumP, r, from, to);
					printf("\n");
					fflush(stdout);
				}
				*/
			}
			
			{ // Splitting
				if (cells[c].v >= sd) {
					int o = c;
					while (o == c) o = randl(PS); // Select offspring daughter cell `o`
					splitCell(c, o);
					splitCount++;
				}
			}

			{ // Fusion
				if ((fusion > 0) && ((PS-sumD) > 1) && (randd() < fusion)) {
					int p = -9, q = -9;
					double wp, wq; 
					
					//printf("%d %le ", t, sumW);
					switch (fdf) {
						case(-1): // fitness-dependent fusion, inversely proportional
							p = selectCellInv(sumW);
							q = selectCellInv(sumW); // Select fusion partners `p` and `q`; dead cells are never selected
							while (q == p) q = selectCellInv(sumW); // Make sure that two distinct cells are selected
							break;
						case(1): // fitness-dependent fusion, directly proportional
							p = selectCell(sumW);
							q = selectCell(sumW); // Select fusion partners `p` and `q`; dead cells are never selected
							while (q == p) q = selectCell(sumW); // Make sure that two distinct cells are selected
							break;
						case(0): // fitness-independent (random) fusion
							p = randl(PS);
							q = randl(PS); // Select fusion partners `p` and `q`
							while (cells[p].d) p = randl(PS); // Dead cells are never fused
							while ((q == p) || (cells[q].d)) q = randl(PS); // Make sure that two distinct cells are selected and the partner is not dead either
							break;
						case(2): { // low-fitness fusion
							int x = 0, fit[MAXC] = {0};
							for(i = 0; i < PS; i++) if(!(cells[i].d) && (cells[i].w <= limit)) { // collect cells matching the fusion condition
								fit[x] = i;
								x++;
							}
							if(x >= 2) { // only find fusion partners if there are at least two cells metching the condition
								p = fit[randl(x)];
								q = fit[randl(x)]; // Select fusion partners `p` and `q`
								while ((cells[p].d) || (cells[p].w > limit)) p = fit[randl(PS)]; // Dead cells are never fused
								while ((q == p) || (cells[q].d) || (cells[q].w > limit)) q = fit[randl(PS)]; // Make sure that two distinct cells are selected and the partner is not dead either
							}
							break;
						}
					}

					if((p >= 0) && (q >= 0)) { // only perform fusion if appropriate partners were found
						fusionCount++;
						wp = cells[p].w;
						wq = cells[q].w;
						for(i = 0; i < N; i++) {
							cells[p].good[i] += cells[q].good[i];
							cells[p].med[i]  += cells[q].med[i];
						}
						cells[p].junk += cells[q].junk;
						cells[p].par  += cells[q].par;
						splitCell(p, q); // Splitting is always happening after fusion
						
						//if((cells[p].w+cells[q].w)/(wp+wq) > 2.5) printf("%d\n", t);
					}
				} 
			}
			
			t++;

		} // if there are living cells
		
	} // time		
	
	
	if (0) { // Event statistics
	  printf("expected mutation probabilities:\n");
		for(i = 0; i < 4; i++) printf("%lf\t", mutg[i]); printf("\n");
		for(i = 0; i < 4; i++) printf("%lf\t", mutm[i]); printf("\n");
		for(i = 0; i < 4; i++) printf("%lf\t", mutj[i]); printf("\n");
		for(i = 0; i < 4; i++) printf("%lf\t", mutp[i]); printf("\n");
		printf("events: "); for(i = 0; i < (2*N)+2; i++) printf("%d\t", event[i]); printf("\n");
		printf("actual mutation probabilities:\n");
		for(i = 0; i < 4; i++) {
			for(j = 0; j < 4; j++) printf("%lf\t",(double)mutC[i][j]/(double)eventC);
			printf("\n");
		}
		printf("split  count (last period): %d\n",  splitCount);
		printf("fusion count (last period): %d\n", fusionCount);
	}
		
	fclose(file);

	//printf("Ready.\n");
	
	return(0);
}