#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <assert.h>
#include <math.h>

#define max_phsp_num 27
#define right(philosopher_id) ((philosopher_id+1)%num_phsp)
#define left(philosopher_id) philosopher_id
#define rightNeighbor(philosopher_id) ((philosopher_id == num_phsp-1)? 0 : philosopher_id+1)
#define leftNeighbor(philosopher_id) ((philosopher_id == 0)? num_phsp-1 : philosopher_id-1)

// states of philosophers 
typedef enum {
	THINKING = 0, HUNGRY = 1, DINING = 2
} State;

//states are defined as free and busy

typedef enum {
	FREE = 0, BUSY = 1
} CHOPSTICKS;

// condition variable is each philosophers

pthread_cond_t cond[max_phsp_num];
State states_philosophers[max_phsp_num];
CHOPSTICKS states_Chopstick[max_phsp_num];

//this is needed for critical section 

pthread_mutex_t mutex;

int num_phsp;

typedef struct {
	int phspID;
	unsigned int think_time;
	unsigned int dine_time;
	int dine_counter;
	unsigned int hungryTimeInTotal;
	int *hungryTimes;
} philosopher;



void assignSituation(void *args) {
	philosopher *phspdispute = args;

	pthread_mutex_lock(&mutex);

	states_Chopstick[left(phspdispute->phspID)] = FREE;
	states_Chopstick[right(phspdispute->phspID)] = FREE;

	states_philosophers[phspdispute->phspID] = THINKING;

	pthread_cond_signal(&cond[leftNeighbor(phspdispute->phspID)]);
	pthread_cond_signal(&cond[rightNeighbor(phspdispute->phspID)]);

	pthread_mutex_unlock(&mutex);

	printf("The philosopher %d is thinking\n", phspdispute->phspID);
	usleep(phspdispute->think_time);
}

// to get chopstick
void getCS(void *args, int diningNumber) {

	philosopher *phspdispute = args;

	pthread_mutex_lock(&mutex);

	states_philosophers[phspdispute->phspID] = HUNGRY;
	printf("The philosopher  %d is hungry\n", phspdispute->phspID);
	clock_t start_time = clock();

	while (states_Chopstick[left(phspdispute->phspID)] == BUSY
			|| states_Chopstick[right(phspdispute->phspID)] == BUSY) {
		pthread_cond_wait(&cond[phspdispute->phspID], &mutex);
	}

	clock_t end_time = clock();

	int starvingDuration = end_time - start_time;

	phspdispute->hungryTimeInTotal += starvingDuration;
	phspdispute->hungryTimes[diningNumber] = starvingDuration;
	states_philosophers[phspdispute->phspID] = DINING;

	states_Chopstick[left(phspdispute->phspID)] = BUSY;
	states_Chopstick[right(phspdispute->phspID)] = BUSY;

	pthread_mutex_unlock(&mutex);

	printf("The philosopher  %d is dining\n", phspdispute->phspID);
	usleep(phspdispute->dine_time);
}

void *lifeTimeOfAPshp(void *args) {

	philosopher *phspdispute = args;
	for (int i = 0; i < phspdispute->dine_counter; i++) {
		getCS(args, i);
		assignSituation(args);
	}

	return 0;
}
// uniform random generation method

unsigned int generateUniformRandomValue(int lowerBound, int upperBound) {
	return (unsigned int) (rand() % (upperBound + 1 - lowerBound) + lowerBound);
}

// exponential random generation method parts taken from internet.

unsigned int generateExponentialRandomValue(double mean, int lowerBound,
		int upperBound) {

	double randomValue = -log((double) rand() / (double) RAND_MAX) / (1 / mean);

	while (randomValue < lowerBound || randomValue > upperBound) {
		randomValue = -log((double) rand() / (double) RAND_MAX) / (1 / mean);
	}

	return (unsigned int) randomValue;
}

// to compute standard deviation
double std(int array[], int dine_count)
{
    double sum = 0.0, mean, standardDeviation = 0.0;
	for (int i = 0; i < dine_count; i++)
		sum += array[i];
	mean = sum / dine_count;

	for (int i = 0; i < dine_count; i++)
		standardDeviation += pow(array[i] - mean, 2);

	return sqrt(standardDeviation / dine_count);
}

int main(int argc, char **argv) {

	num_phsp = atoi(argv[1]);

	assert(num_phsp <= max_phsp_num
);

	const int min_think = atoi(argv[2]) * 100000; //min_think time
	const int max_think = atoi(argv[3]) * 100000; //max_think time
	const int min_dine = atoi(argv[4]) * 100000; //min_dine time
	const int max_dine = atoi(argv[5]) * 100000; //max_dine time
	const char * const distributionType = argv[6]; // distribution can be either “uniform” or “exponential”
	const int counter = atoi(argv[7]);

	printf("number of philosophers : %d\n", num_phsp);
	printf("minimum thinking time : is %d ms\n", min_think / 100); // in milisecond.
	printf("maximum thinking time is : %d ms\n", max_think / 100);
	printf("minimum dining time is: %d ms\n", min_dine / 100);
	printf("maximum dining time is: %d ms\n", max_dine / 100);
	printf("type of the distribution is: %s\n", distributionType);
	printf("dine count value is : %d\n\n", counter);

	// number of phsp maximum 27 and odd numbers
	assert(num_phsp <= 27);
	assert(num_phsp % 2 == 1);

	pthread_t phspThread[num_phsp];
	philosopher *phspArray[num_phsp];
	int error;

	double mean = (min_dine + max_dine) / 2;

	for (int i = 0; i < num_phsp; i++) {
		states_philosophers[i] = THINKING;
		states_Chopstick[i] = FREE;
		error = pthread_cond_init(&cond[i], NULL); // no error, return zero

		assert(!error);
	}

	//the threads for each philosopher

	for (int index = 0; index < num_phsp; index++) {

		phspArray[index] = malloc(sizeof *phspArray[0]); //allocating space for phsp array.

		phspArray[index]->phspID = index; //assigning the properties to philosoperhers.
		phspArray[index]->dine_counter = counter;

		if (strcmp("uniform", distributionType) == 0) {
			phspArray[index]->dine_time = generateUniformRandomValue(min_dine,
					max_dine); //generate random  dine_time
			phspArray[index]->think_time = generateUniformRandomValue(min_dine,
					max_dine); //generate random  think_time
		} else if (strcmp("exponential", distributionType) == 0) {
			phspArray[index]->dine_time = generateExponentialRandomValue(mean,
					min_dine, max_dine); //generate  random  dine_time
			phspArray[index]->think_time = generateExponentialRandomValue(mean,
					min_dine, max_dine); //generate random  think_time
		} else
			assert(0);

		phspArray[index]->hungryTimeInTotal = 0;
		phspArray[index]->hungryTimes = malloc(sizeof(int) * counter);

		error = pthread_create(&(phspThread[index]), NULL, lifeTimeOfAPshp,
				phspArray[index]);  // CREATION OF A NEW PHILOSPHER.

		assert(!error);
	}

	pthread_mutex_init(&mutex, NULL);

	// wait philosophers to be death

	for (int i = 0; i < num_phsp; i++)
		pthread_join(phspThread[i], NULL);

	double averageHungryState = 0;
	for (int i = 0; i < num_phsp; i++) {
		printf("The philosopher  %d duration of hungry state = %u\n",
				phspArray[i]->phspID, phspArray[i]->hungryTimeInTotal);
		averageHungryState += phspArray[i]->hungryTimeInTotal;
	}
	averageHungryState /= num_phsp;

	printf("Avearage hungry state %f \n\n", averageHungryState);

	printf("\nSTANDARD DEVIATION OF HUNGRY STATES OF PHILOSOPHERS \n\n");

	for (int i = 0; i < num_phsp; i++) {
		for (int j = 0; j < counter; j++) {
			printf(
					"The time for philosopher : %d spent for hungry time %d is %d\n",
					phspArray[i]->phspID, j, phspArray[i]->hungryTimes[j]);
		}
		printf("The SD (standard deviation) for philosopher %d is %.1f\n",
				phspArray[i]->phspID,
				std(phspArray[i]->hungryTimes, counter));
	}

	return 0;
}