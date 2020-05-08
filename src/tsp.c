#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>

typedef struct coord {
  float x;
  float y;
} coord;

float dist(coord a, coord b) {
  return (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y);
}

float pathlen(coord *coords, unsigned short *path, size_t n_cities) {
  float len = 0.0;
  for (size_t i = 1; i < n_cities; ++i)
    len += dist(coords[path[i]], coords[path[i - 1]]);
  return len;
}

void swap(unsigned short *a, unsigned short *b) {
  unsigned short tmp = *a;
  *a = *b;
  *b = tmp;
}

void shuffle(unsigned short *arr, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    size_t j = rand() % (i + 1);
    swap(&arr[i], &arr[j]);
  }
}

void InitPath(size_t n_cities, unsigned short* path) {
  for (size_t i = 0; i < n_cities; ++i)
    path[i] = i;
  shuffle(path, n_cities);
}

void mutate(unsigned short* path, size_t n_cities) {
  size_t i = rand() % (n_cities);
  size_t j = rand() % (n_cities);
  swap(&path[i], &path[j]);
}

bool IsInPath(unsigned short *path, size_t city, size_t i_max) {
  for (size_t i = 0; i < i_max; ++i) {
    if (path[i] == city)
      return true;
  }
  return false;
}

void breed(unsigned short* a, unsigned short* b, coord* coords, size_t n_cities, unsigned short* child) {
  size_t next_city = ((rand() % 2) == 1 ? b[0] : a[0]);
  child[0] = next_city;
  for (size_t i = 1; i < n_cities; ++i) {
    // Choose the next city from either parent that is closer to the current city
    float dist_a = dist(coords[a[i]], coords[child[i - 1]]);
    float dist_b = dist(coords[b[i]], coords[child[i - 1]]);
    next_city = (dist_a < dist_b ? a[i] : b[i]);

    // If the chosen city has already been visited, choose the city from the other parent
    if (IsInPath(child, next_city, i))
      next_city = (dist_a < dist_b ? b[i] : a[i]);

    if (IsInPath(child, next_city, i)) {
      // If this city has *also* been visited, try the next city in
      // the chain from a random parnet until an unvisited city is found
      size_t tmp_index = i;
      while (IsInPath(child, next_city, i)) {
    	++tmp_index;

    	// Try a randomly chosen parent
    	size_t random_parent = rand() % 2;
    	next_city = (random_parent == 1 ? b[tmp_index] : a[tmp_index]);

	// If this city has been visited, try the other parent
    	if (IsInPath(child, next_city, i))
    	  next_city = (random_parent == 1 ? a[tmp_index] : b[tmp_index]);
      }
    }
    child[i] = next_city;
  }
}

int CountLines(char* filename) {
  FILE *fp = fopen(filename, "r");
  int lines = 0;
  while(!feof(fp)) {
    int ch = fgetc(fp);
    if(ch == '\n') {
      lines++;
    }
  }
  fclose(fp);
  return lines;
}

coord* ReadCoords(char* filename, size_t n_cities) {
  coord *coords = malloc(n_cities * sizeof(coord));
  FILE *fp;
  fp = fopen (filename, "r");
  int index = 0;
  while(fscanf(fp, "%f %f", &coords[index].x, &coords[index].y) == 2)
    index++;
  fclose(fp);
  return coords;
}

float rand_p() {
  return (float)rand() / (float)RAND_MAX;
}

void MostFit(unsigned short* pops[], coord *coords, size_t pop_size, size_t n_cities) {
  size_t best_index = 0;
  float best_fitness = FLT_MAX;
  for (size_t i = 0; i < pop_size; ++i) {
    float fitness = pathlen(coords, pops[i], n_cities);
    if (fitness < best_fitness) {
      best_fitness = fitness;
      best_index = i;
    }
  }

  printf("%s", "best path: ");
  for (size_t i = 0; i < n_cities; ++i)
    printf("%d", pops[best_index][i]);
  
  printf(" fitness: %f\n", best_fitness);
}

int cmp_ptr(const void *a, const void *b)
{
    const float **left  = (const float **)a;
    const float **right = (const float **)b;

    return (**left < **right) - (**right < **left);
}

size_t * order_float(const float *a, size_t n) {
    const float **pointers = malloc(n * sizeof(const float *));

    for (size_t i = 0; i < n; i++) pointers[i] = a + i;

    qsort(pointers, n, sizeof(const float *), cmp_ptr);

    size_t *indices = malloc(n * sizeof(size_t));

    for (size_t i = 0; i < n; i++) indices[i] = pointers[i] - a;

    free(pointers);

    return indices;
}

void selection(unsigned short* old_pops[], unsigned short* new_pops[], unsigned short pop_size, coord* coords, size_t n_cities) {
  // Calculate fitness for both populationms
  float *fits = malloc(2 * pop_size * sizeof(float));
  for (size_t i = 0; i < pop_size; ++i) {
    fits[2*i] = pathlen(coords, old_pops[i], n_cities);
    fits[2*i + 1] = pathlen(coords, new_pops[i], n_cities);
  }

  // Select the pop_size most fit individuals from the combined old and new population
  size_t* indices = order_float(fits, 2*pop_size);
  for (size_t i = 0; i < pop_size; ++i) {
    size_t pos = 2 * pop_size - (i + 1);
    bool from_old = (indices[pos] % 2) == 0;
    if (from_old) {
      for (size_t k = 0; k < n_cities; ++k) {
	old_pops[i][k] = old_pops[indices[pos]/2][k];
      }
    } else {
      for (size_t k = 0; k < n_cities; ++k) {
	old_pops[i][k] = new_pops[(indices[pos] - 1)/2][k];
      }
    }
  }

  free(fits);
  free(indices);
}

void immigration(unsigned short* pops[], size_t migration_size, size_t n_cities, size_t pop_size, coord *coords) {
  // Calculate fitness for the population
  float *fits = malloc(pop_size * sizeof(float));
  for (size_t i = 0; i < pop_size; ++i) {
    fits[i] = pathlen(coords, pops[i], n_cities);
  }
  size_t* indices = order_float(fits, pop_size);

  // Replace worst 100 individuals with random immigrants
  for (size_t i = 0; i < migration_size; ++i)
    InitPath(n_cities, pops[indices[i]]);

  free(indices);
  free(fits);
}

void check_input(float mutation_prob, size_t pop_size, float migration_prob, size_t migration_size) {
  if (mutation_prob < 0 || mutation_prob >= 1) {
    printf("Mutation probability must be between 0 and 1\n");
    exit(1);
  }
  if (migration_prob < 0 || migration_prob >= 1) {
    printf("Migration probability must be between 0 and 1\n");
    exit(1);
  }
  if (migration_size >= pop_size) {
    printf("Migration size must be less than population size\n");
    exit(1);
  }
}

int w_srand(float *weights, size_t n_weights) {
  float sum = 0.0;
  for (size_t i = 0; i < n_weights; ++i) {
    sum += 1.0/weights[i];
  }
  int max = ceil(sum);
  int rn = rand() % max;
  for (size_t i = 0; i < n_weights; ++i) {
    if (rn < 1.0/weights[i])
      return i;
    rn -= 1.0/weights[i];
  }
  exit(1); // shouldn't get here
}

void w_selection(unsigned short* old_pops[], unsigned short* new_pops[], unsigned short pop_size, coord* coords, size_t n_cities) {
  // Calculate fitness for both populationms
  float *fits = malloc(2 * pop_size * sizeof(float));
  for (size_t i = 0; i < pop_size; ++i) {
    fits[2*i] = pathlen(coords, old_pops[i], n_cities);
    fits[2*i + 1] = pathlen(coords, new_pops[i], n_cities);
  }

  for (size_t i = 0; i < pop_size; ++i) {
    size_t pos = w_srand(fits, 2 * pop_size);
    bool from_old = (pos % 2) == 0;
    if (from_old) {
      for (size_t k = 0; k < n_cities; ++k) {
	old_pops[i][k] = old_pops[pos/2][k];
      }
    } else {
      for (size_t k = 0; k < n_cities; ++k) {
	old_pops[i][k] = new_pops[(pos - 1)/2][k];
      }
    }
  }

  free(fits);
}

int main(int argc, char *argv[]) {
  if (argc < 7) {
    printf("Usage: tsp <mutation probability, float> <n generations, size_t> <population size, size_t> <migration probability, float> <migration size, size_t>\n");
    exit(1);
  }
  char* infile = argv[1];
  float mutation_prob = atof(argv[2]);
  size_t n_generations = atoi(argv[3]);
  size_t pop_size = atoi(argv[4]);
  float migration_prob = atof(argv[5]);
  size_t migration_size = atof(argv[6]);

  check_input(mutation_prob, pop_size, migration_prob, migration_size);

  // Check how many cities there are and read in their coordinates
  size_t n_cities = CountLines(infile);
  coord *coords = ReadCoords(infile, n_cities);

  // Initialize a population from random paths
  unsigned short *pops[pop_size];
  for (size_t i = 0; i < pop_size; ++i)
    pops[i] = (unsigned short*)malloc(n_cities * sizeof(unsigned short));

  srand(time(0));
  for (size_t i = 0; i < pop_size; ++i)
    InitPath(n_cities, pops[i]);

  printf("Generation 0\t");
  MostFit(pops, coords, pop_size, n_cities);

  unsigned short *new_pops[pop_size];
  for (size_t j = 0; j < pop_size; ++j)
    new_pops[j] = (unsigned short*)malloc(n_cities * sizeof(unsigned short));
  
  for (size_t i = 0; i < n_generations; ++i) {
    // Immigrants replace 100 worst fit individuals with prob migration_prob
    if (rand_p() < migration_prob) {
      immigration(pops, migration_size, n_cities, pop_size, coords);
    }

    float *fits = malloc(pop_size * sizeof(float));
    for (size_t i = 0; i < pop_size; ++i) {
      fits[i] = pathlen(coords, pops[i], n_cities);
    }

    for (size_t j = 0; j < pop_size; ++j) {
      // Cross self with a randomly selected pop
      size_t parent_a = j;
      // Can't mate with self.
      while (j == parent_a)
	//	parent_a = w_srand(fits, pop_size);
	parent_a = rand() % pop_size;
      breed(pops[parent_a], pops[j], coords, n_cities, new_pops[j]);

      // Introduce a mutation with probability mutation_p
      float mutate_p = rand_p();
      if (rand_p() < mutation_prob)
	mutate(new_pops[j], n_cities);
    }
    // Select the fit individuals to populate the next generation
    selection(pops, new_pops, pop_size, coords, n_cities);

    printf("Generation %zu\t", i + 1);
    MostFit(pops, coords, pop_size, n_cities);

    free(fits);
  }

  free(coords);
  for (size_t i = 0; i < pop_size; ++i) {
    free(pops[i]);
    free(new_pops[i]);
  }

  return 0;
}
