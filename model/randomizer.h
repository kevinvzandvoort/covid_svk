#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/*
 * Class Randomizer
 * Random number generator using GSL
 * Distributions implemented: binomial and uniform
 */
class Randomizer{
public:
    int seed;
    gsl_rng* r;
    
    Randomizer(int s)
        : seed(s), r(gsl_rng_alloc(gsl_rng_mt19937))
    {
        Reset();
    }
    
    ~Randomizer(){
        gsl_rng_free(r);
    }

    void Reset(){
        gsl_rng_set(r, seed);
    }
    
    unsigned int Binomial(unsigned int n, double p){
        if (p <= 0) return 0;
        return gsl_ran_binomial(r, p, n);
    }

    int UniformInt(int min, int max){
        return min + gsl_rng_uniform(r) * (max - min);
    }
};
