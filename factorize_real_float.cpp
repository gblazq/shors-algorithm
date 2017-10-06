#include <iostream>
#include <cmath>
#include <vector>
#include <gmpxx.h>
#include <cstdlib>
#include <fftw3.h>
#include <time.h>
#include <random>
#include <chrono>

unsigned int N; // Number to factorize
unsigned int M; // Some non-prime factor of N to factorize
unsigned int n;
unsigned long long Q;
const unsigned long MAX_N = std::pow(2, 15) - 1; // Maximum value of N such that the memory needed for the Fourier transform is 16GB

std::vector<unsigned int> prime_factors; // Vector of prime factors of N
bool all_factors_prime = false;
std::vector<unsigned int> factors; // Vector of two factors of M

fftwf_complex *finalState;
fftwf_plan p;

// Initialize the random number generator
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);

// ************************************************************************

int is_prime(unsigned int m){
    /*
    Miller-Rabin primality test
    */
    mpz_t m_mpz;
    mpz_init_set_ui(m_mpz, m);
    return mpz_probab_prime_p(m_mpz, 30);
}

unsigned int gcd(unsigned int m, unsigned int n) {
    /*
    Greatest common divisor for small numbers using Euclid's algorithm
    */
	unsigned int r;

	while ( n != 0 ) {
		r = m % n;
		m = n;
		n = r;
	}

	return m;
}

unsigned int lcm(unsigned int m, unsigned int n) {
    /*
    Least common multiple
    */
    return m*n/gcd(m, n);
}

unsigned int mod_exp(unsigned int base, unsigned long long exp, unsigned int mod) {
    /*
    Modular exponentiation
    */
    unsigned int result = 1;

    base = base % mod;
    while ( exp > 0 ) {
        if ( (exp % 2) == 1) {
            result = ( result*base ) % mod;
        }
        exp >>= 1;
        base = ( base*base ) % mod;
    }
    return result;
}

unsigned int find_period(unsigned int y, unsigned int M) {
    unsigned int* modexpn = new unsigned int [n]; // Array for the modular exponentiations of powers of 2
    unsigned int* modexpq = new unsigned int [Q]; // Array for all the modular exponentiations
    float* modexpq_float = new float [Q];
    unsigned long power; // Temporary variable for powers
    unsigned long m; // Random number for measuring the target
    unsigned int period = 0;
    unsigned long c = 0;
    unsigned int ft_length = Q/2+1;
    unsigned long temp_mod;

    std::cout << "\n\tComputing the modular exponentiations";
    // Compute the modular exponentiations of powers of 2
    for (unsigned int i = 0; i < n; ++i) {
        modexpn[i] = mod_exp(y, pow(2, i), M);
    }

    // Compute all the Q = 2**n modular exponentiations
    // 		Initialize the modexpq array to ones
    std::fill_n(modexpq, Q, 1);

    // 		Multiply mod_exp(y, 2**i, M) if q_i == 1 (in base 2) and copy it where needed
    for (unsigned long i = 0; i < n; ++i) {
        std::cout << '.';
        std::cout.flush();
        power = pow(2, i);
        for (unsigned long k = power; k < 2*power; ++k){
            temp_mod = modexpq[k];
            temp_mod *= modexpn[i];
            temp_mod %= M;
            modexpq[k] = temp_mod;
        }
        for (unsigned long l = 2*power; l < Q/2; l *=2){
            memcpy(&modexpq[l], &modexpq[0], sizeof(int)*l);
        }
//        for (unsigned long j = 3*power; j < Q; j += 2*power){
//            memcpy(&modexpq[j], &modexpq[power], sizeof(int)*power);
//        }
    }
    delete[] modexpn;

    // Measure the target and collapse the state
    //		Select a random result
    m = rand() % Q;
    m = modexpq[m];
    std::cout << "\n\tMeasurement of the target is " << m << ". The state has collapsed." << std::endl;

    // 		Create a complex vector with all the non-normalized coefficients of the states

    // Perform the QFT (Fourier transform of the state)
    // 		Create the vector of the state
    for (unsigned long long i = 0; i < Q; ++i) {
        modexpq_float[i] = (modexpq[i] == m);
    }
    delete[] modexpq;

    //		Compute the FFT of the state
    finalState = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*ft_length);
    p = fftwf_plan_dft_r2c_1d(Q, modexpq_float, finalState, FFTW_BACKWARD | FFTW_ESTIMATE);

    std::cout << "\tComputing the QFT." << std::endl;
    std::cout.flush();
    fftwf_execute(p);

    // Compute the probabilities
    for (unsigned long i = 0; i < ft_length; ++i) {
        auto it = &finalState[0][0] + i;
        *it = finalState[i][0]*finalState[i][0] + finalState[i][1]*finalState[i][1];
    }

    std::discrete_distribution<unsigned long> random_c (&finalState[0][0], &finalState[0][0]+ft_length);
    fftwf_free(finalState);

    while(c == 0){
        c = random_c(generator);
        if (rand() % 2){
            c = ft_length - c;
        }
    }

    std::cout << "\tMeasurement of the register is " << c << "." << std::endl;
    unsigned long pn1 = 1, pn2 = 0, qn1 = 0, qn2 = 0;
    double chi = (double)c/Q;
    long an = floor(chi);
    double chi_n = chi - an;
    double invchi;
    unsigned long pn = an;
    unsigned long qn = 1;
    while (chi_n > 1e-4) {
        pn2 = pn1;
        pn1 = pn;
        qn2 = qn1;
        qn1 = qn;
        invchi = 1/chi_n;
        an = floor(invchi);
        chi_n = invchi - an;
        pn = an*pn1 + pn2;
        if(an*qn1 + qn2 > M){
            break;
        }
        qn = an*qn1 + qn2;
    }
    for(unsigned int i = qn; i<M; i+=qn){
        if(mod_exp(y, i, M) == 1) {
            period = i;
            break;
        }
    }
    return period;
}

std::vector<unsigned int> factorize(unsigned int M){
    /*
    Factorize M in two factors using Shor's algorithm
    */
    std::vector<unsigned int> factors;
    std::cout << "\nFactoring " << M << ".\n" << std::endl;
    // Check if M is even
    if ( M%2 == 0){
        std::cout << '\t' << M << " is even. Factors found.\n" << std::endl;
        factors.push_back(2);
        factors.push_back(M/2);

        return factors;
    } else if (M%3 == 0){
        std::cout << '\t' << M << " is divisible by 3. Factors found.\n" << std::endl;
        factors.push_back(3);
        factors.push_back(M/3);
        return factors;
    }
    n = 2*ceil(log2(M+1)); // Number of bits of M**2
    Q = pow(2, n); // Number of qubits needed
    bool valid_period = false;
    unsigned int period;
    unsigned int y;

    // Initialize a uniform distribution of integers
    std::uniform_int_distribution<long> distribution(2,M-1);

    while(!valid_period) {
        // Find some random integer less than M
        y = distribution(generator);
        std::cout << "\ty = " << y;

        // Check if we were lucky
        if (gcd(y, M) != 1) {
            std::cout << ". gdc(" << y << ", " << M << ") != 1. Factors found.\n" << std::endl;
            factors.push_back(gcd(y, M));
            factors.push_back(M/gcd(y, M));
            return factors;
        } else {
            //If we weren't, find the period of y**q mod M
            std::cout << ". Finding the period of f(x) = " << y << "**x mod " << M << '.';
            period = find_period(y, M);
        }

        // Check if the period is valid
        if(period == 0){
            std::cout << "\tPeriod not found. Trying again.\n" << std::endl;
            continue;
        }
        std::cout << "\tPeriod of f(x) is " << period;
        if (period%2 == 1) {
            std::cout << ". It is odd. Trying again." << std::endl;
        } else if (mod_exp(y, period/2, M) == M-1) {
            std::cout << ". It is not valid. " << y << "**" << period/2 << " = -1 mod " << M << ". Trying again." << std::endl;
        } else {
            std::cout << ". Factors found." << std::endl;
            valid_period = true;
        }
    }

    // Once we have found the valid period, find the factors
    mpz_class power_t;
    mpz_class factor_t;
    mpz_class poweradd_t;
    mpz_ui_pow_ui( power_t.get_mpz_t(), y, period/2);
    mpz_add_ui(poweradd_t.get_mpz_t(), power_t.get_mpz_t(), 1);
    mpz_gcd_ui( factor_t.get_mpz_t(), poweradd_t.get_mpz_t(), M);
    factors.push_back(factor_t.get_ui());
    factors.push_back(M/factor_t.get_ui());

    return factors;
}

int main(int argc, char *argv[]) {
    fftwf_init_threads();
    fftwf_plan_with_nthreads(4);

    // Get N
    N = std::atoi(argv[1]);

    // Check that N is not too high. For higher values we would need at least 16GB of memory
    if (N > MAX_N) {
        std::cout << "N is too high.";
        return 0;
    }

    // Check whether N is prime
    if (is_prime(N)){
        std::cout << N << " is prime." << std::endl;
        return 0;
    }

    prime_factors.push_back(N);
    while(!all_factors_prime) {
        all_factors_prime = true;
        for (std::vector<unsigned int>::iterator it = prime_factors.begin(); it != prime_factors.end(); ++it){
            // Check if the factors found until now are 1 or prime
            if (is_prime(*it) == 0) {
                all_factors_prime = false;
                M = *it;
                prime_factors.erase(it);
                // If one of the factors is composite, factorize it
                factors = factorize(M);
                std::sort(factors.begin(), factors.end());
                std::cout << "\t" << M << " = " << factors[0] << 'x' << factors[1] << std::endl << std::endl;
                // Add the factors of M to the factors of N
                for (std::vector<unsigned int>::iterator it2 = factors.begin(); it2 != factors.end(); ++it2) {
                    if(*it2 != 1){
                        prime_factors.push_back(*it2);
                    }
                }
                break;
            }
        }
    }

    // Sort the factors in ascending order
    std::sort(prime_factors.begin(), prime_factors.end());

    // Print the factorization
    std::cout << N << " = ";
    for (std::vector<unsigned int>::iterator it = prime_factors.begin(); it != prime_factors.end(); ++it){
        if (it == prime_factors.begin()){
            std::cout << *it;
        } else {
            std::cout << 'x' << *it ;
        }
    }
    std::cout << std::endl;

    return 0;
}
