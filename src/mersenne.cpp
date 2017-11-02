/**************************   mersenne.cpp   **********************************
 * Author:        Agner Fog
 * Date created:  2001
 * Last modified: 2008-11-16
 * Project:       randomc.h
 * Platform:      Any C++
 * Description:
 * Random Number generator of type 'Mersenne Twister'
 *
 * This random number generator is described in the article by
 * M. Matsumoto & T. Nishimura, in:
 * ACM Transactions on Modeling and Computer Simulation,
 * vol. 8, no. 1, 1998, pp. 3-30.
 * Details on the initialization scheme can be found at
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 *
 * Further documentation:
 * The file ran-instructions.pdf contains further documentation and
 * instructions.
 *
 * Copyright 2001-2008 by Agner Fog.
 * GNU General Public License http://www.gnu.org/licenses/gpl.html
 *******************************************************************************/

#include "randomc.h"
#include <cmath>
#include <iostream>



void CRandomMersenne::Init0(int seed) {
    // Seed generator
    const uint32_t factor = 1812433253UL;
    mt[0]= seed;
    for (mti=1; mti < MERS_N; mti++) {
        mt[mti] = (factor * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
    }
}

void CRandomMersenne::RandomInit(int seed) {
    // Initialize and seed
    Init0(seed);

    // Randomize some more
    for (int i = 0; i < 37; i++) BRandom();
}


void CRandomMersenne::RandomInitByArray(int const seeds[], int NumSeeds) {
    // Seed by more than 32 bits
    int i, j, k;

    // Initialize
    Init0(19650218);

    if (NumSeeds <= 0) return;

    // Randomize mt[] using whole seeds[] array
    i = 1;  j = 0;
    k = (MERS_N > NumSeeds ? MERS_N : NumSeeds);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + (uint32_t)seeds[j] + j;
        i++; j++;
        if (i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}
        if (j >= NumSeeds) j=0;}
    for (k = MERS_N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i;
        if (++i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}}
    mt[0] = 0x80000000UL;  // MSB is 1; assuring non-zero initial array

    // Randomize some more
    mti = 0;
    for (int i = 0; i <= MERS_N; i++) BRandom();
}


uint32_t CRandomMersenne::BRandom() {
    // Generate 32 random bits
    uint32_t y;

    if (mti >= MERS_N) {
        // Generate MERS_N words at one time
        const uint32_t LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
        const uint32_t UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
        static const uint32_t mag01[2] = {0, MERS_A};

        int kk;
        for (kk=0; kk < MERS_N-MERS_M; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
            mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

        for (; kk < MERS_N-1; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
            mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}

        y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
        mti = 0;
    }
    y = mt[mti++];

    // Tempering (May be omitted):
    y ^=  y >> MERS_U;
    y ^= (y << MERS_S) & MERS_B;
    y ^= (y << MERS_T) & MERS_C;
    y ^=  y >> MERS_L;

    return y;
}


double CRandomMersenne::Random() {
    // Output random float number in the interval 0 <= x < 1
    // Multiply by 2^(-32)
    return (double)BRandom() * (1./(65536.*65536.));
}


int CRandomMersenne::IRandom(int min, int max) {
    // Output random integer in the interval min <= x <= max
    // Relative error on frequencies < 2^-32
    if (max <= min) {
        if (max == min) return min; else return 0x80000000;
    }
    // Multiply interval with random and truncate
    int r = int((double)(uint32_t)(max - min + 1) * Random() + min);
    if (r > max) r = max;
    return r;
}


int CRandomMersenne::IRandomX(int min, int max) {
    // Output random integer in the interval min <= x <= max
    // Each output value has exactly the same probability.
    // This is obtained by rejecting certain bit values so that the number
    // of possible bit values is divisible by the interval length
    if (max <= min) {
        if (max == min) return min; else return 0x80000000;
    }
#ifdef  INT64_SUPPORTED
    // 64 bit integers available. Use multiply and shift method
    uint32_t interval;                    // Length of interval
    uint64_t longran;                     // Random bits * interval
    uint32_t iran;                        // Longran / 2^32
    uint32_t remainder;                   // Longran % 2^32

    interval = uint32_t(max - min + 1);
    if (interval != LastInterval) {
        // Interval length has changed. Must calculate rejection limit
        // Reject when remainder >= 2^32 / interval * interval
        // RLimit will be 0 if interval is a power of 2. No rejection then
        RLimit = uint32_t(((uint64_t)1 << 32) / interval) * interval - 1;
        LastInterval = interval;
    }
    do { // Rejection loop
        longran  = (uint64_t)BRandom() * interval;
        iran = (uint32_t)(longran >> 32);
        remainder = (uint32_t)longran;
    } while (remainder > RLimit);
    // Convert back to signed and return result
    return (int32_t)iran + min;

#else
    // 64 bit integers not available. Use modulo method
    uint32_t interval;                    // Length of interval
    uint32_t bran;                        // Random bits
    uint32_t iran;                        // bran / interval
    uint32_t remainder;                   // bran % interval

    interval = uint32_t(max - min + 1);
    if (interval != LastInterval) {
        // Interval length has changed. Must calculate rejection limit
        // Reject when iran = 2^32 / interval
        // We can't make 2^32 so we use 2^32-1 and correct afterwards
        RLimit = (uint32_t)0xFFFFFFFF / interval;
        if ((uint32_t)0xFFFFFFFF % interval == interval - 1) RLimit++;
    }
    do { // Rejection loop
        bran = BRandom();
        iran = bran / interval;
        remainder = bran % interval;
    } while (iran >= RLimit);
    // Convert back to signed and return result
    return (int32_t)remainder + min;

#endif
}

double CRandomMersenne::normal(double m, double s)
{
    double normal_x1;                   // first random coordinate (normal_x2 is member of class)
    double w;                           // radius

    if (normal_x2_valid) {              // we have a valid result from last call
        normal_x2_valid = 0;
        return normal_x2 * s + m;
    }

    // make two normally distributed variates by Box-Muller transformation
    do {
        normal_x1 = 2. * Random() - 1.;
        normal_x2 = 2. * Random() - 1.;
        w = normal_x1*normal_x1 + normal_x2*normal_x2;
    } while (w >= 1. || w < 1E-30);

    w = std::sqrt(std::log(w)*(-2./w));
    normal_x1 *= w;  normal_x2 *= w;    // normal_x1 and normal_x2 are independent normally distributed variates
    normal_x2_valid = 1;                // save normal_x2 for next call
    return normal_x1 * s + m;
}

/***********************************************************************
 Poisson distribution
 ***********************************************************************/
int CRandomMersenne::Poisson (double L) {
    /*
     This function generates a random variate with the poisson distribution.

     Uses inversion by chop-down method for L < 17, and ratio-of-uniforms
     method for L >= 17.

     For L < 1.E-6 numerical inaccuracy is avoided by direct calculation.
     */

    //------------------------------------------------------------------
    //                 choose method
    //------------------------------------------------------------------
    if (L < 17) {
        if (L < 1.E-6) {
            if (L == 0) return 0;
            if (L < 0) std::cout << "Parameter negative in poisson function\n";

            //--------------------------------------------------------------
            // calculate probabilities
            //--------------------------------------------------------------
            // For extremely small L we calculate the probabilities of x = 1
            // and x = 2 (ignoring higher x). The reason for using this
            // method is to prevent numerical inaccuracies in other methods.
            //--------------------------------------------------------------
            return PoissonLow(L);
        }
        else {
            //--------------------------------------------------------------
            // inversion method
            //--------------------------------------------------------------
            // The computation time for this method grows with L.
            // Gives overflow for L > 80
            //--------------------------------------------------------------
            return PoissonInver(L);
        }
    }
    else {
        if (L > 2.E9) std::cout << "Parameter too big in poisson function\n";

        //----------------------------------------------------------------
        // ratio-of-uniforms method
        //----------------------------------------------------------------
        // The computation time for this method does not depend on L.
        // Use where other methods would be slower.
        //----------------------------------------------------------------
        return PoissonRatioUniforms(L);
    }
}


/***********************************************************************
 Subfunctions used by poisson
 ***********************************************************************/
int CRandomMersenne::PoissonLow(double L) {
    /*
     This subfunction generates a random variate with the poisson
     distribution for extremely low values of L.

     The method is a simple calculation of the probabilities of x = 1
     and x = 2. Higher values are ignored.

     The reason for using this method is to avoid the numerical inaccuracies
     in other methods.
     */
    double d, r;
    d = sqrt(L);
    if (Random() >= d) return 0;
    r = Random() * d;
    if (r > L * (1.-L)) return 0;
    if (r > 0.5 * L*L * (1.-L)) return 1;
    return 2;
}


int CRandomMersenne::PoissonInver(double L) {
    /*
     This subfunction generates a random variate with the poisson
     distribution using inversion by the chop down method (PIN).

     Execution time grows with L. Gives overflow for L > 80.

     The value of bound must be adjusted to the maximal value of L.
     */
    const int bound = 130;              // safety bound. Must be > L + 8*sqrt(L).
    double r;                           // uniform random number
    double f;                           // function value
    int x;                          // return value

    if (L != pois_L_last) {             // set up
        pois_L_last = L;
        pois_f0 = exp(-L);               // f(0) = probability of x=0
    }
    while (1) {
        r = Random();  x = 0;  f = pois_f0;
        do {                             // recursive calculation: f(x) = f(x-1) * L / x
            r -= f;
            if (r <= 0) return x;
            x++;
            f *= L;
            r *= x;                       // instead of f /= x
        }
        while (x <= bound);
    }
}

/***********************************************************************
 Log factorial function
 ***********************************************************************/
double CRandomMersenne::LnFac(int n) {
    // log factorial function. gives natural logarithm of n!

    // define constants
    static const double                 // coefficients in Stirling approximation
    C0 =  0.918938533204672722,      // ln(sqrt(2*pi))
    C1 =  1./12.,
    C3 = -1./360.;
    // C5 =  1./1260.,                  // use r^5 term if FAK_LEN < 50
    // C7 = -1./1680.;                  // use r^7 term if FAK_LEN < 20
    // static variables
    static double fac_table[FAK_LEN];   // table of ln(n!):
    static int initialized = 0;         // remember if fac_table has been initialized

    if (n < FAK_LEN) {
        if (n <= 1) {
            if (n < 0) std::cout << "Parameter negative in LnFac function\n";
            return 0;
        }
        if (!initialized) {              // first time. Must initialize table
            // make table of ln(n!)
            double sum = fac_table[0] = 0.;
            for (int i=1; i<FAK_LEN; i++) {
                sum += log(double(i));
                fac_table[i] = sum;
            }
            initialized = 1;
        }
        return fac_table[n];
    }
    // not found in table. use Stirling approximation
    double  n1, r;
    n1 = n;  r  = 1. / n1;
    return (n1 + 0.5)*log(n1) - n1 + C0 + r*(C1 + r*r*C3);
}


int CRandomMersenne::PoissonRatioUniforms(double L) {
    /*
     This subfunction generates a random variate with the poisson
     distribution using the ratio-of-uniforms rejection method (PRUAt).

     Execution time does not depend on L, except that it matters whether L
     is within the range where ln(n!) is tabulated.

     Reference: E. Stadlober: "The ratio of uniforms approach for generating
     discrete random variates". Journal of Computational and Applied Mathematics,
     vol. 31, no. 1, 1990, pp. 181-189.
     */
    double u;                                          // uniform random
    double lf;                                         // ln(f(x))
    double x;                                          // real sample
    int k;                                             // integer sample

    if (pois_L_last != L) {
        pois_L_last = L;                                // Set-up
        pois_a = L + 0.5;                               // hat center
        int mode = (int)L;                              // mode
        pois_g  = log(L);
        pois_f0 = mode * pois_g - LnFac(mode);          // value at mode
        pois_h = sqrt(SHAT1 * (L+0.5)) + SHAT2;         // hat width
        pois_bound = (int)(pois_a + 6.0 * pois_h);      // safety-bound
    }
    while(1) {
        u = Random();
        if (u == 0) continue;                           // avoid division by 0
        x = pois_a + pois_h * (Random() - 0.5) / u;
        if (x < 0 || x >= pois_bound) continue;         // reject if outside valid range
        k = (int)(x);
        lf = k * pois_g - LnFac(k) - pois_f0;
        if (lf >= u * (4.0 - u) - 3.0) break;           // quick acceptance
        if (u * (u - lf) > 1.0) continue;               // quick rejection
        if (2.0 * log(u) <= lf) break;                  // final acceptance
    }
    return k;
}


CRandomMersenne rndgen(5); //the one random number generator

void set_seed(int seed)
{
    rndgen.RandomInit(seed);
}

double uniform()
{
    return rndgen.Random();
}

int random_number(int n)
{
    return rndgen.IRandom(0,n-1);
}

double normal(double m, double s)
{
    return rndgen.normal(m,s);
}

double poisson(double lambda)
{
    return rndgen.Poisson(lambda);
}






