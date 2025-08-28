//
// Fast Finite Field Transform.
//
// (C) Copyright 2023 Stephen M. Watt.


#include <stdexcept>
#include <type_traits>
#include <initializer_list>
#include <cstdlib>
#include <climits>
#include <cinttypes>


///////////////////////////////////////////////////////////////////////////////
///
/// Misc utilities.
///
template <typename T> T min2(T a, T b) { return a < b ? a : b; }

// Reverse of j considered as integer of nb bits.
template <typename uint_t>
uint_t bitReverse(uint_t j, size_t nb) {
        uint_t r = 0;
        for ( ; nb--; j >>= 1) r = (r << 1) | (j & 1);
        return r;
}

template <typename uint_t>
bool is2pow(uint_t n) {
    bool seen1 = false;
    for ( ; n ; n >>= 1) if (n & 1) { if (seen1) return false; seen1 = true; }
    return seen1;
}

// ceil(log[2](n))   1->0, 2->1, 3..4->2, 5..8->3, etc
template <typename uint_t>
unsigned ceilLg(uint_t n) {
    size_t lastSeen = 0, numSeen = 0;
    for (unsigned i = 0; n; i++, n >>= 1)
        if (n & 1) { lastSeen = i; numSeen++; }
    if (numSeen > 1) lastSeen++;
    return lastSeen;
}



///////////////////////////////////////////////////////////////////////////////
///
/// Unsigned integer type capable of holding all products without overflow.
///
template <typename T> struct uprod_t { using type = __uint128_t; };

template <> struct uprod_t<uint8_t>  { using type = uint16_t; };
template <> struct uprod_t<uint16_t> { using type = uint32_t; };
template <> struct uprod_t<uint32_t> { using type = uint64_t; };


///////////////////////////////////////////////////////////////////////////////
///
/// Integers mod M.  
///
/// The modulus must fit in rep_t, which should be an unsigned type.
///

template <typename _rep_t, _rep_t _modulus>
class zmod_t {

    using ubig_t = typename uprod_t<_rep_t>::type;
    using srep_t = typename std::make_signed<_rep_t>::type;

    _rep_t val;    // Hide rep in case we later use Montgomery multiplication.


public:
    using              rep_t   = _rep_t;
    static const rep_t modulus = _modulus;

    static zmod_t newUnchecked(rep_t v) { zmod_t r; r.val = v; return r; }

    zmod_t(const zmod_t &y) : val(y.val) { }
    zmod_t()                : val(0) { }
    zmod_t(rep_t v)         : val(0 <= v && v < modulus ? v : v % modulus) { }

    bool operator ==(const zmod_t y) const { return val == y.val; }
    bool operator !=(const zmod_t y) const { return val != y.val; }

    zmod_t operator + (const zmod_t y) const {
        ubig_t r = (ubig_t) val + (ubig_t) y.val; // Could check wrap instead
        if (r >= modulus) r -= modulus;
        return newUnchecked(r);
    }
    zmod_t operator - (const zmod_t y) const {
        rep_t r = val;
        if (val < y.val) r += modulus;
        return newUnchecked(r - y.val);
    }
    rep_t  operator + () const { return val; }
    zmod_t operator - () const { return newUnchecked(val ? modulus - val : 0); }

    zmod_t operator * (const zmod_t y) const {
        ubig_t r = (ubig_t) val * (ubig_t) y.val;
        return newUnchecked(r % modulus);
    }
    zmod_t operator / (const zmod_t y) const {
        return (*this) * y.inv();
    }
    zmod_t inv(bool check = false) const {
        rep_t  a = val, b = modulus;
        srep_t s = 1,   t = 0;
        while (b != 0) {
            srep_t q = a/b, r = a % b;
            a = b; b = r; 
            r = s - q*t;
            s = t; t = r; 
        }
        if (check && a != 1) throw std::domain_error("Non-invertible element");
        if (s < 0) s += modulus;
        return newUnchecked(s);
    }
    // pow uses zmod_t arithmetic.
    zmod_t pow(rep_t n) const {
        if (n == 0) return 1;
        zmod_t a = *this, b = 1;
        for ( ; n > 1; n >>= 1) { if (n & 1) b = a*b; a = a*a; }
        return a*b;
    }
};

///////////////////////////////////////////////////////////////////////////////
///
/// Parameters for fft-friendly prime fields.
///

// Choose p < pow(1, bitsizeof(uint_t)) so no arithmetic overflows.
//
// p = k * 2**n + 1.  g**(2**n) = 1 mod p.  g**q != 1 for q < 2**n.
//
// For now simply require base*base < p.` We illustrate a variety below.

#ifdef BINARY_BASE
# define USE_BASE(bin, dec) (bin)
#else
# define USE_BASE(bin, dec) (dec)
#endif

template <typename _uint_t>
struct FftPrime {
    using uint_t = _uint_t;
    static const unsigned n = 0;
    static const unsigned k = 0;
    static const uint_t   p = 0;
    static const uint_t   g = 0;
    static const size_t   base = 0;
};

template <> struct FftPrime<uint64_t> {
    using uint_t = uint64_t;
    static const unsigned n = 57;
    static const unsigned k = 29;
    static const uint_t   p = 4179340454199820289;
    static const uint_t   g = 21;
    static const size_t   base = 100;
};

template <> struct FftPrime<uint32_t> {
    using uint_t = uint32_t;
    static const unsigned n = 30;
    static const unsigned k = 3;
    static const uint_t   p = 3221225473;
    static const uint_t   g = 13;
    static const size_t   base = 1 << 15;
};

template <> struct FftPrime<uint16_t> {
    using uint_t = uint16_t;
    static const unsigned n = 13;
    static const unsigned k = 5;
    static const uint_t   p = 40961;
    static const uint_t   g = 0xc;
    static const size_t   base = 1000000000;
};


///////////////////////////////////////////////////////////////////////////////
///
/// Big integer representation as a vector.
///

template <typename rep_t>
class bigint_t {
    using ubig_t = typename uprod_t<rep_t>::type;
    size_t argc;
    rep_t *argv;
public:
    static const size_t base = FftPrime<rep_t>::base;

    bigint_t(size_t nwords) {
        argc = nwords;
        argv = new rep_t[argc];
        for (unsigned i = 0; i < nwords; i++) argv[i] = 0;
    }
    bigint_t(std::initializer_list<rep_t> v) {
        argc = v.size();
        argv = new rep_t[argc];
        size_t i = 0;
        for (rep_t e : v) argv[i++] = e;
    }

    ~bigint_t() { delete[] argv; }

    size_t size      () const   { return argc; }
    rep_t &operator[](size_t i) { return argv[i]; }

    void copyIn(bigint_t& other) {
        size_t lim = min2(argc, other.argc);
        size_t i = 0;
        for ( ; i < lim;  i++) argv[i] = other.argv[i];
        for ( ; i < argc; i++) argv[i] = 0;
    }
};


///////////////////////////////////////////////////////////////////////////////
//
// In place FFT over arbitrary field by Cooley-Tukey method.
// 
// Generalized from Charles Van Loan "Computational Frameworks 
// for the Fast Fourier Transform" SIAM 1992.
//
// The sin/cos formulation for roots of unity over C has been
// replaced with general fields with 2**n-th roots of unity.
// 
// Stephen M. Watt 2023-11-21.
//
#include <cstdio>
bool db = false;

template <typename _Vec_t, typename _Field_t>
struct InPlaceFFT {

    using Vec_t   = _Vec_t;
    using Field_t = _Field_t;

    // In-place bit reversal permutation.
    // Van Loan Algorithm 1.5.2.
    static void permute(Vec_t & x) {
        size_t   n = x.size();
        unsigned t = ceilLg(n);
        for (size_t k = 0; k < n; k++) {
            size_t j = bitReverse(k, t);
            if (j > k) { Field_t t = x[j]; x[j] = x[k]; x[k] = t; }
        }
    }

    // Check the parameters for FFT.
    static void checkParameters(size_t n, size_t t, Field_t omega) {
if (db) printf("at w\n");
        if (n != (1 << t))
            throw std::invalid_argument("Vector length not power of 2.");

if (db) printf("at x\n");
        Field_t o_nby2{omega.pow(n/2)}, one{1};
if (db) printf("at y\n");
        if (o_nby2 == one || o_nby2 * o_nby2 != one) {
            if (db)
                printf("omega = %llu, omega**%llu = %llu, omega**%llu = %llu\n",
                   (unsigned long long) +omega,
                   (unsigned long long) n/2, (unsigned long long) +o_nby2,
                   (unsigned long long) n, (unsigned long long) +(o_nby2 * o_nby2));
            throw std::invalid_argument("Omega is wrong order root of 2.");
        }
if (db) printf("at z\n");
    }

    // In-place Cooley-Tukey FFT.
    // Modified Van Loan Algorithm 1.6.1, 
    // 1. extended to field with primitive 2**n-th root of unity,
    // 2. replaced sin/cos with one multiply in j (middle) loop.
    //    Can do because we have exact multiplication.
 
    static void forwardFFT1(Vec_t &x, Field_t omega, bool check = false) {
        size_t n = x.size(), t = ceilLg(n);
if (db) printf("at 0\n");
        if (check) checkParameters(n, t, omega);
        permute(x);

if (db) printf("at 1\n");
        for (unsigned q = 1; q <= t; q++) {
            size_t  L = 1 << q, r = n/L;
if (db) printf("at 2\n");
            Field_t omegaPow{1}, omegaStep{omega.pow(r)};
if (db) printf("at 3\n");
            for (size_t j = 0; j < L/2; j++) {
                if (j > 0) omegaPow = omegaPow * omegaStep;
                for (size_t k = 0; k < r; k++) {
                    size_t  kLj  = k*L+j;
                    Field_t tau  = omegaPow * x[kLj + L/2];
                    x[kLj + L/2] = x[kLj] - tau;
                    x[kLj]       = x[kLj] + tau;
                }
            }
        }
    }

    // In-place Cooley-Tukey FFT.
    // Modified Van Loan Algorithm 1.6.2, 
    // 1. extended to field with primitive 2**n-th root of unity,
    // 2. avoided wlong work space by computing omega power in j (inner) loop
    //    and doing one multiply in each j iteration.
    //    The extra multiplication allows a stride of 1 with no extra storage.
 
    static void forwardFFT2(Vec_t &x, Field_t omega, bool check = false) {
        size_t n = x.size(), t = ceilLg(n);
        if (check) checkParameters(n, t, omega);
        permute(x);

        // START LOGGING
        printVec("After Bit-Reversal (C++)", x);
        // END LOGGING

        for (unsigned q = 1; q <= t; q++) {
            size_t  L = 1 << q, r = n/L;
            for (size_t k = 0; k < r; k++) {
                Field_t omegaPow{1}, omegaStep{omega.pow(r)};
                for (size_t j = 0; j < L/2; j++) {
                    if (j > 0) omegaPow = omegaPow * omegaStep;
                    size_t  kLj  = k*L+j;
                    Field_t tau  = omegaPow * x[kLj + L/2];
                    x[kLj + L/2] = x[kLj] - tau;
                    x[kLj]       = x[kLj] + tau;
                }
            }
            //printVec("After Stage (C++)", x);
        }
    }


    // Inverse FFT
    static void inverseFFT(void (*forward)(Vec_t&, Field_t, bool),
                           Vec_t& x, Field_t omega, bool check = false)
    {
        size_t  n    = x.size();
        Field_t nInv = Field_t(n).inv();

        forward(x, omega.inv(), check);
        for (size_t i = 0; i < n; i++) x[i] = nInv * x[i];
    }

    static void inverseFFT1(Vec_t& x, Field_t omega, bool check = false) {
        inverseFFT(forwardFFT1, x, omega, check);
    }
    static void inverseFFT2(Vec_t& x, Field_t omega, bool check = false) {
        inverseFFT(forwardFFT2, x, omega, check);
    }
};

///////////////////////////////////////////////////////////////////////////
///
/// Testing code
///

#include <iostream>
#include <cstdio>

template <typename V>
void printVec(const char *name, V& v, size_t nPerLine = 16) {
    size_t n = v.size();
    printf("%s = [", name);
    for (int i = 0; i < n; i++) {
        printf("%llu", (long long unsigned) +v[i]);
        printf("%s", i==n-1 ? "" : (i+1)%nPerLine == 0 ? ",\n      " : ", ");
    }
    printf("]\n");
}

void checkUtils() {
    for (unsigned n = 0; n < 12; n++)
        printf("n= %2u, ceilLg(n)= %u, is2pow= %d\n", n, ceilLg(n), is2pow(n));
}

template <typename uint_t>
void checkBase() {
    // Use printf so we don't have to change back from hex.
    printf("Info for integers with %d bits:\n", CHAR_BIT * sizeof(uint_t));
    printf("p    = %lld\n",   (long long) FftPrime<uint_t>::p);
    printf("g    = %lld\n",   (long long) FftPrime<uint_t>::g);
    printf("Base = %lld\n",   (long long) FftPrime<uint_t>::base);
    printf("     = %#0llx\n", (long long) FftPrime<uint_t>::base);
}

template <typename uint_t>
void checkZmod() {
    using Zp = zmod_t<uint_t, FftPrime<uint_t>::p>;
    uint_t checkLimit = Zp::modulus + 4;

    std::cout << "Checking Z mod "<< Zp::modulus << "\n";

    // This checks all of the operations on Zp:
    for (uint_t i = 0; i < checkLimit; i++) {
        Zp lol = Zp(101),           one = Zp(1);
        Zp a   = lol.pow(i).inv(),  b = Zp(i);

        if (!(a == lol.inv().pow(i)))
            std::cout << "Error: bad power "<< i << ".\n";
        if (a - b != a+(-b))
            std::cout << "Error: bad sum "  << i <<".\n";
        if (b / a != b*a.inv())
            std::cout << "Error: bad mul "  << i <<".\n";
        if (lol.inv().pow(i)*lol.pow(i) != one)
            std::cout << "Error: bad inv "  << i << ".\n";
        if (a*a != -Zp(1) && (a+a.inv()).inv() != a/(a.pow(2) + 1))
            std::cout << "Error: bad math " << i <<".\n";
    }
}

template <typename uint_t, typename VFp, typename Fp>
void checkFFT_single(void (*fft) (VFp&, Fp, bool), void (*ift) (VFp&, Fp, bool),
                     size_t n, uint_t *v1, uint_t *v2, uint_t *v3)
{
    using Prime = FftPrime<uint_t>;

if (db) printf("at A\n");
    if (Fp::modulus != Prime::p) throw std::logic_error("Moduli don't match.");
if (db) printf("at B\n");

    Fp     g   = Fp(Prime::g).pow((uint_t) 1 << (Prime::n - ceilLg(n)));
    printf("n  = %d\n",   (int) n);
    printf("p  = %llu\n", (unsigned long long) Prime::p);
    printf("g  = %llu\n", (unsigned long long) Prime::g);
    
    // START LOGGING
    printf("Roots of Unity (C++): g = %llu, omega = %llu\n",
       (unsigned long long) Prime::g,
       (unsigned long long) Prime::g);
    // END LOGGING

    VFp f1(n), f2(n), f3(n);
    for (int i = 0; i < n; i++) { f1[i] = Fp(v1[i]); f2[i] = Fp(v2[i]); }
    printVec("f1", f1);
    printVec("f2", f2);

if (db) printf("at C\n");
    fft(f1, g, true);
    fft(f2, g, true);
if (db) printf("at D\n");
    for (int i = 0;i < n; i++) f3[i] = f1[i] * f2[i];
    printVec("f1", f1, 8);
    printVec("f2", f2, 8);
    printVec("f3", f3, 8);

    ift(f3, g, true);
    printVec("f3", f3, 8);
    int errs = 0;
    for (int i = 0; i < n; i++) if (+f3[i] != v3[i]) errs++;
    if (!errs) printf("OK!\n"); else printf("Not OK: %d errors\n", errs);
}

template <typename uint_t>
struct fftData {
    size_t n;
    uint_t *in1, *in2, *out;
}; 

// Base 10**2
const size_t n_16 = 16;
uint16_t in1_16[] = { 38,  0, 44, 87,  6, 45, 22, 93, 0, 0, 0, 0, 0, 0, 0, 0 };
uint16_t in2_16[] = { 80, 18, 62, 90, 17, 96, 27, 97, 0, 0, 0, 0, 0, 0, 0, 0 };
uint16_t out_16[] = { 3040,    684,  5876, 11172,  5420, 16710, 12546, 20555,
                      16730, 15704, 21665,  5490, 13887,  4645,  9021,0 };

// Base 2**15
const size_t n_32 = 16;
uint32_t in1_32[] =
    { 11400, 28374, 23152, 9576, 29511, 20787, 13067, 14015, 0, 0, 0, 0, 0, 0, 0, 0 };
uint32_t in2_32[] =
    { 30268, 20788, 8033, 15446, 26275, 11619, 2494, 7016, 0, 0, 0, 0, 0, 0, 0, 0 };
uint32_t out_32[] =
    { 345055200, 1095807432, 1382179648, 1175142886, 2016084656, 2555168834,
      2179032777, 1990011337, 1860865174, 1389799087, 942120918, 778961552,
      341270975, 126631482, 98329240, 0 };

// Base 10**9
const size_t n_64 = 64;
uint64_t in1_64[] =
    { 33243586, 638827078, 767661659, 778933286, 790244973, 910208076, 425757125,
      478004096, 153380495, 205851834, 668901196, 15731080, 899763115, 551605421,
      181279081, 600279047, 711828654, 483031418, 737709105, 20544909, 609397212,
      201989947, 215952988, 206613081, 471852626, 889775274, 992608567, 947438771,
      969970961, 676943009, 934992634, 922939225, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
uint64_t in2_64[] =
    { 194132110, 219972873, 66644114, 902841100, 565039275, 540721923, 810650854,
      702680360, 147944788, 859947137, 59055854, 288190067, 537655879, 836782561,
      308822170, 315498953, 417177801, 640439652, 198304612, 525827778, 115633328,
      285831984, 136721026, 203065689, 884961191, 222965182, 735241234, 746745227,
      667772468, 739110962, 610860398, 965331182, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
uint64_t out_64[] =
    { 6453647494146460, 131329535698517158, 291767894660778388, 392668443347293259,
      971459521481104784, 1474458811520325621, 1844928110064910283, 2357021332184901128,
      2928892267161886295, 2725517850003984528, 3202505799926570519, 2918543444592941968,
      2772488376791744089, 3248633108357294538, 3254615389814072180, 3638020871734883400,
      55160505208503622, 3969469665294621400, 439789777768675993, 916737048670338429,
      157193402339279849, 1030499289809835368, 534708807109284987, 462608833776141716,
      518270737313306417, 990302136704222252, 862673986833243374, 1706781055673683080,
      2148213235654123180, 4027029548560043607, 3715706394243238489, 966330325631268533,
      724857759400778139, 1014165568394318451, 978244158856038395, 3518954508900415555,
      3481727912868647859, 2905676401026905092, 1913454655595000205, 2281030150295966751,
      2048468707271352286, 1955651308030723278, 1936345891479581000, 2116568874488615349,
      1964776204460631657, 594938508019154838, 665031798826217600, 435270820221219547,
      3944115800695200119, 3877068415832542765, 3375534600145876311, 3739051895812367546,
      3787681810231019302, 3846806706428246918, 215267241912496193, 433277273552403593,
      32647322247915044, 4082693161306839314, 3321007834415954245, 2657237599459774692,
      1906778666014199420, 1466364566853824938, 890942012983413950, 0 };

template <> struct fftData<uint16_t> {
    size_t   n = n_16;
    uint16_t *in1 = in1_16, *in2 = in2_16, *out = out_16;
};
template <> struct fftData<uint32_t> {
    size_t   n = n_32;
    uint32_t *in1 = in1_32, *in2 = in2_32, *out = out_32;
};
template <> struct fftData<uint64_t> {
    size_t   n = n_64;
    uint64_t *in1 = in1_64, *in2 = in2_64, *out = out_64;
};

template <typename uint_t>
void checkFFT() {
    using Fp  = zmod_t<uint_t, FftPrime<uint_t>::p>;
    using VFp = bigint_t<Fp>;
    using FFT = InPlaceFFT<VFp, Fp>;

    fftData<uint_t> d;

    printf("Testing FFT1\n");
    checkFFT_single(FFT::forwardFFT1, FFT::inverseFFT1, d.n, d.in1, d.in2, d.out);

    printf("Testing FFT2\n");
    checkFFT_single(FFT::forwardFFT2, FFT::inverseFFT2, d.n, d.in1, d.in2, d.out);
}


int main(int argc, char **argv) {
    checkUtils();

    checkBase<uint16_t>();
    checkBase<uint32_t>();
    checkBase<uint64_t>();

    checkZmod<uint16_t>();

    //checkFFT <uint16_t>();
    //checkFFT <uint32_t>();
db = true;
    checkFFT <uint64_t>();

    return EXIT_SUCCESS;
}
