#include <iostream>
#include <omp.h>
#include<time.h>
#include <gmp.h>

using namespace std;
void gcd(mpz_t result, const mpz_t a, const mpz_t b) {
    mpz_gcd(result, a, b);
}

bool key(mpz_t result, const mpz_t a, const mpz_t m) {
    return mpz_invert(result, a, m) != 0; // Returns true if the inverse exists
}

void modularExponentiation(mpz_t result, const mpz_t base, const mpz_t exp, const mpz_t mod) {
    mpz_powm(result, base, exp, mod);
}
bool areArraysEqual(mpz_t* a, mpz_t* b, size_t size) {
    for (size_t i = 0; i < size; i++) {
        if (mpz_cmp(a[i], b[i]) != 0) {
            // Not equal, log the difference
            cout << "Difference at index " << i
                << ": a[" << i << "] = " << mpz_get_str(nullptr, 10, a[i])
                << ", b[" << i << "] = " << mpz_get_str(nullptr, 10, b[i]) << endl;
            return false;
        }
    }
    return true;
}
#define size  100000
mpz_t p, q, p_minus_1, q_minus_1, e, iphi, ik, in, end_e, end_ik, chunk_size_e, chunk_size_ik, start_e, start_ik, temp;
int threadNum;

int main(int argc, char* argv[])
{
    mpz_inits(p, q, p_minus_1, q_minus_1, e, iphi, ik, in, end_e, end_ik, chunk_size_e, chunk_size_ik, start_e, start_ik, temp, nullptr);
    mpz_set_str(p, "130774788487146856696749840114493855840785577043415680291425861961210790634818946117017421658663513706550417032913951059982130428363585568661730674064355226564517031250971248140528001265857282971329432432143508030885375826509984410886663147158234378940978995342512108102094436363588805272654370327323356942391", 10);
    mpz_set_str(q, "127263931017426249906210455266268278260277993436627181679325191552231875948823999138536364186356534239716692849674009488231470623987931402867840218483947606279607538956067894856386441370838116441722098641506096736240215861974676626493526591428638929799703854843175699500166285569079126009121746394573961793587", 10);

    omp_set_num_threads(4);
  //  time_t start_t, end_t;
    long long* Message = new long long[size]; 


    //srand(time(NULL));
   // start_t = clock();
    double rsa_start = omp_get_wtime();
    double key_start = omp_get_wtime();

    mpz_mul(in, p, q); // Modulus
    mpz_sub_ui(p_minus_1, p, 1);
    mpz_sub_ui(q_minus_1, q, 1);
    mpz_mul(iphi, p_minus_1, q_minus_1); // Euler's totient function
    mpz_clears(p, q, nullptr);
    mpz_set_ui(e, 65530);
    while (mpz_cmp(e, iphi) < 0) {
        gcd(p_minus_1, e, iphi);
        if (mpz_cmp_ui(p_minus_1, 1) == 0) { // e and phi are coprime
            break;
        }
        mpz_add_ui(e, e, 1);
    }

    if (!key(ik, e, iphi)) {
        cout << "Modular inverse does not exist!" << endl;
        return 1;
    }
    double key_end = omp_get_wtime();

    mpz_clears(p_minus_1, q_minus_1, nullptr);
 
    int max = 4, min = 2;
    int range = max - min + 1;
    for (int i = 0; i < size; ++i) {
        Message[i] = rand() % 4 + 2; // Generate random messages
    }
    mpz_t mpzMessage[size], mpzencrypted[size], mpzdycrypted[size];
    for (int i = 0; i < size; ++i) {
        mpz_inits(mpzMessage[i], mpzencrypted[i], mpzdycrypted[i], nullptr);
        mpz_set_si(mpzMessage[i], Message[i]);

    }
    double compute_start = omp_get_wtime();

#pragma omp parallel for
        for (int i = 0; i < size; i++) {
             threadNum = omp_get_thread_num();
            // RSA encryption on the block
             modularExponentiation(mpzencrypted[i], mpzMessage[i], e, in);
             modularExponentiation(mpzdycrypted[i], mpzencrypted[i], ik, in);
#pragma omp critical
             {
                 long long enc = mpz_get_si(mpzencrypted[i]);
                 long long dec = mpz_get_si(mpzdycrypted[i]);
               //  printf("Thread %d: Message: %lld Encrypted: %lld Decrypted: %lld\n",threadNum, Message[i], enc, dec);
             }
        }
       // printf("Thread %d:  ", threadNum);
        double compute_end = omp_get_wtime();
        double rsa_end = omp_get_wtime();

        if (areArraysEqual(mpzMessage, mpzdycrypted, size)) {
            cout << "Arrays are equal." << endl;
        }
        else {
            cout << "Arrays are not equal." << endl;
        }
       // end_t = clock();
       // double time_taken = (double)(end - start); /// (double)(CLOCKS_PER_SEC);
        double compute_time = (double)(compute_end - compute_start);
        double key_gentime = (double)(key_end - key_start) ;
        double rsa_time = (double)(rsa_end - rsa_start) ;

        printf("compute_time taken is %f \n", compute_time);
        printf("key_gentime taken is %f \n", key_gentime);
        printf("rsa_time taken is %f \n", rsa_time);
      //  printf("time from time.h taken is %f \n", time_taken_t);
        for (int i = 0; i < size; ++i) {
            mpz_clears(mpzMessage[i], mpzencrypted[i], mpzdycrypted[i], nullptr);
        }
        mpz_clears(e, iphi, ik, in, nullptr);
        delete[] Message;
           return 0;

    }
