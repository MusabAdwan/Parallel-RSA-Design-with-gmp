#include <iostream>
#include <omp.h>
#include <time.h>
#include <gmp.h>

using namespace std;

// Helper Functions
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

// Define constants and global variables
#define SIZE 100
mpz_t p, q, p_minus_1, q_minus_1, e, temp, tempMessage, tempEncrypted,  iphi, ik, in, end_e, chunk_e, chunk_remainder_e, end_ik, localResult_k, chunk_k, chunkRemainder_k, start_ik, localResult;

int main() {
    // Initialization
    time_t start_t, end_t;
    mpz_inits(p, q, p_minus_1, q_minus_1, e, iphi, ik, in, end_ik, chunk_e, chunk_remainder_e , end_e, localResult_k, chunk_k, chunkRemainder_k, start_ik, end_ik, temp, localResult, tempMessage, tempEncrypted,  nullptr);

    // Initialize prime numbers
    mpz_set_str(p, "73", 10);
    mpz_set_str(q, "31", 10);

    // Configure threads
    omp_set_num_threads(7);

    // Allocate memory for arrays
    long long* Message = new long long[SIZE];
    // Start time measurement
    //srand(time(NULL));
    //start_t = clock();
    double start_rsa = omp_get_wtime();
    double key_start = omp_get_wtime();
    double encryptyioncritical_end, encryptyioncritical_start, decryptyioncritical_start,decryptyioncritical_end;
    // Modulus and Euler's totient
    mpz_mul(in, p, q); // Modulus
    mpz_sub_ui(p_minus_1, p, 1);
    mpz_sub_ui(q_minus_1, q, 1);
    mpz_mul(iphi, p_minus_1, q_minus_1); // Euler's totient function
    mpz_clears(p, q, nullptr);

    // Find a suitable value for 'e'
    mpz_set_ui(e, 2);
    while (mpz_cmp(e, iphi) < 0) {
        gcd(p_minus_1, e, iphi);
        if (mpz_cmp_ui(p_minus_1, 1) == 0) { // e and phi are coprime
            break;
        }
        mpz_add_ui(e, e, 1);
    }

    // Check if modular inverse exists
    if (!key(ik, e, iphi)) {
        cout << "Modular inverse does not exist!" << endl;
        return 1;
    }
    mpz_clears(p_minus_1, q_minus_1, nullptr);
    double key_end = omp_get_wtime();

    // Generate random messages
    for (int i = 0; i < SIZE; ++i) {
        Message[i] = rand();// % 4 + 2; // Generate random messages
    }
    mpz_t temps;
    mpz_inits(temps, nullptr);
    // Initialize GMP arrays
    mpz_t mpzMessage[SIZE], mpzencrypted[SIZE], mpzdycrypted[SIZE];
    for (int i = 0; i < SIZE; ++i) {
        mpz_inits(mpzMessage[i], mpzencrypted[i], mpzdycrypted[i], nullptr);
        mpz_set_si(mpzMessage[i], Message[i]);
        mpz_set_ui(mpzencrypted[i], 1); // Initialize encrypted to 1
        mpz_set_ui(mpzdycrypted[i], 1); // Initialize encrypted to 1

    }
    double encryptyion_start = omp_get_wtime();

    // Parallel encryption
#pragma omp parallel
    
    {
        mpz_t thread_chunk_e, thread_start_e, thread_end_e, thread_localResult, thread_tempMessage;
        mpz_inits(thread_chunk_e, thread_start_e, thread_end_e, thread_localResult, thread_tempMessage, nullptr);

        // Calculate chunk sizes
        mpz_tdiv_q_ui(thread_chunk_e, e, omp_get_num_threads());
        mpz_tdiv_r_ui(chunk_remainder_e, e, omp_get_num_threads());
        mpz_mul_ui(thread_start_e, thread_chunk_e, omp_get_thread_num());
        mpz_add(thread_end_e, thread_start_e, thread_chunk_e);

        // Add remainder to the last thread
        if (omp_get_thread_num() == omp_get_num_threads() - 1) {
            mpz_add(thread_end_e, thread_end_e, chunk_remainder_e);
        }
        mpz_t i;
        mpz_inits( i, nullptr);
        mpz_set_ui(thread_localResult, 1);
        for (int j = 0; j < SIZE; j++) {
            mpz_set(i, thread_start_e);
            for (; mpz_cmp(i, thread_end_e) <0; mpz_add_ui(i, i, 1)) {
                mpz_t msg;
                mpz_init_set(msg, mpzMessage[j]); // Convert Message[j] to mpz_t
                mpz_mul(thread_localResult, thread_localResult, msg); // Multiply
                mpz_mod(thread_localResult, thread_localResult, in); // Modulo
                mpz_clear(msg);
            }       
           // encryptyioncritical_start = omp_get_wtime();

#pragma omp critical
            {
                mpz_mul(mpzencrypted[j], mpzencrypted[j], thread_localResult); // Multiply
                mpz_mod(mpzencrypted[j], mpzencrypted[j], in); // Modulo

            }    
            mpz_set_ui(thread_localResult, 1);
          //  encryptyioncritical_end = omp_get_wtime();

        }  
        mpz_clears(thread_chunk_e, thread_start_e, thread_end_e, thread_localResult, thread_tempMessage, nullptr);
    }
    double encryptyion_end = omp_get_wtime();
    double decryptyion_start = omp_get_wtime();

#pragma omp parallel
    {
        // Initialize thread-specific variables
        mpz_t thread_chunk_k, thread_chunk_remainder_k, thread_start_ik, thread_end_ik;
        mpz_t thread_localResult_k, thread_tempMessage, thread_tempDecrypted, thread_tempEncrypted, thread_tempIn;
        mpz_inits(thread_chunk_k, thread_chunk_remainder_k, thread_start_ik, thread_end_ik,
            thread_localResult_k, thread_tempMessage, thread_tempDecrypted, thread_tempEncrypted,
            thread_tempIn, nullptr);

        // Compute chunk sizes for the decryption key
        mpz_tdiv_q_ui(thread_chunk_k, ik, omp_get_num_threads());
        mpz_tdiv_r_ui(thread_chunk_remainder_k, ik, omp_get_num_threads());
        mpz_mul_ui(thread_start_ik, thread_chunk_k, omp_get_thread_num());
        mpz_add(thread_end_ik, thread_start_ik, thread_chunk_k);

        // Add remainder to the last thread's chunk
        if (omp_get_thread_num() == omp_get_num_threads() - 1) {
            mpz_add(thread_end_ik, thread_end_ik, thread_chunk_remainder_k);
        }
        // Initialize local result for decryption
        mpz_set_ui(thread_localResult_k, 1);
        mpz_t ii;
        mpz_inits(ii, nullptr);
        // Process each message in the array
        for (int j = 0; j < SIZE; j++) {
            mpz_set(ii, thread_start_ik);
            for (; mpz_cmp(ii, thread_end_ik) <0; mpz_add_ui(ii, ii, 1)) {
                mpz_t msg;
                mpz_init_set(msg, mpzencrypted[j]);
              //  gmp_printf("%Zd\n", msg);
                mpz_mul(thread_localResult_k, thread_localResult_k, msg); // Multiply
                mpz_mod(thread_localResult_k, thread_localResult_k, in); // Modulo
              //  gmp_printf("%Zd\n", thread_localResult_k);
                mpz_clear(msg);
            }
            // Update decrypted array in a critical section
#pragma omp critical
           // decryptyioncritical_start = omp_get_wtime();
            { 
                mpz_mul(mpzdycrypted[j], mpzdycrypted[j], thread_localResult_k); // Multiply
                mpz_mod(mpzdycrypted[j], mpzdycrypted[j], in); // Modulo
            }        mpz_set_ui(thread_localResult_k, 1);
            //decryptyioncritical_end = omp_get_wtime();
        }

        // Clear thread-specific variables
        mpz_clears(thread_chunk_k, thread_chunk_remainder_k, thread_start_ik, thread_end_ik,
            thread_localResult_k, thread_tempMessage, thread_tempDecrypted, thread_tempEncrypted,
            thread_tempIn, nullptr);
    }        double decryptyion_end = omp_get_wtime();

    double end_rsa = omp_get_wtime();

    if (areArraysEqual(mpzMessage, mpzdycrypted, SIZE)) {
        cout << "Arrays are equal." << endl;
    }
    else {
        cout << "Arrays are not equal." << endl;
    }
   
    double rsa_time = (double)(end_rsa - start_rsa);
    double keygen_time = (double)(key_end - key_start);
  //  double encryptyioncritical_time = (double)(encryptyioncritical_end - encryptyioncritical_start);
   // double decryptyioncritical_time = (double)(decryptyioncritical_end - decryptyioncritical_start);
    double encryptyion_time = (double)(encryptyion_end - encryptyion_start);
    double decryptyion_time = (double)(decryptyion_end - decryptyion_start);
    printf("rsa_time taken is %f \n", rsa_time);
    printf("keygen_time taken is %f \n", keygen_time);
    //printf("encryptyioncritical_time taken is %f \n", encryptyioncritical_time);
   // printf("decryptyioncritical_time taken is %f \n", decryptyioncritical_time);
    printf("encryptyion_time taken is %f \n", encryptyion_time);
    printf("decryptyion_time taken is %f \n", decryptyion_time);

    delete[] Message;
    return 0;
}