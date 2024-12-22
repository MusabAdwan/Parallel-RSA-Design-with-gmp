#include <iostream>
#include <time.h>
#include <mpi.h>
#include <gmp.h>

using namespace std;

void modularExponentiation(mpz_t result, const mpz_t base, const mpz_t exp, const mpz_t mod) {
    mpz_powm(result, base, exp, mod);
}

bool key(mpz_t result, const mpz_t a, const mpz_t m) {
    return mpz_invert(result, a, m) != 0;
}

void broadcast_mpz_t(mpz_t value, int root, int myrank, MPI_Comm comm) {
    int size = 0;

    if (myrank == root) {
        size = (mpz_sizeinbase(value, 2) + 7) / 8;
    }

    MPI_Bcast(&size, 1, MPI_INT, root, comm);

    unsigned char* buffer = new unsigned char[size];

    if (myrank == root) {
        mpz_export(buffer, nullptr, 1, sizeof(unsigned char), 0, 0, value);
    }

    MPI_Bcast(buffer, size, MPI_BYTE, root, comm);

    if (myrank != root) {
        mpz_import(value, size, 1, sizeof(unsigned char), 0, 0, buffer);
    }

    delete[] buffer;
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
#define SIZE 1

int main(int argc, char** argv) {
    mpz_t p, q, e, iphi, ik, in, end_e, end_ik, chunk_size_e, chunk_size_ik, start_e, start_ik, temp;
    mpz_inits(p, q, e, iphi, ik, in, end_e, end_ik, chunk_size_e, chunk_size_ik, start_e, start_ik, temp, nullptr);
    mpz_set_str(p, "73", 10);
    mpz_set_str(q, "31", 10);

    mpz_t Message[SIZE], finalencryptionMpz[SIZE], finaldecryptionMpz[SIZE];
    for (int i = 0; i < SIZE; i++) {
        mpz_init(Message[i]);
        mpz_init(finalencryptionMpz[i]);
        mpz_init(finaldecryptionMpz[i]);
    }

    long long* encrypted = (long long*)malloc(SIZE * sizeof(long long));
    long long* decrypted = (long long*)malloc(SIZE * sizeof(long long));
    long long* finalencryption = (long long*)malloc(SIZE * sizeof(long long));
    long long* finaldecryption = (long long*)malloc(SIZE * sizeof(long long));

    if (!encrypted || !decrypted || !finalencryption || !finaldecryption) {
        cerr << "Memory allocation failed!" << endl;
        MPI_Finalize();
        return 1;
    }

  //  srand(time(NULL));
    for (int i = 0; i < SIZE; i++) {
        mpz_set_ui(Message[i], rand() % 4 + 2);
     // gmp_printf("Original: %Zd\n", Message[i]);
        finalencryption[i] = 1;
        finaldecryption[i] = 1;
    }

    //time_t start_t, end_t;
    int myrank, npes;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

 //   start_t = clock();
    double keygen_start = MPI_Wtime();
    double rsa_start = MPI_Wtime();
    if (myrank == 0) {
        mpz_mul(in, p, q);
        mpz_t p_minus_1, q_minus_1;
        mpz_inits(p_minus_1, q_minus_1, nullptr);
        mpz_sub_ui(p_minus_1, p, 1);
        mpz_sub_ui(q_minus_1, q, 1);
        mpz_mul(iphi, p_minus_1, q_minus_1);
        mpz_set_ui(e, 2);
        while (mpz_cmp(e, iphi) < 0) {
            mpz_gcd(temp, e, iphi);
            if (mpz_cmp_ui(temp, 1) == 0) {
                break;
            }
            mpz_add_ui(e, e, 1);
        }

        if (!key(ik, e, iphi)) {
            cout << "Modular inverse does not exist!" << endl;
            return 1;
        }
        mpz_clears(p, q, iphi, p_minus_1, q_minus_1, nullptr);
    }
    double keygen_end = MPI_Wtime();
    double broadcast_start = MPI_Wtime();

    broadcast_mpz_t(e, 0, myrank, MPI_COMM_WORLD);
    broadcast_mpz_t(ik, 0, myrank, MPI_COMM_WORLD);
    broadcast_mpz_t(in, 0, myrank, MPI_COMM_WORLD);
    double broadcast_end = MPI_Wtime();
    double startchunkTime = MPI_Wtime();

    mpz_add_ui(temp, e, npes - 1);
    mpz_fdiv_q_ui(chunk_size_e, temp, npes);
    mpz_add_ui(temp, ik, npes - 1);
    mpz_fdiv_q_ui(chunk_size_ik, temp, npes);

    mpz_mul_ui(start_e, chunk_size_e, myrank);
    mpz_mul_ui(end_e, chunk_size_e, myrank + 1);
    if (mpz_cmp(end_e, e) > 0) {
        mpz_set(end_e, e);
    }

    mpz_mul_ui(start_ik, chunk_size_ik, myrank);
    mpz_mul_ui(end_ik, chunk_size_ik, myrank + 1);
    if (mpz_cmp(end_ik, ik) > 0) {
        mpz_set(end_ik, ik);
    }

    mpz_t partialResult, i;
    mpz_inits(partialResult, i, nullptr);
    mpz_set_ui(partialResult, 1);
    double endchunkTime = MPI_Wtime();
    double computeencry_start = MPI_Wtime();
    for (int j = 0; j < SIZE; j++) {
       // mpz_set_ui(i, 0);
        mpz_set(i, start_e);
        // Reset i to 0 (or any appropriate starting value)
        for (; mpz_cmp(i, end_e) < 0; mpz_add_ui(i, i, 1)) {
            mpz_t msg;
            mpz_init_set(msg, Message[j]);
            mpz_mul(partialResult, partialResult, msg);
            mpz_mod(partialResult, partialResult, in);
            mpz_clear(msg);
        }
        encrypted[j] = mpz_get_ui(partialResult);
        mpz_set_ui(partialResult, 1);
    }
    double computeencry_end = MPI_Wtime();
    double MPI_Allreduceency_start = MPI_Wtime();

    MPI_Allreduce(encrypted, finalencryption, SIZE, MPI_LONG_LONG, MPI_PROD, MPI_COMM_WORLD);
    free(encrypted);
    
    double MPI_Allreduceency_end = MPI_Wtime();


    for (int j = 0; j < SIZE; j++) {
        mpz_init_set_ui(finalencryptionMpz[j], finalencryption[j]);
    }
    double computedecry_start = MPI_Wtime();

    mpz_t partial;
    mpz_init_set_ui(partial, 1);
    for (int j = 0; j < SIZE; j++) {
        mpz_mul(partial, partial, finalencryptionMpz[j]);
        mpz_mod(partial, partial, in);
        mpz_set(finalencryptionMpz[j], partial);
        mpz_set_ui(partial, 1);
    }

    mpz_clears(partial, nullptr);

    mpz_t partialResults, ii;
    mpz_inits(partialResults, ii, nullptr);
    mpz_set_ui(partialResults, 1);

    for (int j = 0; j < SIZE; j++) {
        mpz_t tempFinalEncryption;
        mpz_init_set(tempFinalEncryption, finalencryptionMpz[j]);
        mpz_set(ii, start_ik);
        for (; mpz_cmp(ii, end_ik) < 0; mpz_add_ui(ii, ii, 1)) {
            mpz_mul(partialResults, partialResults, tempFinalEncryption);
            mpz_mod(partialResults, partialResults, in);
        }
        decrypted[j] = mpz_get_ui(partialResults);
        mpz_set_ui(partialResults, 1);
        mpz_clear(tempFinalEncryption);

    }
    double MPI_Allreducedecry_start = MPI_Wtime();

    MPI_Allreduce(decrypted, finaldecryption, SIZE, MPI_LONG_LONG, MPI_PROD, MPI_COMM_WORLD);
    free(decrypted);
    double MPI_Allreducedecry_end= MPI_Wtime();
    double computedecry_end = MPI_Wtime();
    double  rsa_end = MPI_Wtime();
    if (myrank == 0)
    {// Compare the arrays

        cout << "keygen Time: " << keygen_end - keygen_start << endl;
        cout << "broadcast Time: " << broadcast_end - broadcast_start << endl;
        cout << "computeencry Time: " << computeencry_end - computeencry_start << endl;
        cout << "MPI_Allreduceency Time: " <<  MPI_Allreduceency_end - MPI_Allreduceency_start << endl;
        cout << "MPI_Allreducedecry Time: " << MPI_Allreducedecry_end - MPI_Allreducedecry_start << endl;
        cout << "chunk assignment Time: " << endchunkTime - startchunkTime << endl;
        cout << " computedecry Time: " << computedecry_end - computedecry_start << endl;
        cout << "rsa Time: " << rsa_end - rsa_start << endl;

    }
        if (myrank == 0) {
            mpz_t* mpzdecrypted = new mpz_t[SIZE];
            for (size_t i = 0; i < SIZE; i++) {
                mpz_init(mpzdecrypted[i]); // Initialize each mpz_t variable
                mpz_set_si(mpzdecrypted[i], finaldecryption[i]); // Set the value from the long long array
            }
            if (areArraysEqual(Message, mpzdecrypted, SIZE)) {
                cout << "Arrays are equal." << endl;
            }
            else {
                cout << "Arrays are not equal." << endl;
            }
            delete[] mpzdecrypted; // Clean up
        }

        for (int i = 0; i < SIZE; i++) {
            mpz_clear(Message[i]);
            mpz_clear(finalencryptionMpz[i]);
            mpz_clear(finaldecryptionMpz[i]);
        }
    mpz_clears(start_e, end_e, start_ik, end_ik, in, e, ik, temp, nullptr);
    MPI_Finalize();

    return 0;
}
