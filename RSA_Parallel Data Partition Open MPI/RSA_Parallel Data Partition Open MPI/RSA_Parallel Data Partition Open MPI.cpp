
#include <iostream>
#include <mpi.h>
//#include <time.h>
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

void broadcast_mpz_t(mpz_t value, int root, int myrank, MPI_Comm comm) {
    int size = 0;

    if (myrank == root) {
        size = (mpz_sizeinbase(value, 2) + 7) / 8; // Calculate size in bytes
    }

    // Broadcast size to all ranks
    MPI_Bcast(&size, 1, MPI_INT, root, comm);

    // Allocate buffer based on size
    unsigned char* buffer = new unsigned char[size];

    if (myrank == root) {
        mpz_export(buffer, nullptr, 1, sizeof(unsigned char), 0, 0, value);
    }

    // Broadcast the actual data
    MPI_Bcast(buffer, size, MPI_BYTE, root, comm);

    // Import the data back to mpz_t
    if (myrank != root) {
        mpz_import(value, size, 1, sizeof(unsigned char), 0, 0, buffer);
    }

    delete[] buffer; // Clean up
}
bool areArraysEqual(long long arr1[], long long arr2[], int size) {
    for (int i = 0; i < size; i++) {
        if (arr1[i] != arr2[i]) {
            return false;  // Arrays are not equal
        }
    }
    return true;  // Arrays are equal
}
bool areArraysEquals(mpz_t* a, mpz_t* b, size_t size) {
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
#define SIZE 100000

int main(int argc, char** argv) {
    mpz_t p, q, in, iphi, e, ik, p_minus_1, q_minus_1;
    mpz_inits(p, q, in, iphi, e, ik, p_minus_1, q_minus_1, nullptr);
    //long long* Message = new long long[SIZE];
    long long Message[SIZE];
    long long decrypted[SIZE];
    //srand(time(NULL));
  //  clock_t start_t, end_t;
    int myrank, npes;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    // start_t = clock();
    double keygen_start = MPI_Wtime();
    double rsa_start = MPI_Wtime();
    if (myrank == 0) {
        mpz_set_str(p, "130774788487146856696749840114493855840785577043415680291425861961210790634818946117017421658663513706550417032913951059982130428363585568661730674064355226564517031250971248140528001265857282971329432432143508030885375826509984410886663147158234378940978995342512108102094436363588805272654370327323356942391", 10);
        mpz_set_str(q, "127263931017426249906210455266268278260277993436627181679325191552231875948823999138536364186356534239716692849674009488231470623987931402867840218483947606279607538956067894856386441370838116441722098641506096736240215861974676626493526591428638929799703854843175699500166285569079126009121746394573961793587", 10);
        mpz_mul(in, p, q); // Modulus
        mpz_sub_ui(p_minus_1, p, 1);
        mpz_sub_ui(q_minus_1, q, 1);
        mpz_mul(iphi, p_minus_1, q_minus_1); // Euler's totient function
        mpz_clears(p, q, nullptr);
        for (int i = 0; i < SIZE; ++i) {
            Message[i] = rand() % 100 + 1; // Generate random messages
        }
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
        mpz_clears(p_minus_1, q_minus_1, nullptr);
    }
    double keygen_end = MPI_Wtime();

    // Broadcast Message array
    double broadcast_start = MPI_Wtime();
    MPI_Bcast(&Message, SIZE, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    // After broadcasting the mpz_t variables
    broadcast_mpz_t(e, 0, myrank, MPI_COMM_WORLD);
    broadcast_mpz_t(ik, 0, myrank, MPI_COMM_WORLD);
    broadcast_mpz_t(iphi, 0, myrank, MPI_COMM_WORLD);
    broadcast_mpz_t(in, 0, myrank, MPI_COMM_WORLD);
    double broadcast_end = MPI_Wtime();

    // Debug prints after broadcasting
    //long long* local_decrypted = new long long[SIZE / npes]; // Local array for each process
    mpz_t mpzMessage[SIZE], mpzencrypted[SIZE], mpzdycrypted[SIZE];
    for (int i = 0; i < SIZE; ++i) {
        mpz_inits(mpzMessage[i], mpzencrypted[i], mpzdycrypted[i], nullptr);
        mpz_set_si(mpzMessage[i], Message[i]);
    }
    double startchunkTime = MPI_Wtime();
    int chunkSize = SIZE / npes;
    int startpoint = myrank * chunkSize;
    int endpoint = (myrank == npes - 1) ? SIZE : startpoint + chunkSize;
    // Local decryption arrays
    double endchunkTime = MPI_Wtime();
    double compute_start = MPI_Wtime();
    for (int i = startpoint; i < endpoint; ++i) {
        modularExponentiation(mpzencrypted[i], mpzMessage[i], e, in);
        modularExponentiation(mpzdycrypted[i], mpzencrypted[i], ik, in);
        //  local_decrypted[i - startpoint] = mpz_get_si(mpzdycrypted[i]);

          // Output results for each block
    }
    double compute_end = MPI_Wtime();
    double  rsa_end = MPI_Wtime();
    //end_t = clock();
   // double time_taken = (double)(end_t - start_t) / CLOCKS_PER_SEC;
   
    // MPI_Gather(local_decrypted, chunkSize, MPI_LONG_LONG, decrypted, chunkSize, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
   

    if (myrank == 0)
    {// Compare the arrays

        cout << "keygen Time: " << keygen_end - keygen_start << endl;
        cout << "Broadcast Time: " << broadcast_end - broadcast_start << endl;
        cout << "chunk assignment Time: " << endchunkTime - startchunkTime << endl;
        cout << "Computation Time: " << compute_end - compute_start << endl;
        cout << "rsa Time: " << rsa_end - rsa_start << endl;
     
    }

    // Cleanup
    for (int i = 0; i < SIZE; i++) {
        mpz_clears(mpzMessage[i], mpzencrypted[i], mpzdycrypted[i], nullptr);
    }
    mpz_clears(in, iphi, e, ik, nullptr);
    MPI_Finalize();
    return 0;
}