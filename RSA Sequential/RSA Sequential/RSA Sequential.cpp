#include <iostream>
#include<time.h>
#include <gmp.h>

using namespace std;

// Function to calculate modular exponentiation (x^y mod m)

void modularExponentiation(mpz_t result, const mpz_t base, const mpz_t exp, const mpz_t mod) {
    mpz_powm(result, base, exp, mod);
}




// Function to find the modular multiplicative inverse using Extended Euclidean algorithm
bool key(mpz_t result, const mpz_t a, const mpz_t m) {
    return mpz_invert(result, a, m) != 0; // Returns true if the inverse exists
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
// Function to calculate the greatest common divisor (GCD) using Euclidean algorithm
void gcd(mpz_t result, const mpz_t a, const mpz_t b) {
    mpz_gcd(result, a, b);
}

#define size 100
int main() {
   // mpz_t* Message = malloc(size * sizeof(mpz_t));
    mpz_t p, q, n, phi, e, d;
    mpz_inits(p, q, n, phi, e, d, nullptr);

    mpz_t Message[size], encrypted[size], decrypted[size];
    for (int i = 0; i < size; i++) {
        mpz_inits(Message[i], encrypted[i], decrypted[i], nullptr);
    }
    srand(time(NULL));
    time_t start_t, end_t;
    mpz_set_str(p, "10259", 10);  // p = 12345678901234567890
    mpz_set_str(q, "10253", 10);  // q = 98765432109876543210

    mpz_mul(n, p, q); // Modulus
    mpz_t p_minus__1, q_minus__1;
    mpz_inits(p_minus__1, q_minus__1, nullptr);
    mpz_sub_ui(p_minus__1, p, 1);
    mpz_sub_ui(q_minus__1, q, 1);
    mpz_mul(phi, p_minus__1, q_minus__1); // Euler's totient function
    // Ensure e and phi are coprime
    mpz_set_ui(e, 2);
    while (mpz_cmp(e, phi) < 0) {
        gcd(p_minus__1, e, phi);
        if (mpz_cmp_ui(p_minus__1, 1) == 0) { // e and phi are coprime
            break;
        }
        mpz_add_ui(e, e, 1);
    }
    // Calculate private exponent (d)
     // Calculate private exponent (d)
    if (!key(d, e, phi)) {
        cout << "Modular inverse does not exist!" << endl;
        return 1;
    }

     //gmp_printf("e: %Zd\n", e);
    // gmp_printf("n: %Zd\n", n);
     //gmp_printf("phi: %Zd\n", phi);
    gmp_printf("d: %Zd\n", d);
    // Message to be encrypted (as string)
     // Generate random message blocks
    int max = 4, min = 2;
    int range = max - min + 1;
    for (int i = 0; i < size; i++) {
        mpz_set_ui(Message[i], rand());//% range + min);
    }

    start_t = clock();

    // Encrypt and decrypt message blocks
    for (int i = 0; i < size; i++) {
        // Encryption
        modularExponentiation(encrypted[i], Message[i], e, n);

        // Decryption
        modularExponentiation(decrypted[i], encrypted[i], d, n);

        // Output results for each block
      // gmp_printf("Original: %Zd\n", Message[i]);
       //gmp_printf("Encrypted: %Zd\n", encrypted[i]);
       // gmp_printf("Decrypted: %Zd\n", decrypted[i]);
    }

    end_t = clock();
    double time_taken = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    cout << "Time taken: " << time_taken << " seconds" << endl;

    // Compare the arrays
    if (areArraysEqual(Message, decrypted, size)) {
        cout << "Arrays are equal." << endl;
    }
    else {
        cout << "Arrays are not equal." << endl;
    }

    // Clear memory
    mpz_clears(p, q, n, phi, e, d, p_minus__1, q_minus__1, nullptr);
    for (int i = 0; i < size; i++) {
        mpz_clears(Message[i], encrypted[i], decrypted[i], nullptr);
    }

    return 0;

}