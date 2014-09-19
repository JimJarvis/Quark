#ifndef algor_h__
#define algor_h__

#include "qureg.h"

/*
 *	f(x) = (u * x) mod 2
 * nqubit = nbit + 1 output
 * return found secret_u
 */
uint64_t deutsch_josza_parity(int nbit, uint64_t secret_u, bool dense = true);

/*
 *	Find s such that  f(x) = f(x [+] s) 
 * where [+] is the modulo addition operator
 * return both processed Qureg for inspection and found period
 * non-destructive measurement
 */
std::pair<Qureg, uint64_t> simon_period(int nbit, uint64_t period, bool dense = true);

// classical gcd
int gcd(int, int);

/*
 *	Modular exponentiation
 * b^e mod m
 */
uint64_t exp_mod(uint64_t b, uint64_t e, uint64_t m);

// 
/*
 *	produce a shor's algorithm oracle
 * f(x) = b^x mod M, where M is the int to be factored
 */
oracle_function shor_oracle(int b, int M);
/*
 *	Return the factorized result
 */
std::pair<int, int> shor_factorize(int nbit, int prime1, int prime2);

void grover_search(oracle_function);

/*
 *	Teleport: demo Bell state entanglement
 */
void teleportation();

#endif // algor_h__