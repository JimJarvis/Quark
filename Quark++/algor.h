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
 */
void simon_period(uint64_t period);

void shor_factorize(int prime1, int prime2);

void grover_search(oracle_function);

#endif // algor_h__