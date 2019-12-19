import numpy as np

__author__ = 'Pierre Winter'

def sieve_of_eratosthenes(n, ret_list=False):
    """
    Count the number of primes below an integer n using a Sieve of Eratosthenes algorithm.
    If you would like a list of the prime numbers as well, then set ret_list=True.
    """

    n = int(n)
    assert n >= 2

    prime = [True] * (n+1)
    prime[0] = False
    prime[1] = False

    p = 2
    while p**2 <= n:
        if prime[p]:
            for i in range(p**2, n+1, p):
                prime[i] = False
        p += 1
    n_primes = np.count_nonzero(prime)

    if not ret_list:
        return n_primes
    else:
        list_primes = []
        for i in range(2, n):
            if prime[i]:
                list_primes.append(i)
        return n_primes, list_primes

if __name__ == '__main__':
    """
    This is an example of how to use this counting algorithm.
    Here we count the amount of prime numbers below 100000 and we also get a list of those prime numbers.
    """

    max_n = 100000
    count_primes, list_primes = sieve_of_eratosthenes(max_n, ret_list=True)
    print('Number of primes less than or equal to {}: \n{} \n{}'.format(max_n, count_primes, list_primes))
