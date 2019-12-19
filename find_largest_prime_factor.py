import math
from count_primes import sieve_of_eratosthenes

__author__ = 'Pierre Winter'

def largest_prime_factor(n, ret_list=False):
    """
    For a given positive integer n, return the largest prime factor of n.
    If you would also like a list of all prime factors of n, then set ret_list=True.
    """

    n = int(n)
    assert n >= 2
    flist = []
    max_n = math.ceil(math.sqrt(n))
    prime_count, prime_list = sieve_of_eratosthenes(n=max_n, ret_list=True)  # Get list of all prime numbers up to sqrt(n)

    for prime in prime_list:
        if n % prime == 0:  # if the prime number is a factor, add it to the list and update n
            flist.append(prime)
            n //= prime

    if not ret_list:
        return max(flist)
    else:
        return max(flist), flist

if __name__=='__main__':
    """
    This is an example of how to use this prime factorization algorithm.
    Here we find the prime factors of 600851475143 and we also return the largest of these prime factors.
    """

    max_n = 600851475143
    max_prime_factor, prime_factors = largest_prime_factor(n=max_n, ret_list=True)
    print('The prime factors of {} are: \n{}\nof which the largest is: \n{}'.format(max_n, prime_factors, max_prime_factor))
