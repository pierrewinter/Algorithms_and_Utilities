import numpy as np
from itertools import combinations

__author__ = 'Pierre Winter'

class sum_of_multiples(object):
    def __init__(self, n, divisor):
        """
        This module calculates the sum of all the multiples of a divisor below an integer n.
        Divisor can also be a list of divisors and in this case, we avoid overcounting by subtracting multiply counted values.
        One function is a naive counting approach, the other uses an improved Gauss sum.

        :param n: A positive integer which is the upper bound of what multiples are passed into the sum
        :param divisor: A positive integer (or list of integers) which define the amount that a multiple should be divisible by
        """

        self.n = int(n)
        self.divisor = divisor
        assert self.n > 0, 'You need to provide a positive integer for n'
        for div in divisor: # assert divisor array does not contain integer multiples of itself
            num = div/np.array(divisor)
            for fl in num:
                if fl.is_integer():
                    assert fl == 1.0, 'Your divisor contains an integer multiple of itself'

        lcm_prods = []  # need to subtract sums that are due to all least common multiple (LCM) combinations of divisors
        for i in range(2,len(self.divisor)):
            lcm_prods.append([np.prod(x) for x in combinations(self.divisor, i) if np.prod(x) < self.n])  # find LCM for each 2-fold, 3-fold, 4-fold, etc. product combination
        lcm_list = [int(item) for sublist in lcm_prods for item in sublist]  # flatten list into one list of integers
        if int(np.prod(self.divisor)) < self.n:  # add product that is due to LCM of all divisors
            lcm_list.append(int(np.prod(self.divisor)))
        self.lcm_list = lcm_list

    def sum_of_multiples_naive(self):
        """
        Find the naive sum of all the multiples of divisor below an integer n.
        """
    
        if not isinstance(self.divisor, list):  # if divisor is a single item
            sum_list = []
            for i in range(1,self.n):
                if i % self.divisor == 0:
                    sum_list.append(i)
                    continue
            summ = sum(sum_list)
            return summ
    
        else:  # if divisor is a list of divisors
            sum_list1 = []
            for i in range(1,self.n):
                for div in self.divisor:
                    if i % div == 0:
                        sum_list1.append(i)
                        continue
            sum1 = sum(sum_list1)
    
            sum_list2 = []
            for i in range(1,self.n):
                for div in self.lcm_list: # subtract the LCM sums
                    if i % div == 0:
                        sum_list2.append(i)
                        continue
            sum2 = sum(sum_list2)
            summ = sum1 - sum2
            return summ

    def sum_of_multiples_improved(self):
        """
        Find the sum of all the multiples of divisor below an integer n using a Gauss sum.
        """

        if not isinstance(self.divisor, list):  # if divisor is a single item
            p = (self.n-1) // self.divisor
            summ = self.divisor*(p*(p+1)) // 2
            return summ
    
        else:  # if divisor is a list of divisors
            summ = 0
            for div in self.divisor:
                p = (self.n-1) // div
                x = div*(p*(p+1)) // 2
                summ += x
    
            for lcm in self.lcm_list: # subtract the LCM sums
                p = (self.n-1) // lcm
                x = lcm*(p*(p+1)) // 2
                summ -= x
            return summ


if __name__ == '__main__':
    """
    This is an example of how to use this module.
    Both functions should return the same sum, although the naive function is unfeasible for large max_n.
    Here we calculate the sums for all the multiples of a list of divisors [4,7,11,15] below 1000000.
    
    For example, if we set max_n to 50, then we would want to sum up the
    following integers which are multiple of 4, 7, 11, or 15:
    
    SUM of [4,  7,  8, 11, 12, 14, 15, 16, 20, 21, 22, 24, 28, 30, 32, 33, 35, 36, 40, 42, 44, 45, 48, 49] = 636
    """

    max_n = 1000000
    divisors = [4,7,11,15]
    sum_class = sum_of_multiples(n=max_n, divisor=divisors)

    sum_improved = sum_class.sum_of_multiples_improved()
    print('Improved Sum\n{}'.format(sum_improved))

    sum_naive = sum_class.sum_of_multiples_naive()
    print('Naive Sum\n{}'.format(sum_naive))
