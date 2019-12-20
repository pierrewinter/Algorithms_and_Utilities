__author__ = 'Pierre Winter'

def largest_palindrome(n_digits):
    """
    Given a number of length n_digits, return the largest possible palindrome
    made by multiplying two numbers which each have n_digits. If no palindrome exists, return zeros.
    """

    assert n_digits > 0
    max_number = int(n_digits*'9')
    min_number = int('1' + (n_digits-1)*'0')
    max_product = max_number**2
    min_product = min_number**2

    for number in range(max_product, min_product, -1):  # loop over all possible numbers
        if str(number) == str(number)[::-1]:  # if palindrome
            for divisor in range(max_number, min_number+1, -1):  # loop over divisors
                if number % divisor != 0 or len(str(number // divisor)) != n_digits:  # if palindrome has no integer divisors or divisors are not both length n_digits
                        continue  # continue divisor for loop
                else:  # palindrome has integer divisors and divisors are both length n_digits
                    return number, divisor, number // divisor
        else:  # not palindrome, continue looping through numbers
            continue
    return 0, 0, 0


if __name__ == '__main__':
    """
    This is an example of how to use this palindrome search algorithm.
    Here we find the largest palindrome that can be made up of 4-digit integers.
    """

    n_digits = 4
    palindrome, divisor1, divisor2 = largest_palindrome(n_digits=n_digits)
    print('The largest palindrom made up of {}-digit integers is: \n{}'.format(n_digits, palindrome))
    print('The integer divisors are: \n{}, {}'.format(divisor1, divisor2))
