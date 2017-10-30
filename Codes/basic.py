#!/usr/bin/env python

"""An implementation of the Quine McCluskey algorithm.

This implementation of the Quine McCluskey algorithm has no inherent limits
(other than the calculation time) on the size of the inputs.

Also, in the limited tests of the author of this module, this implementation is
considerably faster than other public Python implementations for non-trivial
inputs.

Another unique feature of this implementation is the possibility to use the XOR
and XNOR operators, in addition to the normal AND operator, to minimise the
terms. This slows down the algorithm, but in some cases it can be a big win in
terms of complexity of the output.
"""

from __future__ import print_function
import math


class QuineMcCluskey:
	"""The Quine McCluskey class.

	The QuineMcCluskey class minimises boolean functions using the Quine
	McCluskey algorithm.

	If the class was instantiiated with the use_xor set to True, then the
	resulting boolean function may contain XOR and XNOR operators.
	"""
	__version__ = "0.2"



	def __init__(self):
	    """The class constructor.

	    Kwargs:
	        use_xor (bool): if True, try to use XOR and XNOR operations to give
	        a more compact return.
	    """
	    self.n_bits = 0         # number of bits (i.e. self.n_bits == len(ones[i]) for every i).


	def __num2str(self, i):
	    """
	    Convert an integer to its bit-representation in a string.

	    Args:
	        i (int): the number to convert.

	    Returns:
	        The binary string representation of the parameter i.
	    """
	    x = ['1' if i & (1 << k) else '0' for k in range(self.n_bits - 1, -1, -1)]
	    return "".join(x)



	def simplify(self, ones, dc = []):
	    """Simplify a list of terms.

	    Args:
	        ones (list of int): list of integers that describe when the output
	        function is '1', e.g. [1, 2, 6, 8, 15].

	    Kwargs:
	        dc (list of int): list of numbers for which we don't care if they
	        have one or zero in the output.

	    Returns:
	        see: simplify_los.

	    Example:
	        ones = [2, 6, 10, 14]
	        dc = []

	        This will produce the ouput: ['--10']
	        This means x = b1 & ~b0, (bit1 AND NOT bit0)

	    Example:
	        ones = [1, 2, 5, 6, 9, 10, 13, 14]
	        dc = []

	        This will produce the ouput: ['--^^'].
	        In other words, x = b1 ^ b0, (bit1 XOR bit0).
	    """
	    terms = ones + dc
	    if len(terms) == 0:
	        return None

	    # Calculate the number of bits to use
	    # Needed internally by __num2str()
	    self.n_bits = int(math.ceil(math.log(max(terms) + 1, 2)))

	    # Generate the sets of ones and dontcares
	    ones = set(self.__num2str(i) for i in ones)
	    dc = set(self.__num2str(i) for i in dc)

	    return self.simplify_los(ones, dc)



	def simplify_los(self, ones, dc = []):
	    """The simplification algorithm for a list of string-encoded inputs.

	    Args:
	        ones (list of str): list of strings that describe when the output
	        function is '1', e.g. ['0001', '0010', '0110', '1000', '1111'].

	    Kwargs:
	        dc: (list of str)set of strings that define the don't care
	        combinations.

	    Returns:
	        Returns a set of strings which represent the reduced minterms.  The
	        length of the strings is equal to the number of bits in the input.
	        Character 0 of the output string stands for the most significant
	        bit, Character n - 1 (n is the number of bits) stands for the least
	        significant bit.

	        The following characters are allowed in the return string:
	          '-' don't care: this bit can be either zero or one.
	          '1' the bit must be one.
	          '0' the bit must be zero.
	          '^' all bits with the caret are XOR-ed together.
	          '~' all bits with the tilde are XNOR-ed together.

	    Example:
	        ones = ['0010', '0110', '1010', '1110']
	        dc = []

	        This will produce the ouput: ['--10'].
	        In other words, x = b1 & ~b0, (bit1 AND NOT bit0).

	    Example:
	        ones = ['0001', '0010', '0101', '0110', '1001', '1010' '1101', '1110']
	        dc = []

	        This will produce the ouput: ['--^^'].
	        In other words, x = b1 ^ b0, (bit1 XOR bit0).
	    """
	    self.profile_cmp = 0    # number of comparisons (for profiling)

	    terms = ones | dc
	    if len(terms) == 0:
	        return None

	    # Calculate the number of bits to use
	    self.n_bits = max(len(i) for i in terms)
	    if self.n_bits != min(len(i) for i in terms):
	        return None

	    # First step of Quine-McCluskey method : prime implicants
	    prime_implicants = self.__get_prime_implicants(terms)
	    # Remove essential terms.
	    essential_implicants = self.__get_essential_implicants(prime_implicants)
	    # Second step of Quine McCluskey method : prime implicant chart.
	    prime_implicant_chart = self.__get_prime_implicant_chart(essential_implicants,ones)
	    # Insert here Petrick's Method.
	    final_terms = self.__petricks_method(prime_implicant_chart,ones)
	    #print(final_terms)
	    return essential_implicants

	def __get_prime_implicant_chart(self,terms,ones):
		pi_chart = dict()
		for t in terms:
			perms = set(p for p in self.permutations(t))
			values = set([v for v in perms if v in ones])
			pi_chart[t] = values
		return pi_chart

	def __petricks_method(self, chart, ones):
	    """
	    Use the prime implicants to find the essential prime implicants of the
	    function, as well as other prime implicants that are necessary to cover
	    the function. 
	    chart: the prime implicants chart
	    ones: a list of indices for the minterms for which we want the function to
	    evaluate to 1.
	    """
	    covers = []
	    if len(chart) > 0:
	      covers = [set([i]) for i in chart[0]]
	    for i in xrange(1,len(chart)):
	      new_covers = []
	      for cover in covers:
	        for prime_index in chart[i]:
	          x = set(cover)
	          x.add(prime_index)
	          append = True
	          for j in xrange(len(new_covers)-1,-1,-1):
	            if x <= new_covers[j]:
	              del new_covers[j]
	            elif x > new_covers[j]:
	              append = False
	          if append:
	            new_covers.append(x)
	      covers = new_covers
	    for cover in covers:
	      primes_in_cover = [primes[prime_index] for prime_index in cover]
	    return primes_in_cover

	def __get_prime_implicants(self, terms):
	    """Simplify the set 'terms'.

	    Args:
	        terms (set of str): set of strings representing the minterms of
	        ones and dontcares.

	    Returns:
	        A list of prime implicants. These are the minterms that cannot be
	        reduced with step 1 of the Quine McCluskey method.

	    This is the very first step in the Quine McCluskey algorithm. This
	    generates all prime implicants, whether they are redundant or not.
	    """

	    # Sort and remove duplicates.
	    n_groups = self.n_bits + 1
	    marked = set()

	    # Group terms into the list groups.
	    # groups is a list of length n_groups.
	    # Each element of groups is a set of terms with the same number
	    # of ones.  In other words, each term contained in the set
	    # groups[i] contains exactly i ones.
	    groups = [set() for i in range(n_groups)]
	    for t in terms:
	        n_bits = t.count('1')
	        groups[n_bits].add(t)

	    done = False
	    while not done:
	        # Group terms into groups.
	        # groups is a list of length n_groups.
	        # Each element of groups is a set of terms with the same
	        # number of ones.  In other words, each term contained in the
	        # set groups[i] contains exactly i ones.
	        groups = dict()
	        for t in terms:
	            n_ones = t.count('1')
	            key = n_ones
	            if key not in groups:
	                groups[key] = set()
	            groups[key].add(t)

	        terms = set()           # The set of new created terms
	        used = set()            # The set of used terms

	        # Find prime implicants
	        for key in groups:
	            key_next = key+1
	            if key_next in groups:
	                group_next = groups[key_next]
	                for t1 in groups[key]:
	                    # Optimisation:
	                    # The Quine-McCluskey algorithm compares t1 with
	                    # each element of the next group. (Normal approach)
	                    # But in reality it is faster to construct all
	                    # possible permutations of t1 by adding a '1' in
	                    # opportune positions and check if this new term is
	                    # contained in the set groups[key_next].
	                    for i, c1 in enumerate(t1):
	                        if c1 == '0':
	                            self.profile_cmp += 1
	                            t2 = t1[:i] + '1' + t1[i+1:]
	                            if t2 in group_next:
	                                t12 = t1[:i] + '-' + t1[i+1:]
	                                used.add(t1)
	                                used.add(t2)
	                                terms.add(t12)

	        # Add the unused terms to the list of marked terms
	        for g in list(groups.values()):
	            marked |= g - used

	        if len(used) == 0:
	            done = True

	    # Prepare the list of prime implicants
	    pi = marked
	    for g in list(groups.values()):
	        pi |= g
	    return pi



	def __get_essential_implicants(self, terms):
	    """Simplify the set 'terms'.

	    Args:
	        terms (set of str): set of strings representing the minterms of
	        ones and dontcares.

	    Returns:
	        A list of prime implicants. These are the minterms that cannot be
	        reduced with step 1 of the Quine McCluskey method.

	    This function is usually called after __get_prime_implicants and its
	    objective is to remove non-essential minterms.

	    In reality this function omits all terms that can be covered by at
	    least one other term in the list.
	    """
	    # Create all permutations for each term in terms.
	    perms = {}
	    for t in terms:
	    	perms[t] = set(p for p in self.permutations(t))
	    	#print(perms[t])
	    # Now group the remaining terms and see if any term can be covered
	    # by a combination of terms.
	    ei_range = set()
	    ei = set()
	    groups = dict()
	    for t in terms:
	        n = self.__get_term_rank(t, len(perms[t]))
	        if n not in groups:
	            groups[n] = set()
	        groups[n].add(t)
	    for t in sorted(list(groups.keys()), reverse=True):
	        for g in groups[t]:
	            if not perms[g] <= ei_range:
	                ei.add(g)
	                ei_range |= perms[g]
	    return ei



	def __get_term_rank(self, term, term_range):
	    """Calculate the "rank" of a term.

	    Args:
	        term (str): one single term in string format.

	        term_range (int): the rank of the class of term.

	    Returns:
	        The "rank" of the term.
	    
	    The rank of a term is a positive number or zero.  If a term has all
	    bits fixed '0's then its "rank" is 0. The more 'dontcares' and xor or
	    xnor it contains, the higher its rank.

	    A dontcare weights more than a xor, a xor weights more than a xnor, a
	    xnor weights more than 1 and a 1 weights more than a 0.

	    This means, the higher rank of a term, the more desireable it is to
	    include this term in the final result.
	    """
	    n = 0
	    for t in term:
	        if t == "-":
	            n += 8
	        elif t == "1":
	            n += 1
	    return 4*term_range + n



	def permutations(self, value = ''):
	    """Iterator to generate all possible values out of a string.

	    Args:
	        value (str): A string containing any of the above characters.

	    Returns:
	        The output strings contain only '0' and '1'.

	    Example:
	        from qm import QuineMcCluskey
	        qm = QuineMcCluskey()
	        for i in qm.permutations('1--^^'):
	            print(i)

	    The operation performed by this generator function can be seen as the
	    inverse of binary minimisation methonds such as Karnaugh maps, Quine
	    McCluskey or Espresso.  It takes as input a minterm and generates all
	    possible maxterms from it.  Inputs and outputs are strings.

	    Possible input characters:
	        '0': the bit at this position will always be zero.
	        '1': the bit at this position will always be one.
	        '-': don't care: this bit can be zero or one.

	    Algorithm description:
	        This lovely piece of spaghetti code generates all possibe
	        permutations of a given string describing logic operations.
	        This could be achieved by recursively running through all
	        possibilities, but a more linear approach has been preferred.
	        The basic idea of this algorithm is to consider all bit
	        positions from 0 upwards (direction = +1) until the last bit
	        position. When the last bit position has been reached, then the
	        generated string is yielded.  At this point the algorithm works
	        its way backward (direction = -1) until it finds an operator
	        like '-', '^' or '~'.  The bit at this position is then flipped
	        (generally from '0' to '1') and the direction flag again
	        inverted. This way the bit position pointer (i) runs forth and
	        back several times until all possible permutations have been
	        generated.
	        When the position pointer reaches position -1, all possible
	        combinations have been visited.
	    """
	    n_bits = len(value)
	    res = ['0' for i in range(n_bits)]
	    i = 0
	    direction = +1
	    while i >= 0:
	        # binary constant
	        if value[i] == '0' or value[i] == '1':
	            res[i] = value[i]
	        # dontcare operator
	        elif value[i] == '-':
	            if direction == +1:
	                res[i] = '0'
	            elif res[i] == '0':
	                res[i] = '1'
	                direction = +1
	        # unknown input
	        else:
	            res[i] = '#'

	        i = i + direction
	        if i == n_bits:
	            direction = -1
	            i = n_bits - 1
	            yield "".join(res)
