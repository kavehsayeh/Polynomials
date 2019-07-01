def gcf(a, b):
	"""
	Takes two positive integers and returns the largest factor they have in 
	common. By convention, gcf(a, 0) == a.
	Uses the fact that if a > b, gcf(a, b) == gcf(b, a % b).
	"""
	if b == a or b == 0:
		return a
	elif a == 0:
		return b
	elif a > b:
		return gcf(b, a % b)
	elif b > a:
		return gcf(a, b % a)


def lcm(a, b):
	"""
	Takes two positive integers and returns the smallest number that is a 
	multiple of both. 
	"""
	return (a*b) / gcf(a, b)

	

class Fraction:
	def __init__(self, num, denom=1):
		if denom == 0:
			raise ZeroDivisionError
		elif isinstance(num, int) and isinstance(denom, int):
			self.num = num
			self.denom = denom
			self.__clean()
		else:
			raise TypeError

	def __clean(self):
		"""Cleanup function. Simplifies the fraction and moves the negative
		sign to the numerator."""
		g = gcf(self.num, self.denom)
		if g != 1:
			self.num /= g
			self.denom /= g
		if self.denom < 0:
			self.num = -self.num
			self.denom = -self.denom

	def __str__(a):
		if denom == 1:
			return str(a.num)
		else:
			return "({}/{})".format(a.num, a.denom)

	def __float__(a):
		return a.num / a.denom

	def __int__(a):
		return a.num // a.denom

	def __neg__(a):
		return Fraction(-a.num, a.denom)

	def __add__(a, b):
		return Fraction(a.num*b.denom + a.denom*b.num, a.denom*b.denom)

	def __sub__(a, b):
		return a + (-b)

	def __mul__(a, b):
		return Fraction(a.num * b.num, a.denom * b.denom)

	def __truediv__(a, b):
		return Fraction(a.num * b.denom, a.denom * b.num)




class Monomial:
	def __init__(self, coeff=1, degree=0):
		self.degree = degree
		self.coeff = coeff

	def __str__(a):
		if a.coeff == 1:
			return "x^{}".format(a.degree)
		elif a.coeff == -1:
			return "-x^{}".format(a.degree)
		elif a.coeff == 0:
			return "0"
		else:
			return "{}x^{}".format(a.coeff, a.degree)

	def __neg__(a):
		return Monomial(-a.coeff, a.degree)

	def __add__(a, b):
		# Monomial + Monomial case
		if isinstance(b, Monomial): 
			# if the Monomials are of the same degree, add coefficients
			if a.degree == b.degree:
				return Monomial(a.coeff + b.coeff, a.degree)
			# otherwise, just return a Polynomial made of the Monomials
			else:
				return Polynomial([a, b])
		# Monomial + Polynomial case
		elif isinstance(b, Polynomial):
			# if there is already a term of the same degree, add their coefficients
			if a.degree in [a.degree for a in b.terms]:
				newterm = a + b[a.degree] # recursively call Monomial + Monomial case
				return Polynomial(list(filter(lambda x: x.degree != a.degree, 
					b.terms)) + [newterm])
			# otherwise, just add the Monomial as another term
			else:
				return Polynomial(b.terms + [a])

	def __sub__(a, b):
		return a + (-b)

	def __mul__(a, b):
		return Monomial(a.coeff * b.coeff, a.degree + b.degree)

	def eval(self, x):
		return self.coeff * x**self.degree


class Polynomial:
	def __init__(self, terms):
		# should usually be initializing Polynomials as a sum of Monomials
		# using overloaded addition
		self.terms = terms # list of Monomials
		self.degree = max(a.degree for a in terms)
		self.__clean()

	def __clean(self):
		""" 
		A cleanup helper function.
		Ensures that the terms are ordered by degree, and that
		terms with coefficient 0 are eliminated.
		TODO: consolidate terms of same degree
		"""
		self.terms = list(filter(lambda x: x.coeff != 0, self.terms))
		self.terms.sort(key = lambda x: x.degree, reverse = True)

	# most overloaded operators work by recursing to the Monomial class in 
	# some way
	def __str__(a):
		return ' + '.join(str(term) for term in a.terms)

	def __getitem__(a, degree):
		if degree in [i.degree for i in a.terms]:
			return next(x for x in a.terms if x.degree == degree)
		else:
			return Monomial(0, degree)

	def __neg__(a):
		return Polynomial([-term for term in a.terms])

	def __add__(a, b):
		# Polynomial + Monomial case
		if isinstance(b, Monomial):
			# recursively call Monomial + Polynomial case
			return b + a
		# Polynomial + Polynomial case
		elif isinstance(b, Polynomial):
			out = Polynomial([i for i in a.terms]) # necessary to create a deep copy of a
			# recursively call Monomial + Polynomial case
			for term in b.terms:
				out = term + out
			return out

	def __sub__(a, b):
		return a + (-b)

	def __mul__(a, b):
		return Polynomial([i*j for i in a.terms for j in b.terms])

	def eval(self, x):
		return sum(a.eval(x) for a in self.terms)

	def find_real_roots(self):
		if self.degree == 0: # constant
			# either the Polynomial has no roots or infinitely many roots
			return []
		elif self.degree == 1: # linear
			return [-(self[0].coeff) / self[1].coeff]
		elif self.degree == 2: # quadratic
			delta = (self[1].coeff)**2 - 4 * (self[2].coeff) * (self[0].coeff)
			if delta < 0:
				# no real roots
				return []
			elif delta == 0:
				# two identical real roots
				return 2*[-(self[1].coeff) / self[2].coeff]
			else:
				# two different real roots
				return [(-(self[1].coeff) + delta**(1/2)) / self[2].coeff, 
						(-(self[1].coeff) - delta**(1/2)) / self[2].coeff].sort()
		# there are explicit formulae for the roots of a cubic or quartc, but 
		# they're horrendously long. We'll just use the rational roots theorem 
		# for degree > 2, even though the explicit formulae are more efficient
		# for large coefficients.
		else:
			return []
			# TODO: implement rational roots test, using Fraction class