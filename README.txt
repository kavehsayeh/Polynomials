A basic polynomial class. I know numpy already has one of these, but if I only did stuff nobody else had ever done, I wouldn't do anything.

Also has a helper monomial class. There is probably a way to eliminate this class and treat monomials as polynomials of one term; I'll see if I can swing that.

So far has overloaded operators for addition, subtraction, negation, and multiplication, as well as a conversion to a string representation using overloaded str() and a way to evaluate a polynomial at a specific point.

TODO: Find the zeroes of polynomials, using explicit formulas for degree < 3 and the rational roots test for larger polynomials.