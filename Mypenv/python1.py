import numpy as np 
num=2.0
# We have to find zero of "y=x^2-num" our initial guess is guess
# NewtonRapshon method is used
def sqrt(num):
	guess=1
	for i in range(100):
		guess=0.5*(guess+num/guess) # =guess-(guess^2-num)/(2*guess)
	print(f"The estimated root of {num} is {guess:.5f}.")
	ac_r=np.sqrt(num)
	print(f"The accuracy of prediction is {abs(ac_r-guess)}.")
