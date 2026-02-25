import numpy as np 
num=2.0
# We have to find zero of "y=x^2-num" our initial guess is guess
# NewtonRapshon method is used
def sqrt(num):
	guess=1
	for i in range(100):
		guess=0.5*(guess+num/guess) # =guess-(guess^2-num)/(2*guess)
	print(guess)
	ac_r=np.sqrt(num)
	assert abs(ac_r-guess)/ac_r < 1e-10, "The accuracy of prediction is {abs(ac_r-guess)}."

def main(lis):
	for ele in lis:
		sqrt(ele)

print("This function is run successfully.")
print("It has two functions, sqrt(num), main(lis)")
print("sqrt takes single number as input whereas main takes a list")
