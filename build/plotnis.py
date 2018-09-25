import numpy as np
import csv
import matplotlib.pyplot as plt

nis_value = 7.815
numbers = []

with open('NIS.dat') as file:
	reader = csv.reader(file)
	for line in reader:
		numbers.append(line)
 
#sz = len(numbers)
#del(numbers[sz-10:sz-1])
sz = len(numbers)
print(sz)
nis = np.zeros(sz)
for i in range(0,sz):
    nis[i] = nis_value
x_axis = np.linspace(0,len(numbers),len(numbers))

count = 0; 
for i in numbers:
    num = i[0]
    if( float(num) > nis_value):
        count = count +1
       
#print("Percentage of values > 7.8: ", count/sz, "%", '\n')

plt.figure()
plt.plot(x_axis,numbers, x_axis, nis)
plt.ylabel('NIS Value')
plt.xlabel('Sample Number')
plt.title('Normalized Innovation Squared')
plt.legend(('NIS Values','Chi Squared Value (7.815)'))
plt.show()


