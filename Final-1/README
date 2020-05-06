### In this Project there are 12 files:
The prediction output files are list as follow:
FPp1-1NNoutput.txt
FPp1-3NNoutput.txt
--------------------------
FPp1d-predictions1
FPp1d-predictions3
---------------------------
'find_k.png': An output image 
###  FPp1a-source.py file including all the code in this project

**How to complie  / run code**

This .py file is transfered from .ipynb, but I delete all the print() method in this code, and I think it is clear to check the code.
I did this project by 3 steps:
1. Input the data set and transfer then into the form I want. I read the file as DataFrame type at first, which is easier to find the 'Class" and build the y_train. This part of code is in line 14 to 41.
2. As the project description, I need to write my own knn-function, but I firstly wrote the function to predict every single nodes in X_test, if you want to get a list, you need to run the knn_list function I wrote. This part starts from 43 to 80.
3. Notice that line 73-80 is a way to find the best 'K', I plot a figure ( 'find_k.png' ), as the pic shows, k around 5-10 can get a higher accurancy.
4. Use the KFold Model in sklearn, then count the accurancy for output file FPp1d-predictions1 and FPp1d-predictions3, the results have been added to the last line.
