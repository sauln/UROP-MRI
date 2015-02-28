##
#
# this script will open the  image/stackToR.txt file and load it as data
#


#read in data
mydata = read.table("C:/Users/nathaniel/Documents/Development/UROP-MRI/image/stackToR.txt", header=T)

print(mydata)

firstT = mydata[1:5,1:3]
secondT = mydata[6:10,1:3]
#now figure out how to run hotelling test on this data.
fit = hotelling.test(firstT, secondT)


print(fit)

