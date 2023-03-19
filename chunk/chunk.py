import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import filedialog
from scipy import integrate
from scipy import stats
from scipy import optimize

##########################
### Reading .csv files ###
##########################
print('Choose Channel 2')
root = tk.Tk()
root.withdraw()
file_path = filedialog.askopenfilename()
if file_path == '':
    print('Null')
    input()
f = open(file_path)

xx2 = []
yy2 = []

for x in f:
    t = x.split(',')
    xx2.append(float(t[3])*1000)
    yy2.append(float(t[4])*1000)

f.close()
print('Choose Channel 1')

root = tk.Tk()
root.withdraw()
file_path = filedialog.askopenfilename()
if file_path == '':
    print('Null')
    input()
f = open(file_path)

xx1 = []
yy1 = []

for x in f:
    t = x.split(',')
    xx1.append(float(t[3])*1000)
    yy1.append(float(t[4])*1000)

f.close()


################################################
###          Plot data as mV vs ms           ###
################################################
poly = 6
# Raw data before being truncated, in ms and mV #
# EMF response #
x2b = np.array(xx2)
y2b = np.array(yy2)
# 1.24 Ohm resistor response #
x1b = np.array(xx1)
y1b = np.array(yy1)



#################################################
###   Plot fitted emf function as mV vs ms    ###
#################################################
fig, ax = plt.subplots()
ax.set(xlabel='time (ms)', ylabel='voltage (mV)',
       title='Induced EMF')

ax.grid()
emfFits = []
for n in range(0,2500,100):
    parsx = np.polyfit(x2b[n:n+100],y2b[n:n+100],poly)
    emfFit = np.poly1d(parsx)
    emfFits.append(emfFit)
    plt.plot(x2b[n:n+100],y2b[n:n+100])
    plt.plot(x2b[n:n+100],emfFit(x2b[n:n+100]))



#####################################################
###   Plot fitted current function as mV vs ms    ###
#####################################################
fig, ax = plt.subplots()
ax.set(xlabel='time (ms)', ylabel='voltage (mV)',
       title='Supplied Current')

ax.grid()
currentFits = []
for n in range(0,2500,100):
    parsx = np.polyfit(x1b[n:n+100],y1b[n:n+100],poly)
    currentFit = np.poly1d(parsx)
    currentFits.append(currentFit)

    plt.plot(x1b[n:n+100],y1b[n:n+100])
    plt.plot(x1b[n:n+100],currentFit(x1b[n:n+100]))



#####################################################
###                   Integrate                   ###
#####################################################
temp = []

for emfFit in emfFits:
    # Factor of 100 is conversion to Gauss - cm^2
    emfInt = -np.polyint(emfFit) * 100
    # loops is loops in coil, different based on gap vs pickup.
    loops = 6
    emfInt = emfInt/loops
    temp.append(emfInt)

emfInts = []
for n in range(0,len(temp)):
    xStart = x2b[n * 100]
    emfInts.append(temp[n] - temp[n](xStart))


#####################################################
###          Want all pieces to connect           ###
#####################################################

finalEmfInts = [emfInts[0]]
for n in range(1,len(emfFits)):
    xEnd = x2b[n * 100]
    finalEmfInts.append(emfInts[n] + finalEmfInts[n-1](xEnd))



#####################################################
###      Plot integration as Gauss-cm^2 vs ms     ###
#####################################################
fig, ax = plt.subplots()
ax.set(xlabel='Time (ms)', ylabel='Magnetic Flux (Gauss-cm^2)',
       title='Integrated EMF')

ax.grid()
for n in range(0,len(finalEmfInts)):
    xStart = n * 100
    xEnd = xStart + 100

    plt.plot(x2b[xStart:xEnd], finalEmfInts[n](x2b[xStart:xEnd]))



#####################################################
###              Remove linear offset             ###
#####################################################
# First find minimums of every section
mins = []
minsN = []
for n in range(0,len(finalEmfInts)):
    xStart = n * 100
    xEnd = xStart + 100

    minn = min(finalEmfInts[n](x2b[xStart:xEnd]))
    minnN = np.where(finalEmfInts[n](x2b[xStart:xEnd]) == minn)[0] + xStart
    mins.append(minn)
    minsN.append(x2b[minnN])
    plt.scatter(x2b[minnN],minn)

# Need to find relevant minima

relMin = []
relMinN = []
for n in range(1,len(mins)-1):
    if mins[n] < mins[n-1] and mins[n] < mins[n+1]:
        relMin.append(mins[n])
        relMinN.append(float(minsN[n]))
print(str(relMinN))
linear = np.polyfit(relMinN,relMin,1)
linear = np.poly1d(linear)
plt.plot(x2b, linear(x2b))

fixedEmfInts = []
for emfInt in finalEmfInts:
    fixedEmfInts.append(emfInt - linear)



#####################################################
###      Plot integration as Gauss-cm^2 vs ms     ###
#####################################################
fig, ax = plt.subplots()
ax.set(xlabel='Time (ms)', ylabel='Magnetic Flux (Gauss-cm^2)',
       title='Integrated EMF')

ax.grid()
for n in range(0,len(finalEmfInts)):
    xStart = n * 100
    xEnd = xStart + 100

    plt.plot(x2b[xStart:xEnd], fixedEmfInts[n](x2b[xStart:xEnd]))



#####################################################
###     Plot integration as Gauss-cm^2 vs mA      ###
#####################################################
fig, ax = plt.subplots()
ax.set(xlabel='current (ms)', ylabel='Magnetic Flux (Gauss-cm^2)',
       title='Integrated EMF')

ax.grid()
for n in range(0,len(fixedEmfInts)):
    xStart = n * 100
    xEnd = xStart + 100

    plt.scatter(currentFits[n](x2b[xStart:xEnd]), fixedEmfInts[n](x2b[xStart:xEnd]))
plt.show()
