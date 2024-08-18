import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import statistics
from scipy import stats


### imported metadata from the OGLE IV Catalogue of Delta Scuti Variables into pandaS dataframes ###

#Small Megellanic Cloud Delta Scuti Variable Dataframe
rawsmcds = pd.read_csv ('smcdsdata.csv')
rawsmcds.columns = ['ID', 'mode', 'Ra', 'Decl', 'I', 'V', 'V-I', 'P1', 'P2']

#Large Megellanic Cloud Delta Scuti Variable Dataframe
rawlmcds = pd.read_csv ('lmcdsdata.csv')
rawlmcds.columns = ['ID', 'mode', 'Ra', 'Decl', 'I', 'V', 'V-I', 'P1', 'P2']


print(len(rawsmcds)) #Checking dataset size
print(len(rawlmcds)) #Checking dataset size

### data cleaning/pre-processing ###

#apparent magnitude thresholding - based on OGLE IV telescope saturation and sensitivity limits

rawsmcds = rawsmcds[rawsmcds['I'] > 13]
rawsmcds = rawsmcds[rawsmcds['I'] < 21.5]

rawlmcds = rawlmcds[rawlmcds['I'] > 13]
rawlmcds = rawlmcds[rawlmcds['I'] < 21.5]

print(len(rawsmcds)) #Checking how many stars removed from dataset
print(len(rawlmcds)) #Checking how many stars removed from dataset

#removal of Milkyway Halo Delta Scutis - Equal distance approximation cannot be used with Milkyway Delta Scutis - creating new dataframes

# coefficients of P-L relationship prior to cleansing - required to determine whether to remove star - the loops below take the true apparent magnitude of a star, and compare it with the apparent magnitude the star should have based on the Line of best fit
# if the star's apparent magnitude is lower than the line of best fit apparent magnitude by more than 1.5, the code keeps this star in the raw dataframe (rawsmcds / rawlmcds)
# if the star's apparent magnitude is greater than this thresholded value (1.5 lower than the LOBF apparent magnitude), the loop adds this star with its ID, apparent magnitude and period to a "cleansed" dataframe (cleanlmcds / cleansmcds)

#1. First step is plotting the line of best fit for the stars: Apparent Magnitude (m) vs Log(10) Period and acquiring the coefficients of the gradient and y intercept. The period used is the fundamental mode period.

# Function that calculates the Log10 of the Period of Pulsation

def logperiod(P):
    return(np.log10(P))

# Graphing the m vs logP relation and acquiring the coefficients of the relation: 

# a is the gradient, b is the y intercept for the SMC relationship

y = (((rawsmcds['I'])).tolist())
x = ((logperiod(rawsmcds['P1'])).tolist())

plt.scatter(x, y, marker=".")
a, b = np.polyfit(x, y, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('apparent I-band Mag')
LOBF = [i * a + b for i in x]
LOBF2 = [i * a + (b - 1.5) for i in x]
plt.plot(x, LOBF, color = "red")
plt.plot(x, LOBF2, color = 'blue')
plt.text(-0.9, 13, 'm = ' + format(a.round(3)) + 'logP ' + "+ " + format(b.round(3)))
plt.title('DS SMC raw m vs logP')
plt.ylim(12, 23)
plt.xlim(-1.6,-0.4)
plt.show()

# c is the gradient, d is the y intercept for the LMC relationship

w = (((rawlmcds['I'])).tolist())
q = ((logperiod(rawlmcds['P1'])).tolist())

plt.scatter(q, w, marker=".")
c, d = np.polyfit(q, w, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('apparent I-band Mag')
LOBF = [i * c + d for i in q]
LOBF2 = [i * c + (d - 1.5) for i in q]
plt.plot(q, LOBF, color = "red")
plt.plot(q, LOBF2, color = 'blue')
plt.text(-0.9, 13, 'm = ' + format(c.round(3)) + 'logP ' + "+ " + format(d.round(3)))
plt.title('DS LMC raw m vs logP')
plt.ylim(12, 23)
plt.xlim(-1.6,-0.4)
plt.show()

#2. Compare the true value of the apparent magnitude of stars against the line of best fit apparent magnitude derived from the star's log of Period. If the magnitude is within the accepted range, the star is added to a cleansed dataframe 

# cleansed SMC dataframe generation
cleansmcds = pd.DataFrame(columns = ['ID', 'mode', 'Ra', 'Decl', 'I', 'V', 'V-I', 'P1', 'P2']) #creating a new dataframe for the cleansed data
for index, row in rawsmcds.iterrows(): #looping through raw data dataframe
    if np.float64(row['I']) > (np.float64(np.log10(row['P1'])) * a + b - 1.5): # if a specifc star satisfies this condition, it is added into the new dataframe
        cleansmcds.loc[len(cleansmcds.index)] = [row['ID'], row['mode'], row['Ra'], row['Decl'], row['I'], row['V'], row['V-I'], row['P1'], row['P2']] # adds the star at the end of the dataframe (last index)
    else:
        continue

# cleansed LMC dataframe generation
cleanlmcds = pd.DataFrame(columns = ['ID', 'mode', 'Ra', 'Decl', 'I', 'V', 'V-I', 'P1', 'P2']) #creating a new dataframe for the cleansed data
for index, row in rawlmcds.iterrows(): #looping through raw data dataframe
    if np.float64(row['I']) > (np.float64(np.log10(row['P1'])) * c + d - 1.5): # if a specifc star satisfies this condition, it is added into the new dataframe
        cleanlmcds.loc[len(cleanlmcds.index)] = [row['ID'], row['mode'], row['Ra'], row['Decl'], row['I'], row['V'], row['V-I'], row['P1'], row['P2']] # adds the star at the end of the dataframe (last index)
    else:
        continue

print(len(rawsmcds) - len(cleansmcds)) #Checking how many stars were in pre-processing
print(len(rawlmcds) - len(cleanlmcds)) #Checking how many stars were in pre-processing

#3. Visualisation of the new dataset's m vs logP relationship

y = (((cleansmcds['I'])).tolist())
x = ((logperiod(cleansmcds['P1'])).tolist())

plt.scatter(x, y, marker=".")
m, n = np.polyfit(x, y, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('apparent I-band Mag')
LOBF = [i * m + n for i in x]
LOBF2 = [i * a + (b - 1.5) for i in x]
plt.plot(x, LOBF, color = "red")
plt.plot(x, LOBF2, color = 'blue')
plt.text(-0.9, 13, 'm = ' + format(m.round(3)) + 'logP ' + "+ " + format(n.round(3)))
plt.title('DS SMC cleansed m vs logP')
plt.ylim(12, 23)
plt.xlim(-1.6,-0.4)
plt.show()

w = (((cleanlmcds['I'])).tolist())
q = ((logperiod(cleanlmcds['P1'])).tolist())

plt.scatter(q, w, marker=".")
k, l = np.polyfit(q, w, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('apparent I-band Mag')
LOBF = [i * k + l for i in q]
LOBF2 = [i * c + (d - 1.5) for i in q]
plt.plot(q, LOBF, color = "red")
plt.plot(q, LOBF2, color = 'blue')
plt.text(-0.9, 13, 'M = ' + format(k.round(3)) + 'logP ' + "+ " + format(l.round(3)))
plt.title('DS LMC cleansed m vs logP')
plt.ylim(12, 23)
plt.xlim(-1.6,-0.4)
plt.show()
