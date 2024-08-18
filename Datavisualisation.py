import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import statistics
from scipy import stats
#from testing import cleanblgds

from dsdistances import FundamentalblgdsF, FirstOvertoneblgdsF, FundamentaldiskdsF, FirstOvertonediskdsF
from cephdistances import FundamentalblgcephF, FirstOvertoneblgcephF, FundamentaldiskcephF, FirstOvertonediskcephF

### Data Visualisation Functions ###

def modelvstrue(x1a, y1a, y1b, title="", xlabel="", ylabel="", colour1='#1984c5', colour2='#76c68f', marker1='.', marker2='.'):

    plt.scatter(x1a, y1a, c=colour1, marker=marker1, label="SMC Model", alpha=1)
    plt.scatter(x1a, y1b, c=colour2, marker=marker2, label="LMC Model", alpha=1)
    plt.xlabel('True Distance')
    plt.ylabel('Model Distances')
    plt.title(title)
    plt.legend()
    plt.show()


def errorfunction(x2a, y2a, y2b, title="", xlabel="", ylabel="", colour1='#1984c5', colour2='#76c68f', marker1='.', marker2='.'):
    a, b = np.polyfit(x2a, y2a, 1)
    c, d = np.polyfit(x2a, y2b, 1)
    LOBF = [i * a + b for i in x2a]
    LOBF2 = [i * c + d for i in x2a]
    plt.scatter(x2a, y2a, c=colour1, marker=marker1, label="SMC Model", alpha=1)
    plt.scatter(x2a, y2b, c=colour2, marker=marker2, label="LMC Model", alpha=1)
    plt.xlabel('True Distance')
    plt.ylabel('Error%')
    plt.title(title)
    plt.legend()
    plt.plot(x2a, LOBF, color = colour1)
    plt.plot(x2a, LOBF2, color = colour2)
    plt.show()


### Models Vs True Distance Graphs ###

### Calling upon the functions for all 4 dataframes ###

## Delta Scuti Variables ##

# Fundamental BLG #

x1a = FundamentalblgdsF['Distance']
y1a = FundamentalblgdsF['SMCmodeldist']
y1b = FundamentalblgdsF['LMCmodeldist']

modelvstrue(x1a, y1a, y1b, title="DS BLG Fundamental Modeldist/Truedist %", xlabel="True Dist", ylabel="Model Dist")

# Fundamental DISK #

x1a = FundamentaldiskdsF['Distance']
y1a = FundamentaldiskdsF['SMCmodeldist']
y1b = FundamentaldiskdsF['LMCmodeldist']

modelvstrue(x1a, y1a, y1b, title="DS DISK Fundamental Modeldist/Truedist %", xlabel="True Dist", ylabel="Model Dist")

# FirstOvertone BLG #

x1a = FirstOvertoneblgdsF['Distance']
y1a = FirstOvertoneblgdsF['SMCmodeldist']
y1b = FirstOvertoneblgdsF['LMCmodeldist']

modelvstrue(x1a, y1a, y1b, title="DS BLG FirstOvertone Modeldist/Truedist %", xlabel="True Dist", ylabel="Model Dist")

# FirstOvertone DISK #

x1a = FirstOvertonediskdsF['Distance']
y1a = FirstOvertonediskdsF['SMCmodeldist']
y1b = FirstOvertonediskdsF['LMCmodeldist']

modelvstrue(x1a, y1a, y1b, title="DS DISK FirstOvertone Modeldist/Truedist %", xlabel="True Dist", ylabel="Model Dist")

## Classical Cepheid Variables ##

# Fundamental BLG #

x1a = FundamentalblgcephF['Distance']
y1a = FundamentalblgcephF['SMCmodeldist']
y1b = FundamentalblgcephF['LMCmodeldist']

modelvstrue(x1a, y1a, y1b, title="CEPH BLG Fundamental Modeldist/Truedist %", xlabel="True Dist", ylabel="Model Dist")

# Fundamental DISK #

x1a = FundamentaldiskcephF['Distance']
y1a = FundamentaldiskcephF['SMCmodeldist']
y1b = FundamentaldiskcephF['LMCmodeldist']

modelvstrue(x1a, y1a, y1b, title="CEPH DISK Fundamental Modeldist/Truedist %", xlabel="True Dist", ylabel="Model Dist")

# FirstOvertone BLG #

x1a = FirstOvertoneblgcephF['Distance']
y1a = FirstOvertoneblgcephF['SMCmodeldist']
y1b = FirstOvertoneblgcephF['LMCmodeldist']

modelvstrue(x1a, y1a, y1b, title="CEPH BLG FirstOvertone Modeldist/Truedist %", xlabel="True Dist", ylabel="Model Dist")

# FirstOvertone DISK #

x1a = FirstOvertonediskcephF['Distance']
y1a = FirstOvertonediskcephF['SMCmodeldist']
y1b = FirstOvertonediskcephF['LMCmodeldist']

modelvstrue(x1a, y1a, y1b, title="CPEH DISK FirstOvertone Modeldist/Truedist %", xlabel="True Dist", ylabel="Model Dist")

###

### Percentage Error vs Distance Graphs ###

## Delta Scuti Variables ##

# Fundamental BLG #

x2a = FundamentalblgdsF['Distance']
y2a = (FundamentalblgdsF['SMCmodeldist']/FundamentalblgdsF['Distance'])*100
y2b = (FundamentalblgdsF['LMCmodeldist']/FundamentalblgdsF['Distance'])*100

errorfunction(x2a, y2a, y2b, title="DS BLG Fundamental Percentage Error % Vs True Distance", xlabel="True Distance", ylabel="Error% (Modeldist/Truedist)")

# Fundamental DISK #

x2a = FundamentaldiskdsF['Distance']
y2a = (FundamentaldiskdsF['SMCmodeldist']/FundamentaldiskdsF['Distance'])*100
y2b = (FundamentaldiskdsF['LMCmodeldist']/FundamentaldiskdsF['Distance'])*100

errorfunction(x2a, y2a, y2b, title="DS DISK Fundamental Percentage Error % Vs True Distance", xlabel="True Distance", ylabel="Error% (Modeldist/Truedist)")

# FirstOvertone BLG #

x2a = FirstOvertoneblgdsF['Distance']
y2a = (FirstOvertoneblgdsF['SMCmodeldist']/FirstOvertoneblgdsF['Distance'])*100
y2b = (FirstOvertoneblgdsF['LMCmodeldist']/FirstOvertoneblgdsF['Distance'])*100

errorfunction(x2a, y2a, y2b, title="DS BLG FirstOvertone Percentage Error % Vs True Distance", xlabel="True Distance", ylabel="Error% (Modeldist/Truedist)")

# FirstOvertone DISK #

x2a = FirstOvertonediskdsF['Distance']
y2a = (FirstOvertonediskdsF['SMCmodeldist']/FirstOvertonediskdsF['Distance'])*100
y2b = (FirstOvertonediskdsF['LMCmodeldist']/FirstOvertonediskdsF['Distance'])*100

errorfunction(x2a, y2a, y2b, title="DS DISK FirstOvertone Percentage Error % Vs True Distance", xlabel="True Distance", ylabel="Error% (Modeldist/Truedist)")

## Classical Cepheid Variables ##

# Fundamental BLG #

x2a = FundamentalblgcephF['Distance']
y2a = (FundamentalblgcephF['SMCmodeldist']/FundamentalblgcephF['Distance'])*100
y2b = (FundamentalblgcephF['LMCmodeldist']/FundamentalblgcephF['Distance'])*100

errorfunction(x2a, y2a, y2b, title="CEPH BLG Fundamental Percentage Error % Vs True Distance", xlabel="True Distance", ylabel="Error% (Modeldist/Truedist)")

# Fundamental DISK #

x2a = FundamentaldiskcephF['Distance']
y2a = (FundamentaldiskcephF['SMCmodeldist']/FundamentaldiskcephF['Distance'])*100
y2b = (FundamentaldiskcephF['LMCmodeldist']/FundamentaldiskcephF['Distance'])*100

errorfunction(x2a, y2a, y2b, title="CEPH DISK Fundamental Percentage Error % Vs True Distance", xlabel="True Distance", ylabel="Error% (Modeldist/Truedist)")

# FirstOvertone BLG #

x2a = FirstOvertoneblgcephF['Distance']
y2a = (FirstOvertoneblgcephF['SMCmodeldist']/FirstOvertoneblgcephF['Distance'])*100
y2b = (FirstOvertoneblgcephF['LMCmodeldist']/FirstOvertoneblgcephF['Distance'])*100

errorfunction(x2a, y2a, y2b, title="CEPH BLG FirstOvertone Percentage Error % Vs True Distance", xlabel="True Distance", ylabel="Error% (Modeldist/Truedist)")

# FirstOvertone DISK #

x2a = FirstOvertonediskcephF['Distance']
y2a = (FirstOvertonediskcephF['SMCmodeldist']/FirstOvertonediskcephF['Distance'])*100
y2b = (FirstOvertonediskcephF['LMCmodeldist']/FirstOvertonediskcephF['Distance'])*100

errorfunction(x2a, y2a, y2b, title="CEPH DISK FirstOvertone Percentage Error % Vs True Distance", xlabel="True Distance", ylabel="Error% (Modeldist/Truedist)")

###







