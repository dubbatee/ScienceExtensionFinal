import pandas as pd
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u
import astropy.coordinates as coord
#from intragen import testingset
from intragen import raw1blgceph, raw1diskceph
import statistics
import numpy as np

### Querying the Bailer Jones Catalogue for distances ###

print('Beginning Ceph Query')

def cephquery(dataframe):

    # Create a new empty column for distances
    dataframe['Distance'] = None
    dataframe['SourceID'] = None

    # Iterate through each row in the DataFrame
    for index, row in dataframe.iterrows():
        ra = row['Ra']
        decl = row['Decl']

        try:
            # Connect to Vizier catalog
            v = Vizier(catalog="I/352", columns=['*', 'RA_ICRS', 'DE_ICRS'])

            # Query for nearby objects within 0.0001 degrees radius
            result = v.query_region(
                coord.SkyCoord(ra=ra, dec=decl, unit=(u.deg, u.deg), frame='icrs'),
                radius=0.0001 * u.deg
            )

            # Check for results
            if len(result) > 0:
                # If results found, extract distance from the first object
                distance = result[0]['rgeo'][0]
                dataframe.at[index, 'Distance'] = distance
                id = result[0]['Source'][0]
                dataframe.at[index, 'SourceID']= id
                
            else:
                # If no results, set distance to 0
                dataframe.at[index, 'Distance'] = 0
                dataframe.at[index, 'SourceID'] = 0

        except Exception as e:
            # Handle potential exceptions during Vizier query
            print(f"Error querying Vizier for row {index}: {e}")
            dataframe.at[index, 'Distance'] = 0  # Set distance to 0 for error cases
            dataframe.at[index, 'SourceID'] = 0

    return dataframe

### BULGE Classical Cepheid QUERY / 1st ROUND CLEANSING ### 

rawblgcephdistances = cephquery(raw1blgceph.copy())  # Avoid modifying original DataFrame
#print(rawblgcephdistances)

cleanblgceph = pd.DataFrame(columns = ['ID', 'mode', 'Ra', 'Decl', 'I', 'V', 'V-I', 'P1', 'Distance', 'SourceID']) #creating a new dataframe for the cleansed data
for index, row in rawblgcephdistances.iterrows(): #looping through raw data dataframe
    if np.float64(row['Distance']) != 0:
        cleanblgceph.loc[len(cleanblgceph.index)] = [row['ID'], row['mode'], row['Ra'], row['Decl'], row['I'], row['V'], row['V-I'], row['P1'], row['Distance'], row['SourceID']] # adds the star at the end of the dataframe (last index)
    else:
        continue

#print(cleanblgceph)
print('CHECK 3')


### DISK Classical Cepheid QUERY / 1st ROUND CLEANSING ### 

rawdiskcephdistances = cephquery(raw1diskceph.copy())  # Avoid modifying original DataFrame
print(rawdiskcephdistances)

cleandiskceph = pd.DataFrame(columns = ['ID', 'mode', 'Ra', 'Decl', 'I', 'V', 'V-I', 'P1', 'Distance', 'SourceID']) #creating a new dataframe for the cleansed data
for index, row in rawdiskcephdistances.iterrows(): #looping through raw data dataframe
    if np.float64(row['Distance']) != 0:
        cleandiskceph.loc[len(cleandiskceph.index)] = [row['ID'], row['mode'], row['Ra'], row['Decl'], row['I'], row['V'], row['V-I'], row['P1'], row['Distance'], row['SourceID']] # adds the star at the end of the dataframe (last index)
    else:
        continue

#print(cleandiskceph)
print('CHECK 4')

print('Querying complete, sending dataframes to intragen')

