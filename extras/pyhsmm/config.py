from __future__ import division
from os.path import isdir

# Sets default directories for where to look for data and where to store the results
# datadir: where to load data subdirectories from (see goodtraces.txt)
# resultsdir: where to look for subdirectories in which to store pyhsmm output

def getdirs():
    if isdir('/Users/Steph/Documents/smFRETdata'):
         datadir = '/Users/Steph/Documents/smFRETdata/TracesOutput'
         resultsdir = '/Users/Steph/Documents/smFRETdata/pyhsmmOutput'
    elif isdir('/Users/Steph/Documents/SomethingElse'):
        datadir = '/Users/Steph/Documents/SomethingElse/something'
        resultsdir = '/Users/Steph/Documents/SomethingElse/etc'
    else:
        raise IOError('No default directories found!')

    return datadir, resultsdir


if __name__ == '__main__':
    for dir in getdirs():
        print dir

