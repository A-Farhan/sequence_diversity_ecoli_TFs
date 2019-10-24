# function to print elapsed time in hh:mm:ss format given starting time
from time import time
def timer(b=0):
    from math import modf
    if b==0: b=time()
    rts = round(time() - b)
    if rts < 60:
        out = '00:00:'+str(rts)
    else:
        rtm_frac,rtm_whole = modf(rts/60)
        rtm = rtm_whole
        rtms = round(rtm_frac*60)
        if rtm < 60:
            out = '00:'+str(rtm)+':'+str(rtms)
        else:
            rth_frac,rth_whole = modf(rtm/60)
            rth = rth_whole
            rthm = round(rth_frac*60)
            out = str(rth)+':'+str(rthm)+':00'
    return out

# function to open and load csv files into list
def readcsv(fl,sep=','):
    import csv
    with open(fl) as flob:
        rdob = csv.reader( filter(lambda x: x[0]!='#', flob), delimiter=sep )
        data = [i for i in rdob]
    return data

# function to write a list to csv file
def writecsv(fl,data,header=None,sep=',',comments='',mode='w'):
    import csv
    with open(fl,mode) as flob:
        flob.write(comments)
        wrob = csv.writer(flob,delimiter=sep)
        if header != None: wrob.writerow(header)
        for i in data: wrob.writerow(i)

## Function to check if the file exists and is non zero
import os
def filetest_s(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

# function to split data by a factor
# convert list of lists into a dict where keys come from 
# entries at a particular index in sublists
from itertools import groupby
def split_data( data, ix):
    keyfunc = lambda v:v[ix]
    data = sorted( data, key = keyfunc)  
    out = { k:[ i for i in g] for k, g \
            in groupby( data, keyfunc) }
    return out

# function to remove duplicates from a list
# input: any list ( nested or not) output: list without duplicates
def dedup(inp):
    seen = []
    for item in inp:
        if item not in seen:
            seen.append(item)
    return seen

# function that takes an input a nested list where some sublists
# may contain others and return a nested list where no sublist is
# contained into another
def rm_contained(alist):
    # length of the list 
    L = len(alist)
    out = []
    for x in range(L):
        # list of entries that contain all the 
        # elements of the current entry
        n = sum([ all( el in alist[y] for el in alist[x]) \
                    for y in range(L) if y != x])
        # add to the new list if no other entry contains this entry
        if n == 0:
            out.append(alist[x])
    return out

# function to select set of objects satisfying minimum distance threshold
# given a distance matrix
def selbydist( dm, cut, objects):
    # no. of objects
    R = len(objects)
    # initialize set of selected indices with the first one
    sel = {0}
    # empty set for the indices to be removed
    removed = set()
    # for every index, except the last
    for x in range(0,R-1):
        # if it is present in "removed", skip to next index
        if x in removed: continue
        # for every index following the present one
        for y in range(x+1,R):
            # if it is present in "removed", skip to next index
            if y in removed: continue
            # otherwise, get the distance between runs corresponding to these two indices
            d = dm[x,y]
            # if it is more than threshold, select the second index
            if d > cut:
                sel.add(y)
            # else, mask it from further consideration
            else:
                removed.add(y)
                # if the second index was already selected, remove it
                if y in sel:
                    sel.remove(y)
    # get runs corresponding to selected indices
    out = [ i for x,i in enumerate(objects) if x in sel]
    return out

# decorate log function to return None if input is zero
import math
def log(x):
    if x <= 0: 
        return None
    else:
        return math.log(x)

## Shannon's entropy ##
# require

# define python exception for negative probabilities:
class negative_probability(Exception):
    """Raised when a list of probabilities contains a negative value"""
    pass
class SumNot1(Exception):
    """Raised when sum of probabilities is not one"""
    pass

def shannons_entropy(prob_list,b=2):
    # check that there are no negative values
    if any([i<0 for i in prob_list]):
       raise negative_probability('P must be non-negative')
    # exclude zero-probabilities
    prob_list = [i for i in prob_list if i!=0]
    # check if sum of probabilities is one
    if round(sum(prob_list))!=1:
        raise SumNot1('All probabilities must sum to 1')
    H = -sum([p*math.log(p,b) for p in prob_list])
    return(H)

