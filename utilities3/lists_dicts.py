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

# function to merge two dicts such that the result contains 
# all keys present in either dict and sum of counts of their corresponding values
from collections import Counter
def dict_sum(A,B):
    out = dict( Counter(A) + Counter(B))
    return out

# function to remove duplicates from a list
# input: any list ( nested or not) output: list without duplicates
def dedup(inp):
    seen = []
    for item in inp:
        if item not in seen:
            seen.append(item)
    return seen

### function to merge two dicts of lists of binary 
def merge_dict(A,B):

    # length of binary sequences in each dict
    lA = len(list(A.values())[0])
    lB = len(list(B.values())[0])
    total = lA+lB

    # perform set operations on keys
    g1 = set( A.keys())
    g2 = set( B.keys())
    all_g = g1 | g2
    g1a2 = g1 & g2
    g1n2 = g1 - g2
    g2n1 = g2 - g1

    # initialize dict to hold output
    out = { k:[ 0 for i in range(total)] for k in all_g }
    
    for g in all_g:
        if g in g1n2:
            out[g][:lA] = A[g]
        elif g in g2n1:
            out[g][lA:] = B[g]
        else:
            out[g] = A[g]+B[g]
    return out
    
# example
A = { 'g1':[0,0,1], 'g2':[0,1,2]}
B = { 'g2':[1,0,0,1], 'g3':[1,0,1,1]}
out = merge_dict(A,B)

# Generator to flatten an arbitrarily nested list
# source: https://stackoverflow.com/questions/10823877/what-is-the-fastest-way-to-flatten-arbitrarily-nested-lists-in-python
# author: hexparrot
def flatten(container):
    for i in container:
        if isinstance(i, (list,tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i

## convert a matrix of keys to a matrix of values
import numpy as np
def mat_k2v(mat,adict):
    a = mat
    a_flat = np.ndarray.flatten(a)  
    b = [ adict[x] if x in adict.keys() else '*' for x in a_flat]
    c = np.reshape( b, a.shape)
    return c

# load a multi-FASTA file into a dict
def multifasta_dict(infile):
    dct = {}
    with open(infile) as flob:
        data = [ i.strip('\n') for i in flob]
    # headers indices
    hxs = [ x for x,i in enumerate(data) if '>' in i]
    # no. of lines in data
    L = len(data)
    # no. of sequences in data
    S = len(hxs)
    for x in range(S):
        # use header as key
        key = data[ hxs[x]][1:]
        # skip the sequence if the key is already present in the dict
        if key in dct.keys(): continue
        # starting index of sequence
        xli = hxs[x]+1
        # ending index of sequence
        if x == S-1:
            xlf = L
        else:
            xlf = hxs[x+1]
        # extract sequence
        seq = ''.join(data[xli:xlf])
        # add to the dict
        dct[key] = seq
    return dct


