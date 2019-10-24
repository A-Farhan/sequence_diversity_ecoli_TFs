# require
import os, re
from Bio import AlignIO as align
from Bio.Alphabet import IUPAC
from scipy.special import binom
import numpy as np

class MSA(object):
    def __init__(self,fname,seqtype):
        # gene name
        self.gene = os.path.basename(fname).replace('.fasta','')
        # sequence type: nucletide(n) OR protein(p)
        self.seqtype = seqtype
        if self.seqtype=='p':
            # protein alphabets
            self.alphabets = IUPAC.protein.letters
            # multiple sequence alignment object from biopython
            self.aln = align.read(fname,'fasta',alphabet=IUPAC.protein)
        elif self.seqtype=='n':
            # nucleotide alphabets
            self.alphabets = IUPAC.unambiguous_dna.letters
             # multiple sequence alignment object from biopython
            self.aln = align.read(fname,'fasta',alphabet=IUPAC.unambiguous_dna)
        # if neither of the above, terminate
        else: raise Exception('wrong sequence type') 
        # numpy array form of the MSA
        self.array = np.array( [ list(rec) for rec in self.aln])
        # number of sequences and number of sites
        self.nseq, self.ncol = self.array.shape
        
        # indices of columns with gaps
        self.gaps = [ x for x in range(self.ncol) \
                        if '-' in self.array[:,x]]
        # no. of columns without gaps
        self.nogapcol = self.ncol - len(self.gaps)
        # indices of columns with differences but no gaps
        self.difs = [ x for x in range(self.ncol) if \
                        x not in self.gaps and \
                        len( set( self.array[:,x])) > 1]
        # no. of variant sites
        self.ndif = len(self.difs)

    # method to remove columns with gaps
    def remove_gapcols(self):   
        arr = np.delete( self.array, self.gaps, 1) 
        return arr

    # method to generate a pairwise distance matrix
    def pdistmat(self):
        # array after removing columns with gaps
        arr = self.remove_gapcols()
        # number of sequences
        N = self.nseq
        # number of columns Without gaps
        L = self.nogapcol
        # initialize a 2-D array with -1
        dmat = np.full( (N,N), -1, dtype = float)

        # get p-distances
        for x in range(0,N):
            # distance to self is zero
            dmat[x,x] = 0
            for y in range(x+1,N):
                # number of differing elements between two rows of the array
                nd = sum( arr[x,:] != arr[y,:])
                # p-distance: dividing above by no. of sites
                p = round(float(nd)/L,4)
                dmat[x,y] = p
        return dmat
    
    # method to calculate diversity(nuc or aa) based on above matrix
    def div(self):
        # pairwise distance matrix
        dmat = self.pdistmat()
        # filler entries (-1) replaced by zero
        dmat[ dmat == -1] = 0
        # so diversity could be calculated as below
        D = round( sum(sum(dmat))/binom( self.nseq, 2), 4)
        return D

## Function to compute nucleotide diversity given a 
# biopython MSA object 
def msa_nucdiv(aln):       
    """ This function computes nucleotide diversity given a \
            biopython MSA object """
    # numpy array form of the MSA
    array = np.array( [ list(rec) for rec in aln])
    # indices of columns with gaps
    gaps = [ x for x in range(array.shape[1]) if '-' in array[:,x]]
    # remove above columns
    array = np.delete( array, gaps, 1) 

    # number of sequences and number of sites
    N,L = array.shape
    # initialize a 2-D array with -1
    dmat = np.full( (N,N), -1, dtype = float)

    # get p-distances
    for x in range(0,N):
        # distance to self is zero
        dmat[x,x] = 0
        for y in range(x+1,N):
            # number of differing elements between two rows of the array
            nd = sum( array[x,:] != array[y,:])
            # p-distance: dividing above by no. of sites
            p = round(float(nd)/L,4)
            dmat[x,y] = p
    
    # filler entries (-1) replaced by zero
    dmat[ dmat == -1] = 0
    # so diversity could be calculated as below
    D = round( sum(sum(dmat))/binom( N, 2), 4)
    return D
