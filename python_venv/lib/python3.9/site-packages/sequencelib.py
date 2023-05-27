# By Anders Gorm Pedersen, agpe@dtu.dk, Technical University of Denmark, Bioinformatics, 2012-2022

##############################################################################################################

"""Classes and methods for reading, analyzing, manipulating, and writing DNA and protein sequences"""

from math import ceil
from math import floor
from math import sqrt
from math import log
from math import log10
import math

from collections import Counter
import re
import string
import sys
import time
import random
import itertools
from io import StringIO
import copy
import os
import numpy as np
import Levenshtein as lv

##############################################################################################################
# Various functions used by methods, that do not fit neatly in any class
##############################################################################################################

def find_seqtype(seqsample):
    # Note: for small sequences it can occur that Protein seq contains only DNA_maxambig symbols
    # Perhaps look at frequencies instead? Other approaches?
    letters = set(seqsample) - set("-")   # Accepts most input types (list, string, set)
    if letters <= Const.Standard:
        return "standard"
    elif letters <= Const.DNA_maxambig:
        return "DNA"
    elif letters <= Const.Protein_maxambig:
        return "protein"
    elif letters <= Const.ASCII:
        return "ASCII"     # Maybe drop as a type?
    else:
        raise SeqError("Unknown sequence type. Unrecognized symbols: {}".format(list(letters - Const.ASCII)))

##############################################################################################################

def seqtype_attributes(seqtype):
    """Returns alphabet and ambigsymbols for given seqtype"""

    # Python note: hack. Should rethink logic around seqtypes and refactor
    # They can be set too many places at both container and element level
    if seqtype == "standard":
        return (Const.Standard, set())
    elif seqtype == "DNA":
        return (Const.DNA_maxambig, Const.DNA_maxambig - Const.DNA)
    elif seqtype == "protein":
        return (Const.Protein_maxambig, Const.Protein_maxambig - Const.Protein)
    elif seqtype == "ASCII":
        return (Const.ASCII, set())
    else:
        raise SeqError("Unknown sequence type: {}".format(seqtype))

##############################################################################################################

def indices(mystring, substring):
    """Helper function that finds indices of substring in string. Returns as set"""
    result = set()
    offset = -1
    while True:
        try:
            offset = mystring.index(substring, offset+1)
        except ValueError:
            return result
        result.add(offset)

##############################################################################################################

def remove_comments(text, leftdelim, rightdelim):
    """Takes input string and strips away commented text, delimited by 'leftdelim' and 'rightdelim'.
        Also deals with nested comments."""

    # NOTE: only deals with block comments at present
    # Python note: maybe this is too general. Will I ever use multichar delims?
    def wordsoverlap(w1, w2):
        for i in range(1, len(w2)):
            if w1.startswith(w2[i:]):
                return True
        for i in range(1, len(w1)):
            if w2.startswith(w1[i:]):
                return True
        return False

    if leftdelim == rightdelim:
        raise SeqError("Left and right delimiters are identical")
    elif leftdelim in rightdelim:
        raise SeqError("Left delimiter is substring of right delimiters")
    elif rightdelim in leftdelim:
        raise ExcepSeqErrortion("Right delimiter is substring of left delimiters")
    elif wordsoverlap(leftdelim, rightdelim):
        raise SeqError("Right and left delimiters overlap")

    # Preprocess delims for use in re etc
    leftdelim = re.escape(leftdelim)
    rightdelim = re.escape(rightdelim)

    # Construct sorted list of tuples of the form [(0, 'start'), (5, 'stop'), (7, 'start'), ...]
    delimlist = [(match.start(), match.end(), "start") for match in re.finditer(leftdelim, text)]
    # If text contains no starts (=> no comments): return un-altered text
    if not delimlist:
        return text
    else:
        delimlist.extend([(match.start(), match.end(), "stop") for match in re.finditer(rightdelim, text)])
        delimlist.sort()

    # Traverse text; along the way copy text not inside comment-delimiter pairs.
    # Use stack ("unmatched_starts") to keep track of nesting
    unmatched_starts = 0
    prevpos = 0
    processed_text = []
    for (match_start, match_end, match_type) in delimlist:
        if match_type == "start":
            unmatched_starts += 1
            if unmatched_starts == 1:                               # Beginning of new comment region
                processed_text.append(text[prevpos:match_start])
        elif match_type == "stop":
            unmatched_starts -= 1
            if unmatched_starts == 0:                               # End of comment region
                prevpos = match_end
            elif unmatched_starts == -1:                            # Error: more right delims than left delims
                raise Exception("Unmatched end-comment delimiter. Context: '{}'".format(text[prevpos-10:prevpos+10]))

    # Add final block of text if relevant (i.e., if text does not stop with rightdelim), return processed text
    if prevpos < len(text):
        processed_text.append(text[prevpos:])
    return "".join(processed_text)

#############################################################################################################

def make_sparseencoder(alphabet, padding="X"):
    """Returns function that can sparse-encode strings in specified alphabet"""

    # This function uses a "closure" to create and return a function that can later be used
    # for sparse-encoding strings in the specified alphabet.
    # The purpose is to avoid having to build the translation dictionary every time the
    # encoder function is run (or alternatively to compute it preemptively on module import)

    # Check that padding symbol is not also present in alphabet
    if padding in alphabet:
        raise SeqError("Sparse-encoding error: padding symbol can't also be present in alphabet")

    # Build translation dictionary for specified alphabet.
    # This will be available to enclosed (and returned) function
    alphabet = sorted(alphabet)
    zerolist = [0] * len(alphabet)
    transdict = {}
    for i,residue in enumerate(alphabet):
        vec = zerolist[:]
        vec[i] = 1
        transdict[residue] = vec
    transdict[padding] = zerolist[:]

    # Enclosed function that will be returned and that can do sparse-encoding of input string
    def sparse_encoder(sequence_string):
        """Sparse-encodes input string. Output is numpy array. Padding encoded as all zeroes:
            array([0,0,0,1,1,0,0,0,0,0,0,0])"""
        sparse_list = []
        for residue in sequence_string:
            sparse_list.extend(transdict[residue])
        sparse_seq = np.array(sparse_list)
        return sparse_seq

    return sparse_encoder   # Closure: return pointer to function that knows about transdict

#############################################################################################
#############################################################################################

class SeqError(Exception):
    """Seqlib exceptions"""

    def __init__(self, errormessage):
        self.errormessage = errormessage

    def __str__(self):
        return self.errormessage

#############################################################################################
#############################################################################################

class Const(object):
    """Global constants used by seqlib"""

    # Alphabets
    Standard = set("0123456789")
    DNA = set("ACGT")
    DNA_minambig = set("ACGTN")
    DNA_typicalambig = set("ACGTYRN")
    DNA_maxambig = set("ACGTURYMKWSBDHVN")
    Protein = set("ACDEFGHIKLMNPQRSTVWY")
    Protein_minambig = set("ACDEFGHIKLMNPQRSTVWYX")
    Protein_maxambig = set("ACDEFGHIKLMNPQRSTVWYBZX")
    ASCII = set(string.ascii_uppercase + string.digits + " ,._")

#############################################################################################
#############################################################################################

class Sequence(object):
    """Baseclass for sequence classes"""

    def __init__(self, name, seq, annotation, comments, check_alphabet, degap):

        # NOTE: additional fields can be added later (e.g., feature list, etc.)
        self.name = name
        self.seq = seq.upper()
        self.annotation = annotation
        self.comments = comments
        self.check_alphabet = check_alphabet
        self.degap = degap

        # If requested: remove gap characters
        if self.degap:
            self.remgaps()

        # If requested: use set arithemtic to check whether seq contains any non-alphabet characters
        if check_alphabet:
            seqset = set(self.seq.upper())
            non_alphabet = seqset - self.alphabet - set("-.")

            if len(non_alphabet)>0:
                raise SeqError("Unknown symbols in sequence %s: %s" % (name, list(non_alphabet)))

    #######################################################################################

    def __eq__(self, other):
        """Implements equality check between sequences"""

        # Note: Here I take two sequences to be identical if they have same gapfree seq (name is ignored)

        if self.seq.replace("-", "") == other.seq.replace("-", ""):
            return True
        else:
            return False

    #######################################################################################

    def __ne__(self, other):
        """Implements inequality checking between sequences"""

        # Python note: __eq__ does not cover use of != operator, so this has to be explicitly implemented
        if self.seq.replace("-", "") != other.seq.replace("-", ""):
            return True
        else:
            return False

    #######################################################################################

    def __len__(self):
        """Implements the len() operator for sequences"""
        return len(self.seq)

    #######################################################################################

    def __getitem__(self, index):
        """Implements indexing, slicing, and iteration for sequences"""

        return self.seq[index]

    #######################################################################################

    def __setitem__(self, index, residue):
        """Allows change of single residues by using dict-like syntax: seq[5] = 'A' for instance"""

        newseq = []
        newseq.append(self.seq[:index])     # Sequence string up to right before new residue
        newseq.append(residue.upper())
        newseq.append(self.seq[index + 1:]) # Original sequence string from right after new residue, to end
        self.seq = "".join(newseq)

    #######################################################################################

    def __str__(self):
        """String (fasta) representation of sequence"""

        return self.fasta()

    #######################################################################################

    def copy_seqobject(self):
        """Customized deep copy of sequence object"""

        seqcopy = self.__class__(self.name, self.seq, self.annotation, self.comments)
        return seqcopy

    #######################################################################################

    def rename(self, newname):
        """Changes name of sequence"""

        self.name = newname

    #######################################################################################

    def subseq(self, start, stop, slicesyntax=True, rename=False):
        """Returns subsequence as sequence object of proper type. Indexing can start at zero or one"""

        # Rename if requested. Note: naming should be according to original indices (before alteration)
        if rename:
            name = self.name + "_%d_%d" % (start, stop)
        else:
            name = self.name

        # If slicesyntax is False: indexing starts at 1, and stop is included (and hence left unchanged)
        if not slicesyntax:
            start -= 1

        # Sanity check: are requested subsequence indices within range of seq?
        if start < 0 or stop > len(self):
            raise SeqError("Requested subsequence (%d to %d) exceeds sequence length (%d)" % (start, stop, len(self)))

        seq = self.seq[start:stop]
        comments = self.comments
        annotation = self.annotation[start:stop]

        subseq = self.__class__(name, seq, annotation, comments)    # Create new object of same class as self
                                                                    # (don't know which sequence type we're in)
        return(subseq)

    #######################################################################################

    def subseqpos(self, poslist, namesuffix=None):
        """Returns subsequence containing residues on selected positions as
        sequence object of proper type. Indexing can start at zero or one"""

        subseq = self.__class__(self.name, self.seq,
                                self.annotation, self.comments)    # Create new object of same class as self
                                                                   # (don't know which sequence type we're in)
        subseq.indexfilter(poslist)
        if namesuffix:
            subseq.name = subseq.name + namesuffix

        return(subseq)

    #######################################################################################

    def appendseq(self, other):
        """Appends seq from other to end of self. Name from self is retained"""
        self.seq += other.seq
        self.annotation += other.annotation
        self.comments += " " + other.comments

    #######################################################################################

    def prependseq(self, other):
        """Prepends seq from other before start of self. Name from self is retained"""
        self.seq = other.seq + self.seq
        self.annotation = other.annotation + self.annotation
        self.comments = other.comments + " " + self.comments

    #######################################################################################

    def windows(self, wsize, stepsize=1, l_overhang=0, r_overhang=0, padding="X", rename=False):
        """Returns window iterator object"""

        # Function that returns window iterator object, which can be used to iterate over windows on sequence
        # Iteration returns objects of the proper sequence type (same as parent). The rename option is like for subseq()
        # Main trick here is a class which keeps track of its enclosing Sequence object ("parent"),
        # while at the same time keeping track of which window the iteration has gotten to
        # stepsize: consecutive windows are stepsize residues apart
        # l_overhang, r_overhang: Iteration includes windows that "fall off the edge" of the sequence
        # l_overhang specifies the location of the leftmost window (how many window positions on the left are outside the seq)
        # r_overhang does the same for rightmost window (how many window positions are on the right of seq for rightmost window)
        # padding: use this character to pad the overhangs (all non-sequence positions are set to this)

        class Window_iterator:
            """Window iterator object"""

            def __init__(self, parent, wsize, stepsize, l_overhang, r_overhang, padding, rename):
                if l_overhang > wsize or r_overhang > wsize:
                    raise SeqError("l_overhang and r_overhang must be smaller than wsize")
                if wsize > len(parent):
                    raise SeqError("wsize must be smaller than length of sequence")

                self.i = 0 - l_overhang
                self.wsize = wsize
                self.stepsize = stepsize
                self.l_overhang = l_overhang
                self.r_overhang = r_overhang
                self.padding = padding
                self.length = len(parent)
                self.right_end = self.length + r_overhang
                self.parent = parent
                self.rename = rename

            def __iter__(self):
                return self

            def __next__(self):
                if self.i + self.wsize <= self.right_end:
                    start = self.i
                    stop = start + self.wsize
                    self.i += self.stepsize
                    if start < 0:
                        window = self.parent.subseq(0, stop, rename=self.rename)
                        window.seq = abs(start) * self.padding + window.seq  # Add the required number of padding chars
                    elif start >= 0 and stop <= self.length:
                        window = self.parent.subseq(start, stop, rename=self.rename)
                    elif start > 0 and stop > self.length:
                        window = self.parent.subseq(start, self.length, rename=self.rename)
                        window.seq = window.seq + (stop - self.length) * self.padding # Add the required number of padding chars
                    else:
                        raise SeqError("Execution really should never get to this line")
                    return window
                else:
                    raise StopIteration

        # In call below "self" points to enclosing Sequence object ("parent")
        return Window_iterator(self, wsize, stepsize, l_overhang, r_overhang, padding, rename)

    #######################################################################################

    def remgaps(self):
        """Removes gap characters from sequence"""

        self.seq = self.seq.replace("-", "")

    #######################################################################################

    def shuffle(self):
        """Returns shuffled sequence as object of proper subclass"""
        seqlist = list(self.seq)
        random.shuffle(seqlist)
        shufseq = "".join(seqlist)
        return self.__class__(name=self.name + "_shuffled", seq=shufseq)  # Works regardless of class

    #######################################################################################

    def indexfilter(self, keeplist):
        """Discards symbols at positions that are not in keeplist (seq changed in place)"""

        s = self.seq
        seqlist = [s[i] for i in keeplist]
        self.seq = "".join(seqlist)

        # If Sequence object contains annotation: also apply filter to annotation sequence
        if self.annotation:
            s = self.annotation
            annotlist = [s[i] for i in keeplist]
            self.annotation = "".join(annotlist)


    #######################################################################################

    def seqdiff(self, other, zeroindex=True):
        """Returns list of tuples summarizing difference between seqs in self and other:
            (site, residue1, residue2)
        Option zeroindex=False causes site numbering to start at 1 (otherwise 0)
        """

        difflist = []
        for i, (res1, res2) in enumerate(zip(self.seq, other.seq)):
            if res1 != res2:
                if zeroindex:
                    difflist.append((i, res1, res2))
                else:
                    difflist.append((i+1, res1, res2))
        return difflist
        # Python note: perhaps some of these returned lists should be generators instead

    #######################################################################################

    def hamming(self, other):
        """Directly observable distance between self and other (absolute count of different positions)"""

        # Python note: completely outsourced to Levenshtein library which is blazing fast!
        return lv.hamming(self.seq, other.seq)

    #######################################################################################

    def hamming_ignoregaps(self, other):
        """Like hamming distance, but positions with gaps are ignored"""

        # Note: using this to build distance matrix will mean that for different pairs of sequences
        # a different set of alignment positions are included in the distance measure
        diffs = self.hamming(other)
        gap_indices_self = indices(self.seq, "-")
        gap_indices_other = indices(other.seq, "-")
        n_dont_count = len(gap_indices_self ^ gap_indices_other)

        return diffs - n_dont_count

    #######################################################################################

    def pdist(self, other):
        """Directly observable distance between self and other (in differences per site)"""

        return self.hamming(other)/len(self.seq)

    #######################################################################################

    def pdist_ignoregaps(self, other):
        """Like pdist distance, but positions with gaps are ignored"""

        # Note: using this to build distance matrix will mean that for different pairs of sequences
        # a different set of alignment positions are included in the distance measure

        return self.hamming_ignoregaps(other)/len(self.seq)

    #######################################################################################

    def residuecounts(self):
        """Returns dictionary with counts of residues for single seq. {letter:count}"""

        return Counter(self.seq)

    #######################################################################################

    def composition(self, ignoregaps=True):
        """Returns dictionary with composition for single seq. {letter:(count,freq)}"""

        countdict = self.residuecounts()
        if ignoregaps:
            alphabet = set(countdict.keys()) - set("-")
            length = len(self) - countdict["-"]
        else:
            alphabet = set(countdict.keys())
            length = len(self)

        compdict = {}
        for residue in alphabet:
            count = countdict[residue]
            freq = count/length
            compdict[residue] = (count, freq)
        return compdict

    #######################################################################################

    def fasta(self, width=60, nocomments=False):
        """Return fasta-formatted sequence as a string"""

        if not nocomments and self.comments:
            fasta = [">%s %s\n" % (self.name, self.comments)]
        else:
            fasta = [">%s\n" % (self.name)]

        # Print seq, fold line at "width" characters
        numlines = int(ceil(len(self)/float(width)))
        for i in range(numlines):
            fasta.append(self[i*width:i*width+width])
            fasta.append("\n")
        del fasta[-1]  # Remove trailing newline
        return "".join(fasta)

    #######################################################################################

    def how(self, width=80, nocomments=False):
        """Return how-formatted sequence as a string"""

        # NOTE: HOW requires a fixed width of 80. Should perhaps remove width-option from function

        if self.comments and not nocomments:
            how = ["{length:>6} {name} {com}\n".format(length=len(self.seq), name=self.name, com=self.comments)]
        else:
            how = ["{length:>6} {name}\n".format(length=len(self.seq), name=self.name)]

        # If sequence has no annotation: construct default annotation field consisting of all dots (same len as seq)
        # Note: I should perhaps check that annotation has same length as seq when present?
        if not self.annotation:
            self.annotation = "." * len(self)

        # Print seq and annotation, fold lines at "width" characters
        numlines = int(ceil(len(self)/float(width)))
        for i in range(numlines):
            how.append(self[i*width:i*width+width])
            how.append("\n")
        for i in range(numlines):
            how.append(self.annotation[i*width:i*width+width])
            how.append("\n")

        del how[-1]  # Remove trailing newline
        return "".join(how)

    #######################################################################################

    def gapencoded(self):
        """Returns gap-encoding of sequence (gap=1, nongap=0) as a string"""

        # Note: mostly useful in mrbayes analyses where one wants to model gap changes using binary model
        allsymbols = "".join(Const.DNA_maxambig | Const.Protein_maxambig | Const.ASCII)
        transtable = str.maketrans("-" + allsymbols, "1" + "0" * len(allsymbols))
        return self.seq.translate(transtable)

    #######################################################################################

    def tab(self, nocomments=False):
        """Returns tab-formatted sequence as a string:
            name TAB seq TAB annotation TAB comments
            If no annotation is present, comments are placed after two TABs if present"""

        if (nocomments):
            tabstring = "%s\t%s\t%s" % (self.name, self.seq, self.annotation)
        else:
            tabstring = "%s\t%s\t%s\t%s" % (self.name, self.seq, self.annotation, self.comments)
        return tabstring

    #######################################################################################

    def raw(self):
        """Returns raw-formatted sequence as a string"""

        # Note: RAW format consists of one sequence per line, with no other information
        return(self.seq)

#############################################################################################
#############################################################################################

class DNA_sequence(Sequence):
    """Class containing one DNA sequence, and corresponding methods"""

    # Implementation note: alphabet should really be set in a more principled manner

    def __init__(self, name, seq, annotation="", comments="", check_alphabet=False, degap=False):
        self.seqtype="DNA"
        self.alphabet = Const.DNA_maxambig
        self.ambigsymbols = self.alphabet - Const.DNA
        Sequence.__init__(self, name, seq, annotation, comments, check_alphabet, degap)

##    #######################################################################################
##
## NOTE: Consider whether this is a productive idea for alignments, if so reimplement at some point
##
##
##    def seq_to_intlist(self):
##        """Converts sequence to a list of integers, useful for rapid alignment"""
##
##        trans_dict = {'A':0, 'C':1, 'G':2, 'T':3}
##        intlist = [(trans_dict[letter]) for letter in self.seq]
##        return intlist
##
    #######################################################################################

    def revcomp(self):
        """Returns reverse complement sequence as DNA sequence object - original is unchanged"""

        comptable = str.maketrans("ATUGCYSWKMBDHVN", "TAACGRSWMKVHDBN")
        revcompseq = self.seq.translate(comptable)
        revcompseq = revcompseq[::-1]   # Bit of a hack: Empty slice (returns all) with stride -1 (reverses)
        name = self.name + "_revcomp"
        return DNA_sequence(name, revcompseq)

    #######################################################################################

    def translate(self):
        """Translates DNA sequence to protein. Returns Protein_sequence object"""

        # All unknown triplets are translated as X (ambiguity symbol "N" handled intelligently)
        gencode =  { 'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                     'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                     'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
                     'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                                 'TCN': 'S',

                     'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                     'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                     'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                     'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                     'CTN': 'L', 'CCN': 'P',             'CGN': 'R',

                     'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                     'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                     'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
                     'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
                                 'ACN': 'T',

                     'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                     'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                     'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                     'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G',
                     'GTN': 'V', 'GCN': 'A',             'GGN': 'G'
                }

        protseqlist = []
        for i in range(0, len(self), 3):
            triplet = self[i:i+3]
            aa = gencode.get(triplet, "X")
            protseqlist.append(aa)
        protseq = "".join(protseqlist)
        return Protein_sequence(self.name, protseq)

#############################################################################################
#############################################################################################

class Protein_sequence(Sequence):
    """Class containing one protein sequence, and corresponding methods"""

    def __init__(self, name, seq, annotation="", comments="", check_alphabet=False, degap=False):
        self.seqtype="protein"
        self.alphabet = Const.Protein_maxambig
        self.ambigsymbols = self.alphabet - Const.Protein
        Sequence.__init__(self, name, seq, annotation, comments, check_alphabet, degap)

    #######################################################################################

    def shuffle(self):
        """Returns shuffled sequence as protein sequence object - original is unchanged"""

        shufseq = Sequence.shuffle(self)
        name = self.name + "_shuffled"
        return Protein_sequence(name, shufseq)

#############################################################################################
#############################################################################################

class ASCII_sequence(Sequence):
    """Class containing one sequence containing ASCII letters, and corresponding methods"""

    # Python note: will this ever be used???

    def __init__(self, name, seq, annotation="", comments="", check_alphabet=False, degap=False):
        self.seqtype="ASCII"
        self.alphabet = Const.ASCII
        self.ambigsymbols = set()       # Empty set. Attribute used by some functions. Change?
        Sequence.__init__(self, name, seq, annotation, comments, check_alphabet, degap)

    #######################################################################################

    def shuffle(self):
        """Returns shuffled sequence as ASCII sequence object - original is unchanged"""

        shufseq = Sequence.shuffle(self)
        name = self.name + "_shuffled"
        return ASCII_sequence(name, shufseq)

#############################################################################################
#############################################################################################

class Standard_sequence(Sequence):
    """Sequence containing standard data-type (e.g., morphological)"""

    def __init__(self, name, seq, annotation="", comments="", check_alphabet=False, degap=False):
        self.seqtype="standard"
        self.alphabet = Const.Standard
        self.ambigsymbols = set()       # Empty set. Attribute used by some functions. Change?
        Sequence.__init__(self, name, seq, annotation, comments, check_alphabet, degap)

#############################################################################################
#############################################################################################

class Contig(object):
    """Class corresponding to contig, i.e. assembly of several reads (sequences).
    Has methods for checking Contig overlap, assembly, etc."""

    count = 0   # Used to keep track of number of instances of class, for individual naming of contigs

    #######################################################################################

    def __init__(self, sequence):
        """Initialises Contig object with one Sequence object"""

        self.__class__.count += 1
        self.name = "contig_{:04d}".format(self.__class__.count)
        contigseq = sequence.seq
        self.assembly = sequence.__class__(self.name, contigseq)   # Sequence object containing assembled sequence of entire Contig
        readname = sequence.name
        self.readdict = {}
        self.readdict[readname] = sequence
        self.readdict[readname].startpos = 0       # I am adding "startpos" and "stoppos" attributes to the Sequence object in readdict
        self.readdict[readname].stoppos = len(sequence.seq)

    #######################################################################################

    def findoverlap(self, other, minoverlap=0):
        """Checks if two Contig objects overlap (one contained in other, or 3' of one = 5' of other).
        Minoverlap: overlap has to be at least this large to call a match"""

        alen, blen = len(self.assembly.seq), len(other.assembly.seq)
        minlen = min(alen, blen)

        # Check if one sequence is contained in other sequence
        if alen > blen:
            found = self.assembly.seq.find(other.assembly.seq)
            if found != -1:
                astart, astop = found, found + blen
                bstart, bstop = 0, blen
                size = blen
                return  astart, astop, bstart, bstop, size
        else:
            found = other.assembly.seq.find(self.assembly.seq)
            if found != -1:
                astart, astop = 0, alen
                bstart, bstop = found, found + alen
                size = alen
                return  astart, astop, bstart, bstop, size

        # Check if 3'end of self overlaps with 5' end of other
        bestmatch = None
        bestlen = 0
        k = minlen
        while k >= minoverlap:
            if self.assembly.seq[-k:] == other.assembly[:k]:
                bestlen = k
                astart, astop = alen - k, alen
                bstart, bstop = 0, k
                size = k
                bestmatch = astart, astop, bstart, bstop, size
                break
            else:
                k -= 1

        # Now check if there is a better match where 5' end of self overlaps with 3' end of other
        k = minlen
        while k >= max(minoverlap, bestlen + 1):
            if (self.assembly.seq[:k] == other.assembly[-k:]):
                bestlen = k
                astart, astop = 0, k
                bstart, bstop = blen - k, blen
                size = k
                bestmatch = astart, astop, bstart, bstop, size
                break
            else:
                k -= 1

        return bestmatch

    #######################################################################################

    def merge(self, other, overlap):
        """Merges Contig object with other Contig given overlap info from findoverlap.
           Keeps track of positions of all reads from both Contigs"""

        # "overlap" is result tuple from .findoverlap method
        # b == 0: a is more upstream, so assembly starts with a (vice versa for a == 0)
        # bstop < bhi: b has extra part downstream of overlap that should be appended to assembly
        # Update read positions with respect to assembly when merging
        astart, astop, bstart, bstop, overlaplen = overlap
        alen = len(self.assembly)
        blen = len(other.assembly)

        # First part of assembly is sequence a (self)
        if bstart == 0:
            # Sequence b (other) has extra part downstream that needs to be appended to assembly
            if bstop < blen:
                downstreamseq = other.assembly.subseq(bstop, blen)
                self.assembly.appendseq(downstreamseq)

            # Add reads from b to readlist in a, and update read coordinates:
            # b added to end, so all reads in a are OK, but all reads in b were shifted by astart
            for readname in other.readdict:
                self.readdict[readname] = other.readdict[readname]
                self.readdict[readname].startpos += astart
                self.readdict[readname].stoppos += astart

        # First part of assembly is sequence b (other)
        # Sequence b must have extra part upstream that needs to be prepended to assembly (since bstart != 0)
        elif astart == 0:
            upstreamseq = other.assembly.subseq(0, bstart)
            self.assembly.prependseq(upstreamseq)

            # Correct coordinates for reads in readlist of a: these have been shifted by bstart
            for readname in self.readdict:
                self.readdict[readname].startpos += bstart
                self.readdict[readname].stoppos += bstart

            # Add reads from b to readlist in a. These have correct coordinates since b added at front
            self.readdict.update( other.readdict )

    #######################################################################################

    def regions(self, slicesyntax=True, join_samename=False):
        """Divides contig into regions based on read overlaps:
                region boundaries are created in all locations where a read starts or stops.
        Returns list of regions sorted by start pos. Information for a single region is stored in a special class and includes:
            start, stop pos of region (in contig coordinates)
            list of reads in region
                for each read: name, and start + stop pos (in read-coordinates)
        """

        #######################################################################################

        class Region(object):
            """Struct for keeping region-related info (used only in this method)"""

            def __init__(self, start, stop, read_list, read_dict):
                self.contig_start = start
                self.contig_stop = stop
                self.name_list = []
                self.readinfo_list = []
                for origname in read_list:
                    read = read_dict[origname]
                    read_start = self.contig_start - read.startpos   # Position in read sequence coordinates
                    read_stop = self.contig_stop - read.startpos
                    posname = "{}_{}_{}".format(origname, read_start, read_stop)
                    self.name_list.append(posname)
                    self.readinfo_list.append( (origname, read_start, read_stop) )
                self.name_list.sort()

            def __repr__(self):
                return "Start: {}\tStop:  {}\tReads:  {}".format(self.contig_start, self.contig_stop, self.readinfo_list)

            def get_align_tuple(self):
                return(self.readinfo_list[0])            # original read subsection chosen to represent this region (for seqconverter)

        #######################################################################################

        # First, get an overview over where the boundaries between charsets should be (=where reads start or stop)
        read_boundary_dict = {}
        read_boundary_set = set()

        for read in self.readdict.values():
            read_boundary_set.add(read.startpos)        # set automatically uniquefies values
            read_boundary_set.add(read.stoppos)

        read_boundary_list = list(read_boundary_set)
        read_boundary_list.sort()

        # For each boundary position keep track of which reads that are in the following region
        # => reads where the current boundary position is >= read.startpos AND < read.stoppos
        # Collect info in read_boundary_dict (key=boundary position, value = list of read names)
        for boundary in read_boundary_list[:-1]:
            read_boundary_dict[boundary] = []
            for read in self.readdict.values():
                if read.startpos <= boundary < read.stoppos:
                    read_boundary_dict[boundary].append(read.name)

        # Create Region object for each item in read_boundary_dict
        region_list = []
        for i in range( 0, len(read_boundary_list) - 1 ):
            start = read_boundary_list[i]
            stop = read_boundary_list[i+1]
            readlist = read_boundary_dict[start]
            reg = Region(start, stop, readlist, self.readdict)
            region_list.append(reg)

        return region_list

#############################################################################################
#############################################################################################

class Read_assembler(object):
    """Class for assembling reads into contigs"""

    #######################################################################################

    def __init__(self, seqset=None):
        """Initialises set of reads from Seq_set (or other iterable containing Sequence objects).
        If Seq_set is not provided then individual seqs can be added later"""

        self.contigdict = {}
        self.contiglist = None

        if seqset:
            for seq in seqset:
                c = Contig(seq)
                self.contigdict[c.name] = c

    #######################################################################################

    def addseq(self, seq):
        """For adding one sequence at a time to Read_assembler"""

        c = Contig(seq)
        self.contigdict[c.name] = c

    #######################################################################################

    def assemble(self, minoverlap=-1):
        """Assembles reads into as few contigs as possible, and returns list of Contig objects.
        minoverlap: minimum length of overlap required to merge reads.
        Special value minoverlap = -1: set value automatically based on size of data set"""

        # If reads were already assembled, return again (hmm - possibly wrong way to think of this...)
        if self.contiglist:
            return self.contiglist
        else:
            self.contiglist = []

        # Set minoverlap automatically if requested such that expectation of no. overlap in set is < 1
        if minoverlap == -1:
            l = 0
            for contig in self.contigdict.values():
                l += len( contig.assembly )
            minoverlap = ceil( log( l ) / log( 4 ) )    # Previous self should remind current self why this is clever...
        while len ( self.contigdict ) > 0:

            c1name, c1 = self.contigdict.popitem()
            contigfinished = False

            while not contigfinished:
                foundmatch = False
                for contigname in self.contigdict:
                    c2 = self.contigdict[contigname]
                    overlap = c1.findoverlap(c2, minoverlap)
                    if overlap:
                        foundmatch = True
                        c1.merge(c2, overlap)
                        break

                if foundmatch:
                    self.contigdict.pop(contigname)
                else:
                    contigfinished = True

            self.contiglist.append(c1)

        # Sort list of contigs based on contig name (i.e., based on order in which they were created since contig names are numbered)
        # Using decorate-sort-undecorate
        decorated = [(c.name, c) for c in self.contiglist]
        decorated.sort()
        self.contiglist = [d[1] for d in decorated]
        return self.contiglist

#############################################################################################
#############################################################################################

class Sequences_base(object):
    """Baseclass for classes that contain a number of sequences - don't instantiate!"""

    # Implementation note: This baseclass contains methods and attributes that are
    # common to all sequence collections. Subclasses then implement methods and attributes
    # that are specific to aligned or un-aligned sequence collections.
    # This class should NOT be instantiated.

    #######################################################################################

    def __init__(self, name=None, seqtype=None):
        self.name = name            # Label for set of sequences
        self.seqdict = {}           # seqdict format: {seqname:seqobject, ...}
        self.seqnamelist = []       # seqnamelist format: [name1, name2, ...]
        if seqtype is None:
            self.seqtype = None
        else:
            # Python note: seqtype perhaps should be object with separate attributes for
            # seqtype, alphabet and ambigsymbols!
            # Also: I am setting seqtype in too many different places (Sequence objects
            # for instance). Rethink and refactor!!
            self.alphabet, self.ambigsymbols = seqtype_attributes(seqtype)

    #######################################################################################

    def __len__(self):
        """Implements the len() operator for sequence sets (=number of seqs)"""
        return len(self.seqnamelist)

    #######################################################################################

    def __getitem__(self, index):
        """Implements indexing, slicing, and iteration for sequences in seqlist.
        Integer index returns Sequence object
        Slice returns subset (Seq_set or Seq_alignment object depending on what subclass we are in).
        Tuple ([row, column] indexing) returns subsequence of selected sequence(s)"""

        # Python note: "for loops expect that an IndexError will be raised for illegal indexes
        # to allow proper detection of the end of the sequence."
        # Seems I should have that here?
        if isinstance(index, int):
            return self.seqdict[self.seqnamelist[index]]
        elif isinstance(index, slice):
            return self.subset(self.seqnamelist[index])
        elif isinstance(index, tuple):  # Do I ever use this?
            i, j = index
            s = self[i]
            if isinstance(j, slice):
                return(s.subseq(start=j.start, stop=j.stop))
            elif isinstance(j, int):
                return(s.subseq(start=j, stop=j))
            else:
                raise SeqError("Column index for a set of sequences must be either integer or slice")
        else:
            raise SeqError("A set of sequences must be indexed using either integer, slice, or tuple ([row, column])")

    #######################################################################################

    def __eq__(self, other):
        """Implements equality check between sequence collections"""

        # Two sequence collections are identical if (1) they have same size, and (2) every sequence
        # in collection has a match in collection 2. Note: sequence names can be different
        if len(self) != len(other):
            return False

        else:
            # For each seq in 1, compare to seqs in 2.
            # If match found: move on to next seq in 1. If no match: return False
            for seq1 in self:
                for seq2 in other:
                    if seq1 == seq2:
                        break
                return False

        # If we fell off loop: collections must be identical. Return True
        return True

    #######################################################################################

    def __ne__(self, other):
        """Implements inequality check between sequence collections"""

        # __eq__ does not cover use of != operator, so this explicitly has to be implemented
        # Two sequence collections are identical if (1) they have same size, and (2) every sequence
        # in collection has a match in collection 2. Note: sequence names can be different
        if len(self) != len(other):
            return True

        else:
            # For each seq in 1, compare to seqs in 2.
            # If match found: move on to next seq in 1. If no match: return True
            for seq1 in self:
                for seq2 in other:
                    if seq1 == seq2:
                        break
                return True

        # If we fell off loop: collections must be identical. Return False
        return False

    #######################################################################################

    def __str__(self):
        """String (fasta) representation of sequence collection"""

        seqstrings = []
        for seq in self:
            seqstrings.append(str(seq))
        return "\n".join(seqstrings)

    #######################################################################################

    def sortnames(self, reverse=False):
        """Sorts sequences in object alphabetically by name"""
        self.seqnamelist.sort(reverse=reverse)

    #######################################################################################

    def addseq(self, seq, silently_discard_dup_name=False):
        """Add Sequence object to set"""

        # If name already in set of sequences: silently omit, or raise exception depending on options
        name = seq.name
        if name in self.seqdict:
            if not silently_discard_dup_name:
                raise SeqError("Duplicate sequence names: %s" % name)
        else:
            # Set seqtype, alphabet, and ambiguitysymbols of Seq_set if not already set
            # Check type consistency otherwise.
            if not self.seqtype:
                self.seqtype = seq.seqtype
                self.alphabet = seq.alphabet
                self.ambigsymbols = seq.ambigsymbols

            elif seq.seqtype != self.seqtype:
                raise SeqError("Mismatch between sequence types: %s vs. %s" % (self.seqtype, seq.seqtype))

            # Add seq to Seq_set.
            self.seqnamelist.append(seq.name)
            self.seqdict[seq.name] = seq

    #######################################################################################

    def addseqset(self, other, silently_discard_dup_name=False):
        """Adds all sequences in 'other' to collection in 'self'"""
        for seq in other:
            self.addseq(seq, silently_discard_dup_name)
        return self

    #######################################################################################

    def remseq(self, name):
        """Remove Sequence object from set"""

        if name in self.seqdict:
            self.seqnamelist.remove(name)
            del self.seqdict[name]
        else:
            raise SeqError("No such sequence: %s" % name)

    #######################################################################################

    def remseqs(self, namelist):
        """Remove named sequences from set"""

        for name in namelist:
            self.remseq(name)

    #######################################################################################

    def changeseqname(self, oldname, newname, fix_dupnames=False):
        """Change the name of one sequence object from oldname to newname"""

        if newname in self.seqdict and fix_dupnames:
            i = 2
            fixname = newname + "_" + str(i)
            while fixname in self.seqdict:
                i += 1
                fixname = newname + "_" + str(i)
            newname = fixname

        if oldname in self.seqdict:
            if newname != oldname:
                seq = self.seqdict[oldname]
                seq.name = newname
                self.seqdict[newname] = seq
                del self.seqdict[oldname]
                self.seqnamelist.remove(oldname)
                self.seqnamelist.append(newname)
        # else:
        #     raise SeqError("No such sequence: %s" % oldname)
    #######################################################################################

    def getseq(self, name):
        """Returns Sequence object for named sequence"""

        if name in self.seqdict:
            return self.seqdict[name]
        else:
            raise SeqError("No such sequence: %s" % name)

    #######################################################################################

    def subset(self, namelist):
        """Returns new sequence collection object (of proper type) containing named sequences only"""

        # Sanity check: are requested sequences in sequence collection?
        extranames = set(namelist) - set(self.seqnamelist)
        if extranames:
            raise SeqError("Requested subset contains names that are not in sequence collection: %s" % extranames)
        else:
            subset = self.__class__()       # Create new object of same class as self (don't know which subclass we're in)
            for name in self.seqnamelist:   # Done this way to retain ordering from original sequence collection ("namelist" has random order)
                if name in namelist:
                    subset.addseq(self.getseq(name))     # Should I copy the seq objects before addseq'ing?

        return subset

    #######################################################################################

    def subsample(self, samplesize):
        """Returns new sequence collection object (of proper type) containing samplesize randomly selected sequences"""

        # Sanity check: is samplesize larger than total size of sequence set? Is it less than zero?
        if samplesize > len(self):
            raise SeqError("Requested samplesize larger than full data set")
        if samplesize < 0:
            raise SeqError("Requested samplesize is negative - must be positive integer")
        else:
            subsample = self.__class__()    # Create new object of same class as self (don't know which subclass we're in)
            for seq in random.sample(list(self), samplesize):   # Note: indexing self refers to sequence objects
                subsample.addseq(seq)

        return subsample

    #######################################################################################

    def subseq(self, start, stop, slicesyntax=True, rename=True, aln_name=None, aln_name_number=False):
        """Returns specified positions as Seq_alignment or Seq_set object (depending on class of self)"""

        # If slicesyntax is False: indexing starts at 1, and stop is included
        # If rename or renamealn is True: start and stop indices added to seqnames or alignmentname
        if aln_name is None:
            if aln_name_number:
                aln_name = "{}_{}_{}".format(self.name, start, stop)
            else:
                aln_name = self.name
        subseqs = self.__class__(name = aln_name)       # Create new object of same class as self (don't know which subclass we're in)

        for seq in self:
            newseq = seq.subseq(start, stop, slicesyntax, rename)
            subseqs.addseq(newseq)

        return subseqs

    #######################################################################################

    def getnames(self):
        """Returns list of names of all sequences in set"""

        return copy.copy(self.seqnamelist)

    #######################################################################################

    def range(self, rangefrom, rangeto):
        """Discards all symbols outside of range from sequences in set"""

        # Sanity check: rangeto > rangefrom ?
        if rangefrom > rangeto:
            raise SeqError("End-of-range index is higher than start-of-range index")

        for seq in self:
            if rangeto > len(seq):
                raise SeqError("Range exceeds length of sequence %s: %d" % (seq.name, len(seq)))
            else:
                seq = seq[rangefrom:rangeto]

    #######################################################################################

    def removedupseqs(self):

        # Removes all duplicate sequences (keeping one of each type).
        # Returns lists of identical sequences (first one in list is the one that is kept)

        # Keep track of sets of identical sequences using list and dict:
        # simlist: list containing lists of identical sequence names
        # indexdict: keys are seqnames, value is index in simlist where name occurs
        simlist = []
        indexdict = {}
        for seq1, seq2 in itertools.combinations(self, 2):
            if seq1 == seq2:
                if seq1.name in indexdict and seq2.name not in indexdict:
                    listno = indexdict[seq1.name]
                    simlist[listno].append(seq2.name)
                    indexdict[seq2.name] = listno
                elif seq2.name in indexdict and seq1.name not in indexdict:
                    listno = indexdict[seq2.name]
                    simlist[listno].append(seq1.name)
                    indexdict[seq1.name] = listno
                elif seq1.name not in indexdict and seq2.name not in indexdict:
                    simlist.append([seq1.name, seq2.name])
                    newindex = len(simlist) - 1
                    indexdict[seq1.name] = newindex
                    indexdict[seq2.name] = newindex

        # Remove duplicates (i.e, last seqs in all simlist lists):
        for namelist in simlist:
            for dupname in namelist[1:]:
                self.remseq(dupname)

        # Return list of lists of duplicates (first one in each list was kept)
        return simlist

    #######################################################################################

    def residuecounts(self):
        """Returns dictionary with cumulated counts for set of seqs {symbol:count}"""

        allcounts = Counter()
        for seqobj in self:
            allcounts.update(seqobj.seq)   # Perhaps a bit unclean to rely on implementation
        return allcounts

    #######################################################################################

    def composition(self, ignoregaps=True):
        """Returns dictionary with cumulated counts AND freq for set of seqs {symbol:[count,freq]}"""

        allcounts = self.residuecounts()
        if ignoregaps:
            alphabet = set(allcounts.keys()) - set("-")
            totlength = sum(allcounts.values()) - allcounts["-"]
        else:
            alphabet = set(allcounts.keys())
            totlength = sum(allcounts.values())

        compdict = {}
        for residue in alphabet:
            count = allcounts[residue]
            freq = count/totlength
            compdict[residue] = [count, freq]
        return compdict

    #######################################################################################

    def clean_names(self, illegal=",:;()[]", rep="_"):
        """Change all sequence names such that chars in 'illegal' are replaced by 'rep'.
        For instance before creating phylogenetic tree from sequences"""

        illegal_esc_list = [re.escape(char) for char in illegal]
        illegal_esc_list.append("_") # To avoid two underscores in a row
        illegal_esc_string = "".join(illegal_esc_list)
        regex = f"[{illegal_esc_string}]+"
        for old in self.getnames():
            new = re.sub(regex,"_",old)
            self.changeseqname(old, new)

    #######################################################################################

    def rename_numbered(self, basename, namefile=None):
        """Renames all sequences in collection to this form: basename_001, ..."""

        if namefile:
            transfile = open(namefile, "w")
        nfill = int(log(len(self)) / log(10)) + 1     # Number of digits required
        for i, oldname in enumerate(self.getnames()):
            newname = basename + "_" + str(i + 1).rjust(nfill, "0")       # Pad running index with zeroes
            if namefile:
                transfile.write("%s\t%s\n" % (newname, oldname))
            self.changeseqname(oldname, newname)

    #######################################################################################

    def rename_regexp(self, pattern, namefile=None, silently_discard_dup_name=False, fix_dupnames=False):
        """Renames all sequences in collection by deleting parts of name matching regular expression pattern string"""

        # Note: renaming this way may lead to duplicate names. Need to check this during execution
        if namefile:
            transfile = open(namefile, "w")
        orignames = self.getnames()
        for oldname in orignames:
            newname = re.sub(pattern, "", oldname)
            # Note: output from self.getnames will change during execution as names are changed
            if oldname != newname and newname in self.getnames():
                if fix_dupnames:
                    pass
                elif silently_discard_dup_name:
                    self.remseq(oldname)
                else:
                    raise SeqError("Name clash during renaming by regular expression: two sequences will get name '%s'" % newname)
            if namefile:
                transfile.write("%s\t%s\n" % (newname, oldname))
            if (oldname != newname):
                self.changeseqname(oldname, newname, fix_dupnames)

    #######################################################################################

    def transname(self, namefile):
        """Translate all names using oldname/newname pairs in namefile. E.g., use to restore names after renaming"""

        transfile = open(namefile)
        transdict = {}
        for line in transfile:
            words = line.split()
            oldname = words[0]
            newname = words[1]
            if oldname not in self.seqdict:
                raise SeqError("No sequence with this name: %s" % (oldname))
            else:
                transdict[oldname] = newname

        # Uses getnames() (which returns copy of namelist) to avoid iterating over list and dict while changing them
        for oldname in self.getnames():
            newname = transdict[oldname]
            self.changeseqname(oldname, newname)

    #######################################################################################

    def revcomp(self):
        """Returns new sequence collection (of same type as original) reverse complemented"""

        revcompseqs = self.__class__()    # Create new object of same class as self (don't know which subclass we're in)
        for seq in self:
            revcompseqs.addseq(seq.revcomp())     # Should I copy the seq objects before addseq'ing?
        return revcompseqs

    #######################################################################################

    def translate(self):
        """If seqtype is DNA: Returns new sequence collection with seqtype protein"""

        if not self.seqtype.lower() == "dna":
            raise SeqError("Attempt to translate sequences of wrong type: {}".format(self.seqtype))

        protseqs = self.__class__(seqtype="protein")    # Create new object of same class as self (don't know which subclass we're in)
        for seq in self:
            protseqs.addseq(seq.translate())     # Should I copy the seq objects before addseq'ing?
        return protseqs

    #######################################################################################

    def fasta(self, width=60, nocomments=False):
        """Returns fasta-formatted set of sequences as a string"""

        if len(self) == 0:
            raise SeqError("No sequences in sequence set.  Can't create fasta")

        fastalist = []
        for seq in self:
            fastalist.append(seq.fasta(width, nocomments))
            fastalist.append("\n")

        # Remove last newline
        del fastalist[-1]

        return "".join(fastalist)

    #######################################################################################

    def how(self, width=60, nocomments=False):
        """Returns HOW-formatted set of sequences as a string"""

        if len(self) == 0:
            raise SeqError("No sequences in sequence set.  Can't create HOW")

        howlist = []
        for seq in self:
            howlist.append(seq.how(width, nocomments))
            howlist.append("\n")

        # Remove last newline
        del howlist[-1]

        return "".join(howlist)

    #######################################################################################

    def tab(self, nocomments=False):
        """Returns tab-formatted set of sequences as a string"""

        if len(self) == 0:
            raise SeqError("No sequences in sequence set.  Can't create TAB")

        tablist = []
        for seq in self:
            tablist.append(seq.tab(nocomments))
            tablist.append("\n")

        # Remove last newline
        del tablist[-1]

        return "".join(tablist)

    #######################################################################################

    def raw(self):
        """Returns raw-formatted set of sequences as a string"""

        if len(self) == 0:
            raise SeqError("No sequences in sequence set.  Can't create RAW")

        rawlist = []
        for seq in self:
            rawlist.append(seq.raw())
            rawlist.append("\n")

        # Remove last newline
        del rawlist[-1]

        return "".join(rawlist)

#############################################################################################
#############################################################################################

class Seq_set(Sequences_base):
    """Class containing a set of sequences"""

    # Note: most functionality is in abstract baseclass "Sequences_base"
    # This derived class only contains methods that should never be used on aligned sequences


    #######################################################################################

    def __init__(self, name=None, seqtype=None, seqlist=None):
        Sequences_base.__init__(self, name=name, seqtype=None)
        self.alignment = False
        self.seqpos2alignpos_cache = {}
        self.alignpos2seqpos_cache = {}

        if seqlist:
            for seq in seqlist:
                self.addseq(seq)

    #######################################################################################

    def remgaps(self):
        """Removes gap characters from all sequences in Seq_set"""

        for seq in self:
            seq.remgaps()

#############################################################################################
#############################################################################################

class Seq_alignment(Sequences_base):
    """Class representing aligned set of sequences"""

    #######################################################################################

    def __init__(self, name=None, seqtype=None):
        Sequences_base.__init__(self, name=name, seqtype=None)
        self.alignment = True
        self.seqpos2alignpos_cache = {}
        self.alignpos2seqpos_cache = {}
        self.annotation = None

    #######################################################################################

    def addseq(self, seq, silently_discard_dup_name=False):
        """Add Sequence object to alignment"""
        # Overrides baseclass function. Sets partition info + ensures consistent sequencelengths
        if len(self) == 0:
            self.partitions = [(self.name, 0, len(seq))]  # List of (name, partition-start, partition-length) tuples
        else:
            if len(seq) != len(self[0]):
                raise SeqError("Not an alignment: these sequences have different lengths: %s and %s" % (seq.name, self[0].name))

        # If length was OK, add sequence by calling baseclass method
        Sequences_base.addseq(self, seq, silently_discard_dup_name)

    #######################################################################################

    def appendalignment(self, other):
        """Appends sequences in 'other' to the end of similarly named sequences in 'self'"""

        # Update partition info in self
        partitiontuple = (other.name, self.alignlen(), other.alignlen())    # (name, partition-start, partition-length)
        if len(self) == 0:
            self.partitions = [partitiontuple]  # List of (name, partition-start, partition-length) tuples
        else:
            self.partitions.append(partitiontuple)

        # If alignment object is empty: Initialise with other
        if len(self) == 0:
            self.addseqset(other)
        # Else: match sequences, appending end to end
        else:
            for seq in self:
                try:
                    matchingseq = other.getseq(seq.name)
                except SeqError:
                    # Re-throw exception with more precise description of problem (we know more than just that name was not found)
                    raise SeqError("Sequences in files have different names. No match found for %s" % seq.name)
                seq.appendseq(matchingseq)

    #######################################################################################

    def alignlen(self):
        """Returns length of alignment or zero if no sequences"""
        if len(self) == 0:
            return(0)
        else:
            return len(self[0])

    #######################################################################################

    def getcolumn(self, i):
        """Returns column 'i' as list"""

        # Implementation note: code would be much clearer, and more robust, if I used
        # the indexing interface from seq and seqset, but directly accessing attributes
        # like this is considerably faster, and this function needs to be efficient!
        # Also: for python>=3.6 order in dict is kept
        return [seqobj.seq[i] for seqobj in self.seqdict.values()]

    #######################################################################################

    def columns(self):
        """Returns column iterator object"""

        # Function that returns column iterator object, which can be used to iterate over columns
        # (Note that iterating over Seq_alignment object itself returns one _sequence_ at a time - not column).
        # Main trick here is a class which keeps track of its enclosing alignment object ("parent"),
        # while at the same time keeping track of which column number the iteration has gotten to

        class Col_iterator:
            """Column iterator object"""

            def __init__(self, parent):
                self.i=0
                self.length = len(parent[0])
                self.parent = parent

            def __iter__(self):
                return self

            def __next__(self):
                if self.i<self.length:
                    self.i += 1
                    return self.parent.getcolumn(self.i-1)
                else:
                    raise StopIteration

        return Col_iterator(self)

    #######################################################################################

    def conscols(self):
        """Returns list of columns that are conserved"""
        conscols = []
        for i, col in enumerate(self.columns()):
            if len(set(col)) == 1:        # nvalues == 1 <=> conserved column
                conscols.append(i)
        return conscols

    #######################################################################################

    def varcols(self):
        """Returns list of columns that are variable (not conserved)"""
        varcols = []
        for i, col in enumerate(self.columns()):
            if len(set(col)) > 1:        # nvalues == 1 <=> conserved column
                varcols.append(i)
        return varcols

    #######################################################################################

    def gappycols(self):
        """Returns list of columns that are gappy (contain one or more gaps)"""
        gappycols = []
        for i, col in enumerate(self.columns()):
            if "-" in col:
                gappycols.append(i)
        return gappycols

    #######################################################################################

    def site_summary(self):
        """Returns tuple with summary of sites: (alignlen, n_var_cols, n_gappy_cols)"""
        varcols = self.varcols()
        gappycols = self.gappycols()
        length = self.alignlen()
        return (length, len(varcols), len(gappycols))

    #######################################################################################

    def indexfilter(self, keeplist):
        """Discards all columns whose indices are not in keeplist"""
        for seq in self:
            seq.indexfilter(keeplist)

        # If alignment object contains annotation: also apply filter to annotation sequence
        if self.annotation:
            s = self.annotation
            annotlist = [s[i] for i in keeplist]
            self.annotation = "".join(annotlist)

    #######################################################################################

    def remcols(self, discardlist):
        """Discards all columns whose indices are in discardlist"""
        keeplist = []
        for i in range(self.alignlen()):
            if i not in discardlist:
                keeplist.append(i)
        self.indexfilter(keeplist)

    #######################################################################################

    def remambigcol(self):
        """Removes all columns that contain one or more ambiguity symbols"""

        # Construct list of columns with no ambiguity symbols (these will be kept)
        keeplist = []
        for i in range(self.alignlen()):
            columnset = set(self.getcolumn(i))
            if not (self.ambigsymbols & columnset):
                keeplist.append(i)

        # Apply filter to sequences so only non-ambiguous columns are kept
        self.indexfilter(keeplist)

    #######################################################################################

    def remgapcol(self):
        """Removes all columns that contain one or more gaps"""

        # Construct list of columns with no gaps (these will be kept)
        keeplist = []
        for i in range(self.alignlen()):
            if "-" not in self.getcolumn(i):
                keeplist.append(i)

        # Apply filter to sequences so only non-gapped columns are kept
        self.indexfilter(keeplist)

    #######################################################################################

    def remallgapcol(self):
       """Removes columns that contain only gaps"""

       # Construct keeplist: columns that (1) are not conserved, or (2) contain a non-gap symbol
       keeplist = []
       for i in range(self.alignlen()):
           columnchars = set(self.getcolumn(i))     # set() automatically uniquefies characters
           if len(columnchars) > 1 or "-" not in columnchars:
               keeplist.append(i)

       # Apply filter to sequences so gap-only columns are discarded
       self.indexfilter(keeplist)

    #######################################################################################

    def remfracgapcol(self, frac):
        """Removes all columns where more than "frac" fraction of residues are gaps"""

        # Construct list of columns with <= frac fraction gaps (these will be kept)
        keeplist = []
        nseqs = len(self)
        for i in range(self.alignlen()):
            col = self.getcolumn(i)
            gapfrac = col.count("-") / nseqs
            if gapfrac <= frac:
                keeplist.append(i)

        # Apply filter to sequences so only columns with <= frac fraction gaps are kept
        self.indexfilter(keeplist)

    #######################################################################################

    def remconscol(self):
        """Removes columns where all symbols are identical (conserved columns)"""
        conscols = self.conscols()
        self.remcols(conscols)

    #######################################################################################

    def rem_hmmalign_insertcol(self):
        """For alignments created by HMMer:hmmalign: Removes insert state columns"""

        # In alignments from HMMer's hmmalign: insert states have lower case residues
        # and/or "." for gaps (not really gaps - just lack of insert in other seq).
        # Internally I keep track of these using Seq_alignment.annotation field:
        # "i" = insert state, "m" = match state
        discardlist = []
        if self.annotation:
            for i,char in enumerate(self.annotation):
                if char == "i":
                    discardlist.append(i)
        else:
            raise SeqError("This alignment contains no information about hmmalign insert states")

        self.remcols(discardlist)

    #######################################################################################

    def align_seq_pos_cache_builder(self, seqname):
        """Helper function for seqpos2alignpos and alignpos2seqpos: Builds mapping caches for given seq"""

        self.seqpos2alignpos_cache[seqname] = []
        self.alignpos2seqpos_cache[seqname] = {}
        seq = self.getseq(seqname)
        seq_i = 0

        for align_i in range(self.alignlen()):
            if seq[align_i] != "-":
                gapstatus = False
                self.seqpos2alignpos_cache[seqname].append(align_i)
                self.alignpos2seqpos_cache[seqname][align_i] = (seq_i, gapstatus)
                seq_i += 1
            else:
                # If there is a gap at this position: Return seqpos for preceding residue in sequence
                gapstatus = True
                self.alignpos2seqpos_cache[seqname][align_i] = (seq_i - 1, gapstatus)

    #######################################################################################

    def seqpos2alignpos(self, seqname, seqpos, slicesyntax=True):
        """Input: seqname and position in sequence. Output: alignment position corresponding to seqpos.
            Slicesyntax = True: All numbering starts at 0. Slicesyntax = False: All numbering starts at 1"""
        # NOTE I am not handling the situation where cache is present but user requests info for non-existing site? ALso below?

        # Construct caches with all info for seqname, if it is not in caches already
        if seqname not in self.seqpos2alignpos_cache:
            self.align_seq_pos_cache_builder(seqname)

        # Return requested position.
        # Caches are in slicesyntax (numbering starts at 0), so query and result must be corrected if non-slicesyntax requested
        if slicesyntax:
            return self.seqpos2alignpos_cache[seqname][seqpos]
        else:
            return self.seqpos2alignpos_cache[seqname][seqpos - 1] + 1

    #######################################################################################

    def alignpos2seqpos(self, seqname, alignpos, slicesyntax=True):
        """Input: seqname and position in alignment. Output: tuple of (seqpos, GAPSTATUS)
            If alignpos corresponds to a gap, GAPSTATUS will be True, and seqpos is pos of preceding residue in seq.
            Slicesyntax = True: Numbering starts at 0. Slicesyntax = False: Numbering starts at 1"""

        # Construct caches with all info for seqname, if it is not in caches already
        if seqname not in self.alignpos2seqpos_cache:
            self.align_seq_pos_cache_builder(seqname)

        # Return requested position.
        # Caches are in slicesyntax (numbering starts at 0), so query and result must be corrected if non-slicesyntax requested
        if slicesyntax:
            return self.alignpos2seqpos_cache[seqname][alignpos]
        else:
            restup = self.alignpos2seqpos_cache[seqname][alignpos - 1]
            return (restup[0] + 1, restup[1])

    #######################################################################################

    def shannon(self, countgaps=True):
        """Computes Shannon entropy for all columns in the alignment, returns list of values (index=position in seq)"""

        # NOTE: if countgaps==True, then gaps are treated as an extra symbol. Otherwise they are ignored
        shannonlist = []
        numseqs = len(self)
        for col in self.columns():
            symbolcounts = Counter(col)
            if not countgaps:
                numseqs = len(self) - col.count("-")
                del symbolcounts["-"]       # NOTE: what if column is all gaps???
            entropy = 0.0
            for symbol in symbolcounts:
                symbolfreq = float(symbolcounts[symbol]) / numseqs
                entropy += symbolfreq * log(symbolfreq, 2)
            shannonlist.append(- entropy)

        return shannonlist

    #######################################################################################

    def nucfreq(self, poslist=None):
        """Computes nucleotide frequency distribution for requested columns in the alignment (default: all).
        Returns matrix of values: Rows = A, C, G, T. Column=position in seq."""

        # NOTE: gaps and IUPAC symbols are ignored
        # NOTE 2: this could be a protein alignment - I need to make sequencetype an issue here
        #         make dna and protein versions and set up factory method that selects proper class on construction
        #         (needs some sequences or the type then)
        if poslist:
            nrow = len(poslist)
            positions = poslist
        else:
            nrow = self.alignlen()
            positions = range(nrow)

        freqmat = np.zeros(shape = (nrow, 4))

        for i,pos in enumerate(positions):
            col = self.getcolumn(pos)
            counts = Counter(col)
            nA = counts["A"]
            nC = counts["C"]
            nG = counts["G"]
            nT = counts["T"]
            nseq = nA + nC + nG + nT
            freqmat[i, 0] = nA / nseq
            freqmat[i, 1] = nC / nseq
            freqmat[i, 2] = nG / nseq
            freqmat[i, 3] = nT / nseq

        return(freqmat)

    #######################################################################################

    def consensus(self):
        """Returns consensus sequence (= most frequent symbol in each alignment column) as sequence object"""

        # NOTE: all symbols are counted (also gaps and IUPAC). Should this be changed?
        seqlist = []
        for i in range(self.alignlen()):
            col = self.getcolumn(i)
            nuc_counter = Counter(col)
            [(nuc, count)] = nuc_counter.most_common(1) # Method returns list of one tuple, which is unpacked
            seqlist.append(nuc)
        seq_string = "".join(seqlist)

        seq_class = self[0].__class__   # Hack to find class (DNA, Protein, ...) of Sequence objects in this alignment
        con_seq = seq_class(name=self.name, seq=seq_string)

        return con_seq

    #######################################################################################

    def seqset(self):
        """Removes gaps from Seq_alignment object and returns as Seq_set object"""
        seqset = Seq_set(self.name)
        for seq in self:
            newseq = copy.copy(seq)
            newseq.remgaps()
            seqset.addseq(newseq)

        return seqset

    # #######################################################################################
    #
    # def subseq(self, start, stop, slicesyntax=True, rename=False):
    #     """Returns specified columns from alignment as new Seq_alignment object"""
    #
    #     # If slicesyntax is False: indexing starts at 1, and stop is included
    #     # If rename is True: start and stop indices added to name
    #     subseqs = Seq_alignment(self.name)
    #     for seq in self:
    #         newseq = seq.subseq(start, stop, slicesyntax, rename)
    #         subseqs.addseq(newseq)
    #
    #     return subseqs
    #
    #######################################################################################

    def partitions_as_seqalignments(self):
        """Returns list containing all partitions as Seq_alignment objects"""

        plist = []
        for i in range(len(self.partitions)):
            name = self.partitions[i][0]
            if name == "Partition" or name is None:        # Name has not been set, use number instead
                name = "partition_{:0>2d}".format(i+1)
            pstart = self.partitions[i][1]
            pstop = pstart + self.partitions[i][2]
            part_alignment = self.subseq(pstart, pstop, rename=True, aln_name=name)
            plist.append(part_alignment)

        return plist

    #######################################################################################

    def distdict(self, dist="pdist"):
        """Returns nested dict of pairwise sequence distances: distdict[name1][name2] = dist"""

        # Decide which method to use for computing distance
        methoddict = {"hamming":Sequence.hamming,
                      "pdist":Sequence.pdist,
                      "hamming_ignoregaps":Sequence.hamming_ignoregaps,
                      "pdist_ignoregaps":Sequence.pdist_ignoregaps}
        try:
            distmethod = methoddict[dist]
        except KeyError:
            raise SeqError("Unknown distance measure: %s" % dist)

        # Initialize empty nested 2D dictionary with sequence names as keys.
        distdict = dict.fromkeys(self.seqnamelist)
        for name in self.seqnamelist:
            distdict[name] = dict.fromkeys(self.seqnamelist)

        # Fill dictionary with values
        for seq1, seq2 in itertools.combinations(self, 2):
            dist = distmethod(seq1,seq2)
            distdict[seq1.name][seq2.name] = dist
            distdict[seq2.name][seq1.name] = dist

        return distdict

    #######################################################################################

    def sequence_diversity(self, ignoregaps=False):
        """
        Compute pairwise sequence diversity (pi).
        Return mean, and standard deviation, as tuple: (mean, std)
        Discard pairwise gappy positions if requested.
        """
        # Note: std should perhaps be computed according to Nei, Molecular Evolutionary Genetics, equation 10.7?
        # Here std is computed using Welfords 1962 single-pass algorithm
        # Correlations among pairwise distances are ignored
        mean_prev = var_prev = num_vals = 0.0
        nseqs = len(self)

        if ignoregaps:
            distmethod = Sequence.pdist_ignoregaps
        else:
            distmethod = Sequence.pdist

        # Online (single-pass) computation of mean and variance
        # Also keep track of max and min
        minpi = math.inf
        maxpi = -math.inf
        for s1, s2 in itertools.combinations(self, 2):
            num_vals += 1.0
            dist = distmethod(s1, s2)
            if dist < minpi:
                minpi = dist
            if dist > maxpi:
                maxpi = dist
            diff = dist - mean_prev
            mean_cur = mean_prev + diff / num_vals
            var_cur = var_prev + diff * (dist - mean_cur)
            mean_prev = mean_cur
            var_prev = var_cur

        # NOTE: need to doublecheck this:
        variance = var_cur / (num_vals)
        std = sqrt(variance)
        return (mean_cur, std, minpi, maxpi)

    #######################################################################################

    def overlap(self, other, discardgaps=True):
        """Computes overlap between two alignments as fraction of identically paired residues"""

        # Implementation note: measure is based on representing each pairwise alignment (implied
        # by the multiple alignment) as a "set of index_tuples" of this form: set((1,1), (2,2), (3,4), ...)
        # Here each residue is represented by its index in the sequence (not counting gaps).
        # Pairs where one or both symbols are gaps are discarded, unless discardgaps=False.
        # Thus the following alignment:
        #
        #  AC---GT
        #  A-T--GT
        #
        # Would be: set((1,1), (3,3), (4,4)), or: set((1,1), (2,"-"), ("-",2), ("-","-"), (3,3), (4,4))
        # Note that set() automatically removes duplicate entries (hence only one ("-","-") tuple).
        # Overlap for this seqpair is then found by first finding the set of tuples representation of the same
        # seqpair in the second alignment, and then using set arithmetic to find the union of the
        # set of tuples for the two alignments. overlapfrac = size of the union divided by avg set length
        #
        # Note 2: it might be relevant to represent gaps by "_residueindex" to differentiate gaps in different places

        ###################################################################################

        def seq2index(seq):

            # Converts sequence to indexformat:
            # residues are represented by residue number, gaps are represented by "-"
            index = 0
            indexlist = []
            for char in seq:
                if char == "-":
                    indexlist.append("-")
                else:
                    index += 1
                    indexlist.append(index)
            return indexlist

        ###################################################################################

        if set(self.getnames()) != set(other.getnames()):
            raise SeqError("Alignments do not contain same sequences - not possible to compute overlap")
        numseqs = len(self)
        overlapfrac = 0.0

        # Convert each sequence in each alignment to indexformat
        align1dict = {}
        align2dict = {}
        for seq in self:
            align1dict[seq.name] = seq2index(seq.seq)
        for seq in other:
            align2dict[seq.name] = seq2index(seq.seq)

        # For each pair of sequences in each alignment:
        # construct list of tuples representation of mapping, compute overlap by set arithmetic
        for s1, s2 in itertools.combinations(self, 2):
            mapping1 = list(zip(align1dict[s1.name], align1dict[s2.name]))
            mapping2 = list(zip(align2dict[s1.name], align2dict[s2.name]))

            if discardgaps:
                mapping1 = [tup for tup in mapping1 if "-" not in tup]
                mapping2 = [tup for tup in mapping2 if "-" not in tup]

            avglen = (len(mapping1) + len(mapping2)) / 2.0
            overlapfrac += len(set(mapping1) & set(mapping2)) / avglen

        overlapfrac = overlapfrac / ((numseqs * (numseqs - 1)) / 2)
        return overlapfrac

    #######################################################################################

    def phylip(self, width=60):
        """Returns phylip-formatted, interleaved alignment as string"""

        if len(self) == 0:
            raise SeqError("No sequences in sequence set.  Can't create phylip")

        alignlen = self.alignlen()
        numseq = len(self)
        numseqblocks = int(ceil(self.alignlen()/float(width)))

        # Find longest name (NOTE: could also check this during construction of alignment)
        lenlist = [len(seq.name) for seq in self]
        frontspace = max(lenlist) + 2

        # Header
        phylip = ["%6d   %d\n" % (numseq, alignlen)]

        # First block, with names
        for seq in self:
            phylip.append(seq.name.ljust(frontspace))
            phylip.append(seq[0:width])
            phylip.append("\n")
        phylip.append("\n")

        # Remaining sequence blocks
        for i in range(1,numseqblocks):
            for seq in self:
                phylip.append(" ".rjust(frontspace))
                phylip.append(seq[i*width:i*width+width])
                phylip.append("\n")
            phylip.append("\n")
        del phylip[-2:] # added two newlines too many

        return "".join(phylip)

    #######################################################################################

    def clustal(self, width=60):
        """Returns clustal-formatted alignment as string"""

        if len(self) == 0:
            raise SeqError("No sequences in sequence set.  Can't create clustal")

        alignlen = self.alignlen()
        numseqblocks = int(ceil(alignlen/float(width)))

        # If format is not protein (presumably "DNA") then conservation line is limited to scoring conserved
        # versus non-conserved. If format is "protein" then conservation line should score
        # strong and weak conservation also.

        # Format is "protein"
        if self.seqtype == "protein":

            # Predefined conserved and semiconserved amino acid sets for constructing conservation line
            # Information is kept in form of list of conserved sets
            # Info taken from http://bioportal.cgb.indiana.edu/docs/tools/clustalw/clustalw.html

            cons_sub_strings = ["STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW"]
            semicons_sub_strings = ["CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK", "NEQHRK", "FVLIM", "FVLIM"]

            cons_sub = []
            for string in cons_sub_strings:
                cons_sub.append(set(string))

            semicons_sub = []
            for string in semicons_sub_strings:
                semicons_sub.append(set(string))

            # Function for checking whether "query_set" is a subset of one of the symbol sets in "set_list"
            def containsset(query_set, set_list):
                for symbol_set in set_list:
                    if query_set <= symbol_set:
                        return True

                # If we fell off for loop then set is not contained
                return False

            # Compute conservation line
            conservation = []
            for i in range(alignlen):
                col = set(self.getcolumn(i))    # set() automatically uniquefies list of symbols
                if len(col) ==1:
                    conservation.append("*")    # Only one symbol: full conservation = "*"
                elif containsset(col, cons_sub):
                    conservation.append(":")    # Conserved substitutions = ":"
                elif containsset(col, semicons_sub):
                    conservation.append(".")    # Semiconserved substitutions = "."
                else:
                    conservation.append(" ")    # High diversity = " "

        # Format is "DNA" or "other"
        else:
            # Compute conservation line
            # Note: how about ambiguity symbols? Currently "W" is counted as different from both "A" and "T", etc.
            conservation = []
            for i in range(alignlen):
                col = set(self.getcolumn(i))    # set() automatically uniquefies list of symbols
                if len(col) ==1:
                    conservation.append("*")    # Only one symbol: full conservation = "*"
                else:
                    conservation.append(" ")    # Non conserved column = " "

        # Convert to string
        conservation = "".join(conservation)

        # Find longest name (NOTE: could also check this during construction of alignment)
        lenlist = [len(seq.name) for seq in self]
        frontspace = max(lenlist) + 6

        # Header
        clustal = ["CLUSTAL W (1.82) multiple sequence alignment\n\n\n"]

        # Add blocks
        for i in range(numseqblocks):

            # Sequence block, with names
            for seq in self:
                clustal.append(seq.name.ljust(frontspace))
                clustal.append(seq[i*width:i*width+width])
                clustal.append("\n")

            # Conservation line
            clustal.append("".ljust(frontspace))
            clustal.append(conservation[i*width:i*width+width])
            clustal.append("\n\n")

        # Added one newline too many
        del clustal[-1:]

        return "".join(clustal)

    #######################################################################################

    def nexus(self, width=60, print_partitioned=False):
        """Returns nexus-formatted alignment as string"""

        if len(self) == 0:
            raise SeqError("No sequences in sequence set.  Can't create nexus")

        def add_partition(partitionstart, partitionlen):
            numseqblocks = int(ceil(partitionlen/float(width)))
            for i in range(numseqblocks):
                for seq in self:
                    blockstart = partitionstart + (i * width)
                    blockstop = partitionstart + (i * width) + width
                    if blockstop > partitionstart + partitionlen:
                        blockstop = partitionstart + partitionlen
                    nexus.append("    " + seq.name.ljust(frontspace))
                    nexus.append(seq[blockstart:blockstop])
                    nexus.append("\n")
                nexus.append("\n")
            del nexus[-1:]              # Added one newline too many

        ####################################################################################

        alignlen = self.alignlen()
        numseq = len(self)
        seqtype = self.seqtype.lower()
        alphabet = "".join(sorted(list(self.alphabet)))

        # Find longest name (NOTE: could also check this during construction of alignment)
        lenlist = [len(seq.name) for seq in self]
        frontspace = max(lenlist) + 6
        nexus = []

        # Header
        nexus.append("#NEXUS\n\n")
        nexus.append("begin data;\n")
        nexus.append("    dimensions ntax=%d nchar=%d;\n" % (numseq, alignlen))
        nexus.append("    format datatype=%s interleave=yes gap=-;\n\n" % (seqtype))
        nexus.append("    matrix\n")

        # Sequence blocks. Divides nexus file into one subsection per partition if requested
        if print_partitioned:
            for i in range(len(self.partitions)):
                name = self.partitions[i][0]
                if name == "Partition":        # Name has not been set, use number instead
                    name = str(i + 1)
                nexus.append("\n    [Partition: %s]\n" % name)
                add_partition(self.partitions[i][1], self.partitions[i][2])
        else:
            add_partition(0, alignlen)

        # Footer
        nexus.append(";\nend;")

        return "".join(nexus)

    #######################################################################################

    def nexusgap(self, width=60):
        """Returns nexus-formatted alignment with binary gap encoding, mixed mrbayes format, as string"""

        alignlen = self.alignlen()
        numseq = len(self)
        numseqblocks = int(ceil(alignlen/float(width)))
        seqtype = self.seqtype.lower()
        alphabet = "".join(sorted(list(self.alphabet)))

        # Find longest name (NOTE: could also check this during construction of alignment)
        lenlist = [len(seq.name) for seq in self]
        frontspace = max(lenlist) + 6

        # Header
        # Note: alignment is prepended with binary encoding of same length
        nexus = ["#NEXUS\n\n"]
        nexus.append("begin data;\n")
        nexus.append("dimensions ntax=%d nchar=%d;\n" % (numseq, 2 * alignlen))
        nexus.append('format datatype=mixed(restriction:1-%d,%s:%d-%d) interleave=yes gap=-;\n\n' %
                                                 (alignlen, seqtype, alignlen + 1, 2 * alignlen))
        nexus.append("matrix\n")

        # Add binary formatted blocks, with names
        for i in range(numseqblocks):
            for seq in self:
                nexus.append(seq.name.ljust(frontspace))
                nexus.append(seq.gapencoded()[i*width:i*width+width])
                nexus.append("\n")
            nexus.append("\n")

        # Add actual sequence blocks, with names
        for i in range(numseqblocks):
            for seq in self:
                nexus.append(seq.name.ljust(frontspace))
                nexus.append(seq[i*width:i*width+width])
                nexus.append("\n")
            nexus.append("\n")

        # Added one newline too many
        del nexus[-1:]
        nexus.append(";\nend;")

        return "".join(nexus)


    #######################################################################################

    def charsetblock(self):
        """Returns MrBayes block with charset commands for partitioned analysis as string"""

        numpartitions = len(self.partitions)

        # Find longest name and longest number to enable properly aligned output
        namewidth = 0
        for i in range(numpartitions):
            name = self.partitions[i][0]
            namewidth = max(namewidth, len(name))
        largestnum = self.partitions[numpartitions-1][1] + self.partitions[numpartitions-1][2] - 1
        numwidth = floor( log10( largestnum ) + 1 )

        # Header
        block = ["begin mrbayes;\n"]

        # Charset commands
        block.append("    [Charset commands:]\n")
        for i in range(numpartitions):
            name = self.partitions[i][0]
            start = self.partitions[i][1]
            stop = start + self.partitions[i][2] - 1
            block.append("    charset {n:{naw}} = {b:d}-{e:d};\t[partition no. {p}]\n".format(naw=namewidth,
                                                                                                  n=name,
                                                                                                  b=start + 1,
                                                                                                  e=stop + 1,
                                                                                                  p=i + 1)) # +1: not list syntax
        block.append("    partition allgenes = %d: " % numpartitions)
        for i in range(numpartitions):
            name = self.partitions[i][0]
            block.append("%s" % name)
            block.append(", ")
            # print 4 names per line
            if all([i%4==0, i>0, i<numpartitions - 1]):
                block.append("\n        ")
        del block[-1]               # Added one comma too many
        block.append(";\n")
        block.append("    set partition = allgenes;\n")

        # Footer
        block.append("end;\n")
        return "".join(block)


    #######################################################################################

    def mbpartblock(self):
        """Returns MrBayes block for partitioned analysis as string"""

        numpartitions = len(self.partitions)

        # Header
        block = ["begin mrbayes;\n    set autoclose=yes nowarn=yes;\n"]

        block.append("\n\n    [######################################################################################]\n\n")

        # Outgroup command
        block.append("    [Outgroup command:]\n")
        block.append("    [An outgroup may be specified. Replace MYOUTGROUP with actual value]\n\n")
        block.append("    [outgroup MYOUTGROUP;]\n")

        block.append("\n\n    [######################################################################################]\n\n")

        # Charset commands
        block.append("    [Charset commands:]\n")
        block.append("    [The location of each gene is defined by CHARSET. For example, CHARSET gene1 = 1 - 500;\n")
        block.append("    says that the first 500 nucleotides belong to the gene 'gene1'. The command 'partition'\n")
        block.append("    divides sequences into genes specified by CHARSET. The partition must be activated by\n")
        block.append("    the command 'set partition=gene'.]\n\n")
        for i in range(numpartitions):
            name = self.partitions[i][0]
            start = self.partitions[i][1]
            stop = start + self.partitions[i][2] - 1
            block.append("    charset %s = %6d - %6d;\t\t[partition no. %d]\n" % (name, start + 1, stop + 1, i + 1)) # +1: not list syntax
        block.append("    partition allgenes = %d: " % numpartitions)
        for i in range(numpartitions):
            name = self.partitions[i][0]
            block.append("%s" % name)
            block.append(", ")
            # print 10 names per line
            if all([i%10==0, i>0, i<numpartitions - 1]):
                block.append("\n        ")
        del block[-1]               # Added one comma too many
        block.append(";\n")
        block.append("    set partition = allgenes;\n")

        block.append("\n\n    [######################################################################################]\n\n")

        # Model specification commands
        block.append("    [Substitution model commands:]\n")
        block.append("    [Specify substitution models for each partition. Different partitions may all use the \n")
        block.append("    same type of model or different types. Nucleotide models are specified using 'lset', \n")
        block.append("    while amino-acid models are specified using 'prset' and 'lset'.\n\n")
        block.append("    Separate partitions are assumed to share tree topology, but differ in all other model aspects.\n")
        block.append("    The command 'unlink' is used to unlink relevant parameters across partitions.]\n\n")

        if self.seqtype == "DNA":
            block.append("    [DNA commands - alter as needed]\n")
            for i in range(numpartitions):
                block.append("    lset applyto=({}) nst=mixed rates=invgamma;\n".format(i+1))
            block.append("    unlink shape=(all) statefreq=(all) tratio=(all) revmat=(all) pinvar=(all);\n")

        elif self.seqtype == "protein":
            block.append("    [Protein commands - alter as needed]\n")
            for i in range(numpartitions):
                block.append("    prset applyto=({}) aamodelpr=mixed; lset applyto=({}) rates=invgamma;\n".format(i+1, i+1))
            block.append("    unlink shape=(all) pinvar=(all);\n")

        block.append("    prset applyto=(all) ratepr=variable;")

        block.append("\n\n    [######################################################################################]\n\n")

        # MCMC commands
        block.append("    [MrBayes commands:]\n")
        block.append("    [Specify how to run MCMC and possibly sump and sumt as for a regular MrBayes run]\n")
        block.append("    mcmcp ngen=1000000 nchain=3 samplefreq=1000 diagnfreq=50000 printfreq=5000;\n")

        block.append("\n\n    [######################################################################################]\n\n")

        # Add final end statement and return mrbayes block as string
        block.append("end;\n")
        return "".join(block)


    #######################################################################################

    def bestblock(self):
        """Returns MrBayes block for BEST (species tree) analysis as string"""

        numpartitions = len(self.partitions)

        # Header
        bestblock = ["begin mrbayes;\n    set autoclose=yes nowarn=yes;\n"]

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # Outgroup command
        bestblock.append("    [Outgroup command:]\n")
        bestblock.append("    [An outgroup has to be specified. Replace MYOUTGROUP with actual value]\n\n")
        bestblock.append("    outgroup MYOUTGROUP;\n")

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # Taxset commands
        bestblock.append("    [Taxset commands:]\n")
        bestblock.append("    [The command 'TAXSET' tells the program which sequences belong to which species.\n")
        bestblock.append("    The commands below indicate that each species has only one sequence (single allele data).\n")
        bestblock.append("    Although the species names below coincide with the sequences' names, it is totally valid\n")
        bestblock.append("    to use other names for species. For multiple allele data, multiple sequences may belong\n")
        bestblock.append("    to the same species. For example, 'taxset Homo = 1-4;' indicates that sequences 1 to 4\n")
        bestblock.append("    belong to the species Homo.]\n\n")
        seqnames = self.getnames()
        for i in range(numpartitions):
            bestblock.append("    taxset %10s = %6d;\n" % (seqnames[i], i + 1))

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # Charset commands
        bestblock.append("    [Charset commands:]\n")
        bestblock.append("    [The location of each gene is defined by CHARSET. For example, CHARSET gene1 = 1 - 500;\n")
        bestblock.append("    says that the first 500 nucleotides belong to the gene 'gene1'. The command 'partition'\n")
        bestblock.append("    divides sequences into genes specified by CHARSET. The partition must be activated by\n")
        bestblock.append("    the command 'set partition=gene'.]\n\n")
        for i in range(numpartitions):
            name = self.partitions[i][0]
            start = self.partitions[i][1]
            stop = start + self.partitions[i][2] - 1
            bestblock.append("    charset %s = %6d - %6d;\t\t[partition no. %d]\n" % (name, start + 1, stop + 1, i + 1)) # +1: not list syntax
        bestblock.append("    partition allgenes = %d: " % numpartitions)
        for i in range(numpartitions):
            name = self.partitions[i][0]
            bestblock.append("%s" % name)
            bestblock.append(", ")
            # print 10 names per line
            if all([i%10==0, i>0, i<numpartitions - 1]):
                bestblock.append("\n        ")
        del bestblock[-1]               # Added one comma too many
        bestblock.append(";\n")
        bestblock.append("    set partition = allgenes;\n")

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # Model specification commands
        bestblock.append("    [Substitution model commands:]\n")
        bestblock.append("    [Specify substitution models for each partition. Different partitions may all use the \n")
        bestblock.append("    same type of model or different types (as in the examples below). Nucleotide models are\n")
        bestblock.append("    specified using 'lset', while amino-acid models are specified using 'prset' and 'lset'.\n")
        bestblock.append("    Modify and un-comment examples below as relevant.]\n\n")
        if numpartitions <= 2:
            set1size = 1
            set2size = 1
        else:
            set1size = 2
            set2size = numpartitions - 2
        if self.seqtype == "DNA":
            bestblock.append("    [DNA example]\n")
            bestblock.append("    [lset applyto=(")
            for i in range(set1size):
                bestblock.append(str(i + 1))
                bestblock.append(",")
            del bestblock[-1]       # Added one comma too many
            bestblock.append(") nst=6 rates=invgamma;\n")
            bestblock.append("    lset applyto=(")
            for i in range(set1size, set1size + set2size):
                bestblock.append(str(i + 1))
                bestblock.append(",")
            del bestblock[-1]       # Added one comma too many
            bestblock.append(") nst=2 rates=equal;]\n")
        elif self.seqtype == "protein":
            bestblock.append("    [Protein example]\n")
            bestblock.append("    [prset applyto=(")
            for i in range(set1size):
                bestblock.append(str(i + 1))
                bestblock.append(",")
            del bestblock[-1]       # Added one comma too many
            bestblock.append(") aamodelpr=fixed(jones);\n")
            bestblock.append("    prset applyto=(")
            for i in range(set1size, set1size + set2size):
                bestblock.append(str(i + 1))
                bestblock.append(",")
            del bestblock[-1]       # Added one comma too many
            bestblock.append(") aamodelpr=fixed(wag);\n    lset applyto=(all) rates=invgamma]\n")

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # BEST specific prset commands
        bestblock.append("    [BEST specific prset commands:]\n")
        bestblock.append("    [BEST: this parameter initiates the Bayesian analysis for estimating species trees when\n")
        bestblock.append("    setting BEST=1. If BEST=0, the regular MrBayes program is implemented. thetapr: this parameter\n")
        bestblock.append("    sets the prior distribution for population sizes. There is only one option, inverse gamma\n")
        bestblock.append("    distribution with parameters alpha and beta. The mean of the inverse gamma distribution is\n")
        bestblock.append("    beta/(alpha - 1). Users should choose reasonable values for alpha and beta such that the\n")
        bestblock.append("    prior mean of the population size theta is in a reasonable range. In the example below\n")
        bestblock.append("    thetapr=invgamma(3,0.003) implies that the prior mean of the theta is 0.0015.\n")
        bestblock.append("    genemupr: this parameter sets the prior distribution for the mutation rates across genes.\n")
        bestblock.append("    Two options: genemupr=uniform or genemupr=fixed(a).]\n\n")
        bestblock.append("    prset thetapr=invgamma(3,0.003) GeneMuPr=uniform(0.5,1.5) BEST=1;\n")

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # Unlink commands
        bestblock.append("    [Unlink commands:]\n")
        bestblock.append("    [Gene trees are assumed to be independent given the species tree. Therefore, the command\n")
        bestblock.append("    'unlink' must be used to unlink topology, branch lengths, and mutation rates across loci.]\n\n")
        bestblock.append("    unlink topology=(all) brlens=(all) genemu=(all);\n")
        bestblock.append("\n    [Users may also unlink other parameters associated with the substitution models for\n")
        bestblock.append("    each loci. Modify and un-comment example below as necessary.]\n\n")
        if self.seqtype == "DNA":
            bestblock.append("    [unlink shape=(all) statefreq=(all) tratio=(all) revmat=(all) pinvar=(all);]\n")
        elif self.seqtype == "protein":
            bestblock.append("    [unlink shape=(all) pinvar=(all);]\n")

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # MCMC commands
        bestblock.append("    [MrBayes commands:]\n")
        bestblock.append("    [Specify how to run MCMC and possibly sump and sumt as for a regular MrBayes run]\n")

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # Add final end statement and return mrbayes block as string
        bestblock.append("end;\n")
        return "".join(bestblock)


    #######################################################################################
    #######################################################################################

    def nexuspart(self):
        """Returns MrBayes block for partitioned data set as string"""

        numpartitions = len(self.partitions)

        # Header
        bestblock = ["begin mrbayes;\n    set autoclose=yes nowarn=yes;\n"]

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # Outgroup command
        bestblock.append("    [Outgroup command:]\n")
        bestblock.append("    [An outgroup has to be specified. Replace MYOUTGROUP with actual value]\n\n")
        bestblock.append("    outgroup MYOUTGROUP;\n")

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # Taxset commands
        bestblock.append("    [Taxset commands:]\n")
        bestblock.append("    [The command 'TAXSET' tells the program which sequences belong to which species.\n")
        bestblock.append("    The commands below indicate that each species has only one sequence (single allele data).\n")
        bestblock.append("    Although the species names below coincide with the sequences' names, it is totally valid\n")
        bestblock.append("    to use other names for species. For multiple allele data, multiple sequences may belong\n")
        bestblock.append("    to the same species. For example, 'taxset Homo = 1-4;' indicates that sequences 1 to 4\n")
        bestblock.append("    belong to the species Homo.]\n\n")
        seqnames = self.getnames()
        for i in range(numpartitions):
            bestblock.append("    taxset %10s = %6d;\n" % (seqnames[i], i + 1))

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # Charset commands
        bestblock.append("    [Charset commands:]\n")
        bestblock.append("    [The location of each gene is defined by CHARSET. For example, CHARSET gene1 = 1 - 500;\n")
        bestblock.append("    says that the first 500 nucleotides belong to the gene 'gene1'. The command 'partition'\n")
        bestblock.append("    divides sequences into genes specified by CHARSET. The partition must be activated by\n")
        bestblock.append("    the command 'set partition=gene'.]\n\n")
        for i in range(numpartitions):
            name = self.partitions[i][0]
            start = self.partitions[i][1]
            stop = start + self.partitions[i][2] - 1
            bestblock.append("    charset %s = %6d - %6d;\t\t[partition no. %d]\n" % (name, start + 1, stop + 1, i + 1)) # +1: not list syntax
        bestblock.append("    partition allgenes = %d: " % numpartitions)
        for i in range(numpartitions):
            name = self.partitions[i][0]
            bestblock.append("%s" % name)
            bestblock.append(", ")
            # print 10 names per line
            if all([i%10==0, i>0, i<numpartitions - 1]):
                bestblock.append("\n        ")
        del bestblock[-1]               # Added one comma too many
        bestblock.append(";\n")
        bestblock.append("    set partition = allgenes;\n")

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # Model specification commands
        bestblock.append("    [Substitution model commands:]\n")
        bestblock.append("    [Specify substitution models for each partition. Different partitions may all use the \n")
        bestblock.append("    same type of model or different types (as in the examples below). Nucleotide models are\n")
        bestblock.append("    specified using 'lset', while amino-acid models are specified using 'prset' and 'lset'.\n")
        bestblock.append("    Modify and un-comment examples below as relevant.]\n\n")

        if self.seqtype == "DNA":
            bestblock.append("    [DNA example]\n")
            bestblock.append("    [lset applyto=(")
            for i in range(numpartitions):
                bestblock.append(str(i + 1))
                bestblock.append(",")
            del bestblock[-1]       # Added one comma too many
            bestblock.append(") nst=mixed rates=invgamma;\n")

        elif self.seqtype == "protein":
            bestblock.append("    [Protein example]\n")
            bestblock.append("    [prset applyto=(")
            for i in range(numpartitions):
                bestblock.append(str(i + 1))
                bestblock.append(",")
            del bestblock[-1]       # Added one comma too many
            bestblock.append(") aamodelpr=mixed;\n    lset applyto=(all) rates=invgamma]\n")

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # BEST specific prset commands
        bestblock.append("    [BEST specific prset commands:]\n")
        bestblock.append("    [BEST: this parameter initiates the Bayesian analysis for estimating species trees when\n")
        bestblock.append("    setting BEST=1. If BEST=0, the regular MrBayes program is implemented. thetapr: this parameter\n")
        bestblock.append("    sets the prior distribution for population sizes. There is only one option, inverse gamma\n")
        bestblock.append("    distribution with parameters alpha and beta. The mean of the inverse gamma distribution is\n")
        bestblock.append("    beta/(alpha - 1). Users should choose reasonable values for alpha and beta such that the\n")
        bestblock.append("    prior mean of the population size theta is in a reasonable range. In the example below\n")
        bestblock.append("    thetapr=invgamma(3,0.003) implies that the prior mean of the theta is 0.0015.\n")
        bestblock.append("    genemupr: this parameter sets the prior distribution for the mutation rates across genes.\n")
        bestblock.append("    Two options: genemupr=uniform or genemupr=fixed(a).]\n\n")
        bestblock.append("    prset thetapr=invgamma(3,0.003) GeneMuPr=uniform(0.5,1.5) BEST=1;\n")

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # Unlink commands
        bestblock.append("    [Unlink commands:]\n")
        bestblock.append("    [Gene trees are assumed to be independent given the species tree. Therefore, the command\n")
        bestblock.append("    'unlink' must be used to unlink topology, branch lengths, and mutation rates across loci.]\n\n")
        bestblock.append("    unlink topology=(all) brlens=(all) genemu=(all);\n")
        bestblock.append("\n    [Users may also unlink other parameters associated with the substitution models for\n")
        bestblock.append("    each loci. Modify and un-comment example below as necessary.]\n\n")
        if self.seqtype == "DNA":
            bestblock.append("    [unlink shape=(all) statefreq=(all) tratio=(all) revmat=(all) pinvar=(all);]\n")
        elif self.seqtype == "protein":
            bestblock.append("    [unlink shape=(all) pinvar=(all);]\n")

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # MCMC commands
        bestblock.append("    [MrBayes commands:]\n")
        bestblock.append("    [Specify how to run MCMC and possibly sump and sumt as for a regular MrBayes run]\n")

        bestblock.append("\n\n    [######################################################################################]\n\n")

        # Add final end statement and return mrbayes block as string
        bestblock.append("end;\n")
        return "".join(bestblock)


    #######################################################################################

##    def alignformat(self, width=60):
##        """Only valid for two aligned sequences: Returns alignment with formatting a la fasta align programs"""
##
##        numseq = len(self)
##        if numseq != 2:
##            raise SeqError("Alignformat can only be used with exactly two sequences")
##        alignlen = self.alignlen()
##
##        # If format is not protein (presumably "DNA") then conservation line is limited to scoring conserved
##        # versus non-conserved. If format is "protein" then conservation line should score
##        # strong and weak conservation also.
##
##        # Format is "protein"
##        if self.seqtype == "protein":
##
##            similar_pair_strings = ["CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK", "NEQHRK", "FVLIM", "FVLIM"]
##            semicons_sub = []
##            for string in semicons_sub_strings:
##                semicons_sub.append(set(string))
##
##            # Function for checking whether "query_set" is a subset of one of the symbol sets in "set_list"
##            def containsset(query_set, set_list):
##                for symbol_set in set_list:
##                    if query_set <= symbol_set:
##                        return True
##
##                # If we fell off for loop then set is not contained
##                return False
##
##
##    ('W', 'F') : 1
##    ('V', 'M') : 1
##    ('Z', 'K') : 1
##    ('D', 'N') : 1
##    ('Q', 'R') : 1
##    ('Y', 'F') : 3
##    ('V', 'L') : 1
##    ('K', 'R') : 2
##    ('E', 'D') : 2
##    ('M', 'I') : 1
##    ('Z', 'D') : 1
##    ('B', 'D') : 4
##    ('H', 'N') : 1
##    ('S', 'N') : 1
##    ('E', 'Q') : 2
##    ('Y', 'W') : 2
##    ('K', 'Q') : 1
##    ('M', 'L') : 2
##    ('K', 'E') : 1
##    ('Z', 'E') : 4
##    ('Z', 'Q') : 3
##    ('B', 'E') : 1
##    ('S', 'A') : 1
##    ('Y', 'H') : 2
##    ('T', 'S') : 1
##    ('L', 'I') : 2
##    ('Z', 'B') : 1
##    ('B', 'N') : 3
##
##

#############################################################################################
#############################################################################################

class Seqfile_reader(object):
    """Baseclass for sequence file readers. Don't instantiate"""

    def __init__(self, filename, seqtype, check_alphabet, degap, nameishandle):

        # If nameishandle==True, directly assign "filename" as filehandle (should I do error checking?)
        if nameishandle:
            self.filename = "handle"
            self.seqfile = filename

        # Special filename "-" indicates stdin.
        # Read all data from stdin, place it in virtual (RAM) file, and use that instead of real file
        elif filename == "-":
            self.filename = "stdin"
            data = sys.stdin.read()
            self.seqfile = StringIO(data)

        # Regular file
        else:
            self.filename = filename
            self.seqfile = open(self.filename, "r")

        self.seqtype = seqtype
        self.check_alphabet = check_alphabet
        self.degap = degap

        # Translation tables for removing whitespace and numbering
        # Note: this is an empty table, since I will only be using the "deletechars" functionality in .translate
        self.spacetrans = str.maketrans("", "", string.whitespace)
        self.alltrans = str.maketrans("", "", string.whitespace + "0123456789")

    #######################################################################################

    def makeseq(self, name, seq, annotation="", comments=""):
        """Takes name, sequence, annotation, and comments, returns sequence object of correct type.
        Called by subclass"""

        seq = seq.upper()
        if self.seqtype == "autodetect":
            self.seqtype = find_seqtype(seq)
        if self.seqtype == "standard":
            return Standard_sequence(name, seq, annotation, comments, self.check_alphabet, self.degap)
        elif self.seqtype == "DNA":
            return DNA_sequence(name, seq, annotation, comments, self.check_alphabet, self.degap)
        elif self.seqtype == "protein":
            return Protein_sequence(name, seq, annotation, comments, self.check_alphabet, self.degap)
        elif self.seqtype == "ASCII":
            return ASCII_sequence(name, seq, annotation, comments, self.check_alphabet, self.degap)
        else:
            raise SeqError("Unknown sequence type: %s" % self.seqtype)

    #######################################################################################

    def readseq(self):
        """Reads single sequence from file, returns as sequence object"""
        return next(self)

    #######################################################################################

    def read_seqs(self, silently_discard_dup_name=False):
        """Reads all sequences, returns them as Seq_set object"""
        if self.filename == "stdin" or self.filename == "handle":
            name = "Partition"
        else:
            name = os.path.split(self.filename)[1].split(".")[0]    # Discard path and extension from file name
        seqset = Seq_set(name)
        for seq in self:
            seqset.addseq(seq, silently_discard_dup_name)
        return seqset

    #######################################################################################

    def read_alignment(self, silently_discard_dup_name=False):
        """Reads all sequences, returns them as Seq_alignment object"""

        if self.filename == "stdin" or self.filename == "handle":
            name = "Partition"
        else:
            name = os.path.split(self.filename)[1].split(".")[0]    # Discard path and extension from file name
        alignment = Seq_alignment(name)
        for seq in self:
            alignment.addseq(seq, silently_discard_dup_name)
        return alignment

#############################################################################################
#############################################################################################

class Fastafilehandle(Seqfile_reader):
    """Reader class for fasta files"""

    def __init__(self, filename, seqtype="autodetect", check_alphabet=False, degap=False, nameishandle=False):
        Seqfile_reader.__init__(self, filename, seqtype, check_alphabet, degap, nameishandle)

        # Perform "magic number" check of whether file appears to be in fasta format
        line = self.seqfile.readline()
        if not line.startswith(">"):
            raise SeqError("File '%s' does not appear to be in FASTA format" % self.filename)

        # If everything was kosher: move filepointer back to beginning of file
        else:
            self.seqfile.seek(0)

    #######################################################################################

    def __iter__(self):
        return self

    #######################################################################################

    def __next__(self):
        """Parser function, returns next seq as Sequence object"""

        # First time this function is entered, we have to read next line from file
        # All other times, we have next line in buffer (because function reads one line too far)
        try:
            line = self.buffer
        except AttributeError:
            line = self.seqfile.readline()

        # If next line is empty, then EOF has been reached, and iteration should stop
        if len(line)==0:
            self.seqfile.close()
            raise StopIteration()
        else:
            nameline = line

        # Get seq.: read until next nameline or EOF.
        seqlist = []
        line = self.seqfile.readline()
        while (">" not in line) and line:
            seqlist.append(line)
            line = self.seqfile.readline()

        # We have read one line too far, save in buffer for next time
        self.buffer = line

        # Sequence has been read, now parse text and return Sequence object,
        words = nameline.split()
        name = words[0].replace(">","")
        if len(words) > 1:
            comments = " ".join(words[1:])
        else:
            comments = ""
        annotation = ""         # Note: annotation not possible in fasta format, set to empty string
        seq = "".join(seqlist)
        seqletters = set(seq)
        if seqletters <= set(["0","1","\n"]):   # Only 0 and 1: this is gapencoding => do NOT remove numbers!!!
            seq = seq.translate(self.spacetrans)            # Remove whitespace
        else:
            seq = seq.translate(self.alltrans)              # Remove whitespace and numbering

        # Use method in baseclass to determine seqtype and construct proper sequence object
        seqobject = Seqfile_reader.makeseq(self, name, seq, annotation, comments)
        return seqobject

#############################################################################################
#############################################################################################

class Howfilehandle(Seqfile_reader):
    """Reader class for HOW files"""

    #######################################################################################

    def __init__(self, filename, seqtype="autodetect", check_alphabet=False, degap=False, nameishandle=False):
        Seqfile_reader.__init__(self, filename, seqtype, check_alphabet, degap, nameishandle)

        # Perform minimal "magic number" check of whether file appears to be in HOW format
        line = self.seqfile.readline()
        words = line.split()
        # Note: testing for line[6]!=" " could fail for seqlen>=1E6 (but perhaps not problem)
        if any([len(words) != 2, not words[0].isdigit(), line[6] != " "]):
            raise SeqError("File '{}' does not appear to be in HOW format".format(self.filename))

        # If everything was kosher: move filepointer back to beginning of file
        else:
            self.seqfile.seek(0)

    #######################################################################################

    def __iter__(self):
        return self

    #######################################################################################

    def __next__(self):
        """Parser function, returns next seq as Sequence object"""

        # First time this function is entered, we have to read next line from file
        # All other times, we have next line in buffer (because function reads one line too far)
        try:
            line = self.buffer
        except AttributeError:
            line = self.seqfile.readline()

        # If next line is empty, then EOF has been reached, and iteration should stop
        if len(line)==0:
            self.seqfile.close()
            raise StopIteration()
        else:
            nameline = line

        # Get seq + annotation: read until next nameline or EOF.
        seqlist = []
        line = self.seqfile.readline()
        while not line.startswith(" ") and line:
            seqlist.append(line)
            line = self.seqfile.readline()

        # We have read one line too far, save in buffer for next time
        self.buffer = line

        # Sequence and annotation has been read, now parse text and return Sequence object,
        words = nameline.split()
        seqlen = int(words[0])
        name = words[1]
        comments = " ".join(words[2:])
        seqann = "".join(seqlist)
        seqann = seqann.translate(self.alltrans)        # Remove whitespace and numbering
        seq = seqann[:seqlen]
        annotation = seqann[seqlen:]

        # Use method in baseclass to determine seqtype and construct proper sequence object
        seqobject = Seqfile_reader.makeseq(self, name, seq, annotation, comments)
        return seqobject

#############################################################################################
#############################################################################################

class Genbankfilehandle(Seqfile_reader):
    """Reader class for GenBank files"""

    def __init__(self, filename, seqtype="autodetect", check_alphabet=False, degap=False, nameishandle=False,
                    namefromfields=None):

        self.namefromfields = namefromfields
        Seqfile_reader.__init__(self, filename, seqtype, check_alphabet, degap, nameishandle)

        # Perform "magic number" check of whether file appears to be in GenBank format
        # Note: In genbank databasefiles LOCUS is not on first line. Alter?
        line = self.seqfile.readline()
        if not line.startswith("LOCUS"):
            raise SeqError("File '%s' does not appear to be in GenBank format" % self.filename)

        # If everything was kosher: move filepointer back to beginning of file
        # initialise line to placeholder
        else:
            self.seqfile.seek(0)
            self.line = "placeholder"

    #######################################################################################

    def __iter__(self):
        return self

    #######################################################################################

    def __next__(self):
        """Parser function, returns next seq as Sequence object"""

        # Locate next LOCUS line (= start of sequence).
        # If EOF reached before finding LOCUS: raise StopIteration and close seqfile
        # After this function filepointer is on line starting with LOCUS
        self.locusname,self.seqlen = self.find_LOCUS()

        # Read all non-sequence info into list of lines, for subsequent parsing
        # After this function filepointer is on line starting with ORIGIN
        metadata = self.read_metadata()

        # Parse metadata to extract annotation info
        annotation,comments = self.extract_annotation(metadata)

        # Parse metadata to extract sequence name
        seqname = self.extract_name(metadata)
        del metadata #probably not worth it to save memory, even if full genome...

        # Read sequence
        seq = self.read_seq()

        # Use method in baseclass to determine seqtype and construct proper sequence object
        seqobject = Seqfile_reader.makeseq(self, seqname, seq, annotation, comments)
        return seqobject

    #######################################################################################

    def find_LOCUS(self):
        # Find start of next sequence record, halt iteration i EOF encountered first
        while not self.line.startswith("LOCUS"):
            self.line = self.seqfile.readline()
            if len(self.line)==0:
                self.seqfile.close()
                raise StopIteration()

        # We are now at LOCUS line at top of sequence record
        # derive sequence type if possible, store ID for possible use as name
        words = self.line.split()
        locusname = words[1]
        seqlen = int(words[2])
        molecule = words[4]
        if molecule == "DNA":
            self.seqtype = "DNA"
        elif molecule == "PRT":
            self.seqtype = "protein"
        else:
            self.seqtype = "autodetect"

        return locusname,seqlen

    #######################################################################################

    def read_metadata(self):
        metadatalist = []
        while not self.line.startswith("ORIGIN"):
            self.line = self.seqfile.readline()
            metadatalist.append(self.line)
        return metadatalist

    #######################################################################################

    def extract_annotation(self, metadata):
        # TO BE IMPLEMENTED
        annotation = ""
        commments= ""
        return annotation,commments

    #######################################################################################

    def extract_name(self, metadata):
        if self.namefromfields:
            seqname = ""
            fieldlist = self.namefromfields.split(",")
            for field in fieldlist:
                seen = False
                for line in metadata:
                    words = line.split()
                    if words[0] == field:
                        if seqname:
                            seqname += "_"
                        seqname += "_".join(words[1:])
                        seen = True
                        break
                if not seen:
                    raise SeqError("This field was not seen in GenBank file: '{}'".format(field))
        else:
            seqname = self.locusname

        return seqname

    #######################################################################################

    def read_seq(self):
        # File pointer is currently at ORIGIN: sequence starts on next line
        seqlist = []
        line = self.seqfile.readline()
        while not line.startswith("//"):
            words = line.split()
            seqlist.append("".join(words[1:]))
            line = self.seqfile.readline()
        seq = "".join(seqlist)
        return seq

#############################################################################################
#############################################################################################

class Tabfilehandle(Seqfile_reader):
    """Reader class for TAB-files"""

    def __init__(self, filename, seqtype="autodetect", check_alphabet=False, degap=False, nameishandle=False):
        Seqfile_reader.__init__(self, filename, seqtype, check_alphabet, degap, nameishandle)

        # Perform "magic number" check of whether file appears to be in TAB format
        line = self.seqfile.readline()
        words = line.split("\t")        # Split on TAB - also handles multiple, consecutive TABS
        if not len(words)>=2:
            raise SeqError("File '%s' does not appear to be in TAB format" % self.filename)

        # If everything was kosher: move filepointer back to beginning of file
        else:
            self.seqfile.seek(0)

    #######################################################################################

    def __iter__(self):
        return self

    #######################################################################################

    def __next__(self):
        """Parser function, returns next seq as Sequence object"""

        # Note: This function is required in combination with baseclass __iter__() for iteration

        # If next line is empty, then EOF has been reached, and iteration should stop
        line = self.seqfile.readline()
        if len(line)==0:
            self.seqfile.close()
            raise StopIteration()

        line = line.rstrip()
        words = line.split("\t")        # Splitting on tab: strict parsing
        numwords = len(words)

        if numwords == 1:
            raise SeqError("This line does not appear to be in TAB format: %s" % line)
        else:
            name = words[0]
            seq = words[1]
            annotation = comments = ""
            if numwords >= 3:
                annotation = words[2]
            if numwords == 4:
                comments = words[3]

        # Use method in baseclass to determine seqtype and construct proper sequence object
        seqobject = Seqfile_reader.makeseq(self, name, seq, annotation, comments)
        return seqobject

#############################################################################################
#############################################################################################

class Rawfilehandle(Seqfile_reader):
    """Reader class for RAW-files"""

    def __init__(self, filename, seqtype="autodetect", check_alphabet=False, degap=False, nameishandle=False):

        # Counter used to assign continuously numbered names to RAW sequences
        self.cur_seq_no = 1

        Seqfile_reader.__init__(self, filename, seqtype, check_alphabet, degap, nameishandle)

    #######################################################################################

    def __iter__(self):
        return self

    #######################################################################################

    def __next__(self):
        """Parser function, returns next seq as Sequence object"""

        # Note: This function is required in combination with baseclass __iter__() for iteration

        # If next line is empty, then EOF has been reached, and iteration should stop
        line = self.seqfile.readline()
        if len(line)==0:
            self.seqfile.close()
            raise StopIteration()

        # Parser constructed to accept practcally anything with no checking:
        # All non-whitespace on a line is concatenated and interpreted as a sequence
        # Names are: seq1, seq2, seq3, ...
        words = line.split()
        seq = "".join(words)
        name = "seq%d" % (self.cur_seq_no)
        self.cur_seq_no += 1

        # Use method in baseclass to determine seqtype and construct proper sequence object
        seqobject = Seqfile_reader.makeseq(self, name, seq)
        return seqobject

#############################################################################################
#############################################################################################

class Alignfile_reader(object):
    """Baseclass for alignment-type sequence file readers"""

    # Note: Unlike non-alignment sequence files, alignment files do not have iteration defined:
    # You typically will have to read the entire alignment before processing

    def __init__(self, filename, seqtype, check_alphabet, degap, nameishandle):

        # If nameishandle==True, directly assign this as filehandle (should I do error checking?)
        if nameishandle:
            self.seqfile = filename

        # Special filename "-" indicates stdin.
        # Read all data from stdin, place it in virtual (RAM) file, and use that instead of real file
        elif filename == "-":
            self.filename = "stdin"
            data = sys.stdin.read()
            self.seqfile = StringIO(data)

        # Regular file
        else:
            self.filename = filename
            self.seqfile = open(self.filename, "r")

        self.seqtype = seqtype
        self.check_alphabet = check_alphabet
        self.degap = degap

        # Read all data into object, so it can be parsed for single sequences, close file
        self.seqdata = self.seqfile.read()
        self.seqfile.close()

        # Translation tables for removing whitespace and numbering
        # Note: this is an empty table, since I will only be using the "deletechars" functionality in .translate
        self.spacetrans = str.maketrans("", "", string.whitespace)
        self.alltrans = str.maketrans("", "", string.whitespace + "0123456789")

    #######################################################################################

    def makeseq(self, name, seq, annotation="", comments=""):
        """Takes name, sequence, annotation, and comments, returns sequence object of correct type"""

        seq = seq.upper()
        if self.seqtype == "autodetect":
            self.seqtype = find_seqtype(seq)
            if self.seqtype == "standard":
                return Standard_sequence(name, seq, annotation, comments, self.check_alphabet, self.degap)
            elif self.seqtype == "DNA":
                return DNA_sequence(name, seq, annotation, comments, self.check_alphabet, self.degap)
            elif self.seqtype ==  "protein":
                return Protein_sequence(name, seq, annotation, comments, self.check_alphabet, self.degap)
            elif self.seqtype == "ASCII":
                return ASCII_sequence(name, seq, annotation, comments, self.check_alphabet, self.degap)
            else:
                unknown_symbols = list(seqletters - Const.ASCII)
                msg = "Unknown symbols encountered during seqtype autodetection: {}".format(unknown_symbols)
                raise SeqError(msg)
        elif self.seqtype == "standard":
            return Standard_sequence(name, seq, annotation, comments, self.check_alphabet, self.degap)
        elif self.seqtype == "DNA":
            return DNA_sequence(name, seq, annotation, comments, self.check_alphabet, self.degap)
        elif self.seqtype == "protein":
            return Protein_sequence(name, seq, annotation, comments, self.check_alphabet, self.degap)
        elif self.seqtype == "ASCII":
            return ASCII_sequence(name, seq, annotation, comments, self.check_alphabet, self.degap)
        else:
            raise SeqError("Unknown sequence type: %s" % self.seqtype)

    #######################################################################################

    def read_seqs(self, silently_discard_dup_name=False):

        aligned_seqs = self.read_alignment(silently_discard_dup_name)
        return aligned_seqs.seqset()

#############################################################################################
#############################################################################################

class Clustalfilehandle(Alignfile_reader):
    """Reader class for Clustal alignment files"""

    def __init__(self, filename, seqtype="autodetect", check_alphabet=False, degap=False, nameishandle=False):
        Alignfile_reader.__init__(self, filename, seqtype, check_alphabet, degap, nameishandle)

        # Split input data into list of lines
        self.seqdata = self.seqdata.split("\n")

        # Perform "magic number" check of whether file appears to be in clustal format
        line = self.seqdata[0]
        if not line.startswith("CLUSTAL"):
            raise SeqError("File '%s' does not appear to be in Clustal format" % self.filename)

    #######################################################################################

    def read_alignment(self, silently_discard_dup_name=False):
        """Reads all sequences, returns them as Seq_alignment object"""

        if self.filename == "stdin" or self.filename == "handle":
            name = "Partition"
        else:
            name = os.path.split(self.filename)[1].split(".")[0]    # Discard path and extension from file name
        alignment = Seq_alignment(name)
        # line no. 0 contains CLUSTAL header. Find first non-blank line after this
        lineno = 1
        while  len(self.seqdata[lineno].strip()) == 0:
            lineno += 1

        # Iterate over remaining lines, adding sequences to dictionary as we go
        seqdict = {}
        for i in range(lineno,len(self.seqdata)):
            line = self.seqdata[i]

            # If line is not empty and does not start with whitespace, then it must contain sequence info
            if len(line) and not line[0].isspace():
                words = line.split()
                name = words[0]
                seq = "".join(words[1:])        # join: in case sequence is divided into blocks
                seq = re.sub("[0-9]*", "", seq) # remove numbering if present

                # If name has been seen before: add next seq fragment to existing list of strings
                if name in seqdict:
                    seqdict[name].append(seq)

                # If this is first sighting: add entry to dictionary
                else:
                    seqdict[name] = [seq]

        # For each entry in sequence dictionary:  Join list of strings to single string,
        # convert to Sequence object of proper type, add Sequence object to alignment object
        for name in seqdict:
            seq = "".join(seqdict[name])
            seqobject = Alignfile_reader.makeseq(self, name, seq)
            alignment.addseq(seqobject, silently_discard_dup_name)

        return alignment

#############################################################################################
#############################################################################################

class Phylipfilehandle(Alignfile_reader):
    """Reader class for Phylip alignment files"""

    def __init__(self, filename, seqtype="autodetect", check_alphabet=False, degap=False, nameishandle=False):
        Alignfile_reader.__init__(self, filename, seqtype, check_alphabet, degap, nameishandle)

        # Split input data into list of lines
        self.seqdata = self.seqdata.split("\n")

        # Perform "magic number" check of whether file appears to be in phylip format
        line = self.seqdata[0]
        words = line.split()
        if not (line[0].isspace() and len(words) >= 2):
            raise SeqError("File '%s' does not appear to be in Phylip format" % self.filename)

    #######################################################################################

    def read_alignment(self, silently_discard_dup_name=False):
        """Reads all sequences, returns them as Seq_alignment object"""

        if self.filename == "stdin" or self.filename == "handle":
            name = "Partition"
        else:
            name = os.path.split(self.filename)[1].split(".")[0]    # Discard path and extension from file name
        alignment = Seq_alignment(name)

        # Line no. 0 contains Phylip header. Extract number of sequences and alignment length for later checks
        words = self.seqdata[0].split()
        nseqs = int(words[0])
        alignlen = int(words[1])

        # Read sequences one line at a time, discarding blank lines, constantly keeping track of names
        # NOTE: could make this more efficient if I made more assumptions about file structure
        #   (e.g., same order of names in consecutive blocks)
        i = 0
        seqdict = {}
        for line in self.seqdata[1:]:
            if line.strip():   # If line is not empty
                words = line.split()
                name = words[0]
                seq = "".join(words[1:])
                seqletters = set(seq)
                if seqletters <= set(["0","1","\n"]):   # Only 0 and 1: this is gapencoding => do NOT remove numbers!!!
                    seq = seq.translate(self.spacetrans)            # Remove whitespace
                else:
                    seq = seq.translate(self.alltrans)              # Remove whitespace and numbering
                if name in seqdict:
                    seqdict[name].append(seq)
                else:
                    seqdict[name] = [seq]

        # For each entry in sequence dictionary:  Join list of strings to single string,
        # convert to Sequence object of proper type, add Sequence object to alignment object
        for name,seqlist in seqdict.items():
            seq = "".join(seqlist)
            seqobject = Alignfile_reader.makeseq(self, name, seq)
            alignment.addseq(seqobject, silently_discard_dup_name)

        if len(alignment) != nseqs:
            raise SeqError("Number of sequences in phylip file ({}) does not match header info ({})".format(len(alignment) , nseqs))
        if len(alignment[0]) != alignlen:
            raise SeqError("Number of sequences in phylip file ({}) does not match header info ({})".format(len(alignment[0]) , alignlen))

        return alignment

#############################################################################################
#############################################################################################

class Nexusfilehandle(Alignfile_reader):
    """Reader class for NEXUS alignment files"""

    def __init__(self, filename, seqtype="autodetect", check_alphabet=False, degap=False, nameishandle=False):
        Alignfile_reader.__init__(self, filename, seqtype, check_alphabet, degap, nameishandle)

        # Perform "magic number" check of whether file appears to be in NEXUS format
        if not self.seqdata[:6].lower() == "#nexus":
            raise SeqError("File '%s' does not appear to be in NEXUS format" % self.filename)

        # Remove comments
        self.seqdata =  remove_comments(self.seqdata, leftdelim = "[", rightdelim="]")

        # Split input data into list of lines
        self.seqdata = self.seqdata.split("\n")

    #######################################################################################

    def read_alignment(self, silently_discard_dup_name=False):
        """Reads all sequences, returns them as Seq_alignment object"""

        if self.filename == "stdin" or self.filename == "handle":
            name = "Partition"
        else:
            name = os.path.split(self.filename)[1].split(".")[0]    # Discard path and extension from file name
        alignment = Seq_alignment(name)

        # Ignore all header information: read until after "matrix" statement, then parse sequence blocks
        # NOTE: it is possible this will crash if NEXUS file also contains a matrix of non-sequence characters
        lineno = 1
        while self.seqdata[lineno].lower().strip() != "matrix":
            lineno += 1
        lineno += 1

        # Iterate over remaining lines, adding sequences to dictionary as we go
        seqdict = {}
        for i in range(lineno,len(self.seqdata)):
            line = self.seqdata[i]

            # If semicolon reached, then entire matrix block has been read: stop iteration
            if ";" in line:
                break

            # Skip blank lines, parse lines with text
            words = line.split()
            if len(words)>0:
                name = words[0]
                seq = "".join(words[1:])        # join: in case sequence is divided into blocks
                seq = seq.replace("?", "-")     # "?" sometimes used for gap or unknown in NEXUS
                seqletters = set(seq)
                seq = seq.translate(self.spacetrans)            # Remove whitespace

                # If name has been seen before: add next seq fragment to existing list of strings
                if name in seqdict:
                    seqdict[name].append(seq)

                # If this is first sighting: add entry to dictionary
                else:
                    seqdict[name] = [seq]

        # For each entry in sequence dictionary:  Join list of strings to single string,
        # convert to Sequence object of proper type, add Sequence object to alignment object
        for name in seqdict:
            seq = "".join(seqdict[name])
            seqobject = Alignfile_reader.makeseq(self, name, seq)
            alignment.addseq(seqobject, silently_discard_dup_name)

        return alignment

#############################################################################################
#############################################################################################

class Stockholmfilehandle(Alignfile_reader):
    """Reader class for STOCKHOLM format alignment files"""

    # NOTE: currently ignores metadata

    def __init__(self, filename, seqtype="autodetect", check_alphabet=False, degap=False, nameishandle=False):
        Alignfile_reader.__init__(self, filename, seqtype, check_alphabet, degap, nameishandle)

        # Perform "magic number" check of whether file appears to be in STOCKHOLM format
        if not self.seqdata[:11].lower() == "# stockholm":
            raise SeqError("File '{}' does not appear to be in STOCKHOLM format".format(self.filename))

        # Split input data into list of lines
        self.seqdata = self.seqdata.split("\n")

    #######################################################################################

    def read_alignment(self, silently_discard_dup_name=False):
        """Reads all sequences, returns them as Seq_alignment object"""

        # Meta data is present on lines starting with #. Currently these lines are ignored
        # Sequence data can be written over several lines (even though original standard
        # specified one sequence per line I think?), with each line being: <name> <seq>
        # Sequence block ends with "//"

        # In alignments from HMMer's hmmalign insert states will have lower case residues
        # and/or "." for gaps (not really gaps - just lack of insert in other seq).
        # I keep track of these in alignment's annotation field (not seq annotation)

        if self.filename == "stdin" or self.filename == "handle":
            name = "Partition"
        else:
            name = os.path.split(self.filename)[1].split(".")[0]    # Discard path and extension from file name
        alignment = Seq_alignment(name)
        seqdict = {}
        for line in self.seqdata:
            if line.strip():                            # If line is not empty:
                words = line.rstrip().split()
                if words[0].startswith("#"):
                    pass                               # Skip metadata
                elif words[0].startswith("//"):
                    break                              # End loop on //
                else:
                    name = words[0]
                    seq = "".join(words[1:])           # join: in case sequence is divided into blocks
                    if name in seqdict:
                        seqdict[name].append(seq)
                    else:
                        seqdict[name] = [seq]

        # For each entry in sequence dictionary: Join list of strings to single string,
        for name,seqlist in seqdict.items():
            seq = "".join(seqlist)
            seqdict[name] = seq

        # Extract information about sites corresponding to insert states
        # Note: I only have to look at one sequence:
        # insert states will be either lowercase or "." in each individual sequence
        firstseq = next(iter(seqdict.values()))
        annotation = ["m"] * len(firstseq)
        for i,char in enumerate(firstseq):
            if char == "." or char.islower():
                annotation[i] = "i"
        annotation = "".join(annotation)

        # convert to Sequence object of proper type, add Sequence object to alignment object
        # Note: alignment annotation is added as annotation field to eacg seq. Consider alternatives...
        for name,seq in seqdict.items():
            seq = seq.replace(".", "-")         # Replace "." gaps with regular gapsymbols
            seqobject = Alignfile_reader.makeseq(self, name, seq, annotation=annotation)
            alignment.addseq(seqobject, silently_discard_dup_name)

        alignment.annotation = annotation
        return alignment

#############################################################################################
#############################################################################################

class Seqfile(object):
    """Factory for constructing seqfilehandle objects. Type selected based on fileformat string"""

    # This is handy for avoiding long "elif" statements when selecting among several filehandles
    # Also implements filetype autodetection (should that go somewhere else?)
    def __new__(klass, filename, filetype="autodetect", seqtype="autodetect", check_alphabet=False,
                degap=False, nameishandle=False):

        known_handles = {"raw": Rawfilehandle, "tab": Tabfilehandle, "fasta": Fastafilehandle,
                   "nexus": Nexusfilehandle, "phylip": Phylipfilehandle, "clustal": Clustalfilehandle,
                    "genbank": Genbankfilehandle, "how": Howfilehandle, "stockholm": Stockholmfilehandle}

        # Return requested handle if among known formats
        if filetype in known_handles:
            return known_handles[filetype](filename, seqtype, check_alphabet, degap, nameishandle)

        # If autodetection of "filetype" is requested:
        # Perform "magic number" test based on first line in file. Likely to not always work...
        elif filetype == "autodetect":

            # Special filename "-" indicates stdin. Autodetection on stdin is complicated, since stdin can only be read once
            # Solution: read all data into memory, place it in virtual (RAM) file and directly pass handle to that.
            if filename == "-":
                data = sys.stdin.read()
                stdinfile = StringIO(data)
                firstline = stdinfile.readline()
                stdinfile.seek(0)
                filename = stdinfile
                nameishandle = True
            else:
                tempfile = open(filename)
                firstline = tempfile.readline()
                tempfile.close()

            if firstline.startswith(">"):
                filetype = "fasta"
            elif firstline.lower().startswith("#nexus"):
                filetype = "nexus"
            elif firstline.lower().startswith("# stockholm"):
                filetype = "stockholm"
            elif (firstline[0].isspace() and len(firstline.split()) >= 2):
                filetype = "phylip"
            elif firstline.startswith("CLUSTAL"):
                filetype = "clustal"
            elif len(firstline.split("\t")) >= 2:
                filetype = "tab"
            elif firstline.startswith("LOCUS"):
                filetype = "genbank"
            # NOTE: something wrong with how vs phylip autodetection - repair!
            elif (len(firstline.split()) >= 2) and (firstline.split()[0].isdigit) and firstline[6] == " ":
                filetype = "how"
            else:
                # If all other tests did not work: raise SeqError
                raise SeqError("Cannot autodetect type of sequence file: %s" % filename)
            return known_handles[filetype](filename, seqtype, check_alphabet, degap, nameishandle)

        else:
            raise SeqError("Uknown sequence file type: %s" % filetype)

#############################################################################################
#############################################################################################

class vcf(object):

    # NOTE: rewrite to base parsing on actual info in INFO header.
    # Or perhaps simply to take into account most common fields in INFO:
    #   https://en.wikipedia.org/wiki/Variant_Call_Format

    def __init__(self):
        self.snpdict = {}
        self.refdict = {}
        self.snpposcache = None
        self.refseq = None

    #######################################################################################

    # VCF file format (header lines + example of SNP)
    # ##INFO=<ID=AF,Number=1,Type=Float,Description="allele frequency for ALT allele">
    # ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
    # ##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality (i.e. after filtering) ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
    # ##INFO=<ID=SB,Number=1,Type=Integer,Description="Phred-scaled strand bias at this position">
    # ##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description="Indicates that the variant is a consensus variant (as opposed to a low frequency variant).">
    # #CHROM                   POS    ID  REF  ALT  QUAL  FILTER     INFO
    # SA13-JFH1-recombinant    327    .    G    A    136    .        AF=0.001661;DP=109594;DP4=104565,4793,155,6;SB=1

    @classmethod
    def from_vcffile(cls, filename, pos1=1, include_min=None, include_max=None):
        """Constructor 1: constructs vcf object from standard VCF file"""
        # pos1: for mapping between VCF file POS and some external numbering system.
        #       pos1 = VCF POS that is position 1 in outside system
        # include_min, include_max: optional cutoffs in sequence below and above which to discard SNP info

        self = cls()
        with open(filename) as myf:
            for line in myf:
                if not line.startswith("#"):
                    words = line.split()

                    # Compute position according to outside numbering scheme
                    # Discard positions if requested
                    pos = int(words[1]) - pos1 + 1          # NOTE: should I use python numbering here? need to think..
                    if (include_min and pos < include_min) or (include_max and pos > include_max):
                        pass
                    else:
                        refnuc = words[3]
                        altnuc = words[4]
                        if pos not in self.snpdict:
                            self.snpdict[pos] = {   "A":{"count":0, "freq":0.0},
                                                    "C":{"count":0, "freq":0.0},
                                                    "G":{"count":0, "freq":0.0},
                                                    "T":{"count":0, "freq":0.0}}
                            self.snpdict[pos][refnuc]["freq"] = 1.0     # Reference frequency has to be computed by subtracting ALT frequency/ies from 1.0
                            self.refdict[pos] = refnuc

                        infowords = words[7].split(";")
                        countstr = infowords[2].replace("DP4=", "")
                        countwords = countstr.split(",")
                        nref = int(countwords[0]) + int(countwords[1])
                        nalt = int(countwords[2]) + int(countwords[3])
                        self.snpdict[pos][refnuc]["count"] = nref       # Note: if pos is repeated, then ref count is given again, so no harm done here
                        self.snpdict[pos][altnuc]["count"] = nalt       # Each alternative nuc is only mentioned on one line, so again no need to check previous value

                        falt = float(infowords[0].replace("AF=", ""))
                        self.snpdict[pos][altnuc]["freq"] = falt
                        self.snpdict[pos][refnuc]["freq"] -= falt
        return self

    #######################################################################################

    # Homegrown, pre-parsed format looking as follows:
    # #POS  A   C   G   T
    # 6 0   48517   0   119
    # 9 541 0   93756   1275
    # 27    49090   0   163 0
    # 46    49705   0   323 0
    # 49    156 50469   0   0

    @classmethod
    def from_vcfcountfile(cls, filename):
        """Constructor 2: constructs vcf object from VCF counts file"""
        self = cls()
        with open(filename) as myf:
            for line in myf:
                if not line.startswith("#"):
                    words = line.split()
                    pos = int(words[0])
                    self.snpdict[pos] = {"A":int(words[1]), "C":int(words[2]), "G":int(words[3]), "T":int(words[4])}

        return self

    #######################################################################################

    def addrefseq(self, seq):
        # Adding reference sequence allows object to return snpvectors also for constant sites
        # All counts are then set to REF nuc, and count is set to minimum value observed
        countlist = []
        for pos in self.snpposlist():
            countlist.append(sum(self.snpvector(pos, ascounts=True)))
        self.mincount = min(countlist)
        self.refseq = seq
        self.nuc2int = {"A": 0, "C": 1, "G": 2, "T": 3}
        # Should have error check for whether sequence is shorter than max pos

    #######################################################################################

    def snpvector(self, pos, ascounts=False):
        # Return vector of nucleotide frequencies (or counts) + indication of REF nucleotide
        if pos in self.snpposlist():
            s = self.snpdict[pos]
            if ascounts:
                t = "count"
            else:
                t = "freq"
            return [s["A"][t], s["C"][t], s["G"][t], s["T"][t]]
        elif self.refseq:
            refnuc = self.refseq[pos-1]
            refindex = self.nuc2int[refnuc]
            fakevec = [0, 0, 0, 0]
            if ascounts:
                fakevec[refindex] = self.mincount
            else:
                fakevec[refindex] = 1.0
            return fakevec
        else:
            raise SeqError("The VCF file contains no information about position {}".format(pos))

    #######################################################################################

    def snpvectorlist(self, poslist, ascounts=False):
        # Return list of SNP vectors for requested positions
        poslist.sort()
        countarray = []
        for pos in poslist:
            countarray.append(self.snpvector(pos, ascounts))
        return(countarray)

    #######################################################################################

    def snpposlist(self):
        if not self.snpposcache:
            snpposlist = list(self.snpdict.keys())
            snpposlist.sort()
            self.snpposcache = snpposlist
        return self.snpposcache

    #######################################################################################

    def consensus(self, pos):
        maxfreq = 0
        for nuc in "ACGT":
            if self.snpdict[pos][nuc]["freq"] > maxfreq:
                maxfreq = self.snpdict[pos][nuc]["freq"]
                connuc = nuc
        return(connuc)

    #######################################################################################

    def maf(self, pos):
        # Note: sum of frequencies for ALT nucs (or of single ALT if only one)
        refnuc = self.refdict[pos]
        majorfreq = self.snpdict[pos][refnuc]["freq"]
        maf = 1.0 - majorfreq
        return(maf)

    #######################################################################################

    def nucs(self, pos):
        nucset = set()
        for nuc in "ACGT":
            if self.snpdict[pos][nuc]["freq"] > 0:
                nucset.add(nuc)
        return(nucset)

    #######################################################################################

    def vcftxt(self, poslist=None, ascounts=False):
        # Returns string ready for printing of summary of nucleotide frequencies (or counts) for requested positions
        # Format (tab delimited):
        #   #   POS     A           C           G           T
        #       23      0.0413      0           0           0.9587
        #
        #   or:
        #
        #   #   POS     A           C           G           T
        #       23      230         0           0           5342

        outstringlist = ["{}\t{}\t{}\t{}\t{}\n".format("POS", "A", "C", "G", "T")]
        # If no specific positions are requested: return info for variable sites in VCF
        if not poslist:
            poslist = self.snpposlist()
        for pos in poslist:
            a, c, g, t = self.snpvector(pos, ascounts)
            outstringlist.append("{}\t{}\t{}\t{}\t{}\n".format(pos, a, c, g, t))
        outstring = "".join(outstringlist)
        return(outstring)

#############################################################################################
#############################################################################################


# Execute Test code (run module in standalone mode)
# def main():

#    start = time.time()

    # for i in range(100):
    #     seqfile = Fastafilehandle("/Users/gorm/Documents/python/example_data/ATS.aligned.fasta",seqtype="protein")
    #     seqs = seqfile.read_alignment()
    #     phy = seqs.phylip()
##
##    for i in range(1):
    # dmat = seqs.distmatrix("pdist_ignoregaps")
    # print dmat
        #physeqs = seqs.phylip()

#    print "Number of seqs: %d" % len(seqs)

#    print seqs[1]
##    print seqs.phylip()
##    print myalignment.tab()
##    print seqs[1].fasta()
##    print seqs[2].fasta()

    # seqfile = Seqfile("/Users/gorm/Documents/python/example_data/ATS.aligned.fasta")
    # seqs = seqfile.read_alignment()
    # for ali in range(1,61):
    #     print(ali, seqs.alignpos2seqpos("ATS_D1_ITseq19", ali, False))
    # print("#######")
    #
    # for seqi in range(1,48):
    #     print(seqi, seqs.seqpos2alignpos("ATS_D1_ITseq19", seqi, False))
    #

   # print seqs.distmatrix()
##
##    difs = gapgap = 0
##    for i in range(len(seqs[0])):
##        if seqs[0][i] != seqs[1][i]:
##            difs += 1.0
##            print "different: %s %s" % (seqs[0][i], seqs[1][i])
##        else:
##            print "same: %s %s" % (seqs[0][i], seqs[1][i])
##            if seqs[0][i] == "-":
##                gapgap +=1
##
##    print "Manual dist: %f" % (difs)
##    print "gapgap: %d  tot length: %d" % (gapgap, len(seqs[0]))

##    phylipfile = Phylipfilehandle("example-data/test.phy")
##    seqs = phylipfile.read_alignment()

##    seqfile = Fastafilehandle("example-data/small_alignment.fasta")
##    seqs = seqfile.read_alignment()
##    print seqs.fasta()
##    dmat = seqs.distmatrix("pdist")
##    for (name1,name2) in dmat:
##        print "%s %s  %.4f" % (name1, name2, dmat[(name1, name2)])


##    import psyco
##    psyco.bind(Aligner.fill_DP)

##    start = time.time()
##
##    mytab = Tabfilehandle("example-data/mhc.tab", "DNA", degap=True)
##    seq1 = mytab.readseq()
##    seq2 = mytab.readseq()

#    tim = time.time()-start
#    print "Time used: %.2f" % tim

# if __name__ == "__main__":

    # main()

    # import cProfile
    # cProfile.run('main()', 'mainprofile')
##
##    In separate window do:
##    >>> import pstats
##    >>> p = pstats.Stats('mainprofile')
##    >>> p.sort_stats('time').print_stats(20)
