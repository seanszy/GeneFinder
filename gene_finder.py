# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
dna_earliest = load_seq("./data/X73525.fa") #import DNA
print(dna_earliest)

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide
        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'

    This creates a test to make sure T is converted to A
    >>> get_complement('T')
    'A'

    This test makes sure that G is converted to C
    >>> get_complement('G')
    'C'
    """
    # TODO: implement this
    #switches every character, C to G, G to C, A to T, and T to A
    #Each if statement converts a letter to its pair.
    if nucleotide == "C":
        return("G") #return statements return the opposite dna pair
    if nucleotide == "G":
        return("C")
    if nucleotide == "A":
        return("T")
    if nucleotide == "T":
        return ("A")


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO: implement this
    dna_list = list(dna) #converts dna to a list
    dna_length = len(dna) #length integer to calculate how long to run while
    repetition = len(dna)
    reverse = "" #initialize blank string

    while (repetition > 0):
        dna_index = repetition-1 #takes index starting from end moving backwards
        current_letter = dna_list[dna_index] #the letter at the curren index
        reverse_letter = get_complement(current_letter) #finds the pair
        reverse = reverse + reverse_letter #add current letter to string
        repetition = repetition-1 #used as a count to end the loop

    return reverse


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'

    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'

    This function should return the whole string if there is no stop codon.
    The other tests don't include this.
    This test checks to make sure this is the case.
    >>> rest_of_ORF("ATGATATTCG")
    'ATGATATTCG'
    """
    # TODO: implement this
    length = len(dna) #number of characters in the dna string
    dna_list = list(dna)
    reach_end = False #boolea which is used to tell the loop to stop
    current_index = 0 #keeps track of where in the string it is
    while reach_end == False: #while loop used to search for end codons
        current_index = current_index + 3
        current_dna = dna_list [current_index-3:current_index] #this finds the current codon
        #loops exit if stop codon or end of list is reached
        if current_dna == ['T', 'A', 'G'] or current_dna == ['T','A','A'] or current_dna == ['T','G','A']:
            reach_end = True
        if current_index > len(dna):
            return (dna)
            reach_end = True

    rejoin_string = ''.join(dna_list[0:current_index-3]) #converts list of characters into string

    return rejoin_string

def find_all_ORFs_oneframe(dna):

    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']

    #This test makes sure that if function starts with something other than a start codon, it still works
    >>> find_all_ORFs_oneframe("CCCATGTAG")
    ['ATG']
    """
    # TODO: implement this
    dna_frames = [] #initialize variable that stores ORFs
    length = len(dna)
    dna_list = list(dna)
    start = False #boolean used to run inner while loop. Changes if the whole strand is searched
    current_index = -3
    full_dna_boolean = 1
    while full_dna_boolean < len(dna): #out while loop which finds multiple ORFs
        start = False
        input_boolean = False
        full_dna_boolean = full_dna_boolean + 1
        if len(dna_list[:]) == 0:
            full_dna_boolean = len(dna) + 3

        while start ==  False: #inner while loop takes an ORF and adds it dna_frames
            current_index = current_index+3 #moves one codon at a time
            current_dna = dna_list [current_index:current_index+3]
            if current_dna == ['A', 'T', 'G']: #if an ORF is starting it changes a boolean to send that string through rest_of_orf
                start = True
                input_boolean = True
            if current_index > len(dna_list[:])-1: #escapes if no more start codon
                return (dna_frames) #returns the ORFs
                full_dna_boolean = full_dna_boolean + 10
                start = True

        rejoin_string = ''.join(dna_list[current_index:])

        if input_boolean == True:
            newest_dna_frame = rest_of_ORF(rejoin_string) #adds an ORF
        else:
            break
        dna_list = dna_list[current_index+len(newest_dna_frame):] #calculates new dna string without current ORF
        dna_frames.append(newest_dna_frame) #adds the ORF to the frame
        current_index = -3

find_all_ORFs_oneframe("123ATGTAG123ATG456TAG123")

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']

    This test makes sure that if the dna doesn't start with a start codon, it
    still runs
    >>> find_all_ORFs("AAAATGCCCTAG")
    ['ATGCCC']
    """
    # TODO: implement this
    frame_one = find_all_ORFs_oneframe(dna)
    frame_two = find_all_ORFs_oneframe(dna[1:]) #finds frame shifted by 1 steps
    frame_three = find_all_ORFs_oneframe(dna[2:]) #finds frame shifted by 2 steps
    frame = frame_one+frame_two+frame_three
    return frame


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    original = find_all_ORFs(dna) #finds ORFS of original
    new = get_reverse_complement(dna)
    new = find_all_ORFs(new) #finds ORFs of switched
    together = original + new #combines them
    return together


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    dna_list = find_all_ORFs_both_strands(dna)
    maximum_index = max(dna_list, key=len) #finds the longest ORF
    return maximum_index


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    longest_ORF_shuffle = []
    for i in range(num_trials): #shuffles for num_trials
        dna = shuffle_string(dna)
        maximum_shuffled_dna = longest_ORF(dna) #finds maximum of the current stirng
        longest_ORF_shuffle = longest_ORF_shuffle + [maximum_shuffled_dna] #adds the maximum in the strand to a list of maximums in each strand

    maximum_index = max(longest_ORF_shuffle, key=len) #finds the index of the longest ORF overall
    return len(maximum_index)


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    count = 0
    amino_acids = [] #amino acid empty string
    while count < len(dna):
        count = count + 3
        current_codon = dna[count-3:count] #finds next codon
        #print(current_codon)
        if len(current_codon) == 3: #if the last codon is not a complete codon it won't break it
            amino_acid = aa_table[current_codon] #calculates amino acid from the table
            amino_acids = amino_acids + [amino_acid]

    rejoin_string = ''.join(amino_acids) #rejoins them
    return rejoin_string

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
        gene_finder("ATGAAACCCTTTGGGTAG")
        'AAAA'
    """
    # TODO: implement this
    threshold = longest_ORF_noncoding(dna, 1500) #find length of maximum ORF
    real_dna = []
    Aminos = find_all_ORFs_both_strands(dna) #convert DNA to aminos
    for ORF in Aminos:
        if len(ORF) > threshold:
            real_dna = real_dna + [coding_strand_to_AA(ORF)] #converts ORF to amino and adds to list
    print("threshold",threshold)
    print("hello",real_dna)
    return real_dna

gene_finder(dna_earliest)
if __name__ == "__main__":
    import doctest
    doctest.testmod()
