#!/usr/local/bin/python
# Author: Eli Draizen
# Date: 16-3-2014
# File: PyFrame.py

#Standard Libraries
import argparse
import os, sys
from collections import defaultdict
import string
import random
import math

#Custom Libraries
from yahmm import *
from read_fasta import read_fasta, Sequence

"""Correct frameshift mutations using an HMM, soon this will be conbined
with the similarity-based correction program HSP-Tiler
"""

class HighOrderNucleotide(Distribution):
    """Higher order markov chains using the discrete distribution, made up of 
    nucleotides and their probabilities, assuming that these probabilities 
    will sum to 1.0. 
    """

    """
    This is the name that should be used for serializing this distribution.
    """
    name = "HighOrderNucleotideDistribution"

    def __init__(self, nucleotides=None, pseudocounts=True):
        """
        Make a new discrete distribution with a dictionary of discrete
        nucleotides and their probabilities, checking to see that these
        sum to 1.0. Each discrete nucletotide can be modelled as a
        Bernoulli distribution.

        Parameters:
        ___________
        nucleotides : dictionary with key sequence and value emission probability.
                      The condition will be the seq[:-1],
        pseudocounts : Bool. Add psuedocount if the nucleotide does not exist.
                       Default is False.
        """
        
        pseudocount = 1e-10 if pseudocounts else float("-inf")
        self.nucleotides = defaultdict(lambda: pseudocount)

        if nucleotides is  None:
            nucleotides = {"A": 0.25, "C": 0.25, "G": 0.25, "T":0.25}

        self.nucleotides.update(nucleotides.iteritems())
        self.order = len(self.nucleotides.keys()[0])-1
        self.parameters = [self.nucleotides, self.order]
    
    def log_probability(self, nucleotide):
        """
        What's the probability of the given nucleotide under this distribution?
        Simply the log probability value given at initiation. If the nucleotide
        is not part of the discrete distribution, return 0 or a pseudocount
        of .001.

        Parameters:
        ___________
        nucleotide : nucleotide to get prob for. Conditions are nucleotide[:-1]
                     to get the previous nucleotides for higher order Markov models.
        """
        return log(self.parameters[0][nucleotide])

    def sample(self):
        """
        Sample randomly from the discrete distribution, returning the character
        which was randomly generated.
        """
        
        rand = random.random()
        for key, value in self.parameters[0]:
            if value >= rand:
                return key[-1]
            rand -= value

HighOrderNucleotide.register()

def Stop(frame=None):
    """Build stop codon nucleotides. Probabilities based of Prof. Haussler's 
    toy HMM example of an E. coli genome

    Parameters:
    ___________
    frame : int. Frame the codon is in. frame=1,2,3

    Returns:
    ________
    nt1 - first nucleotide
    nt2 - second nucleotide
    nt3 - third nucleotide
    """
    name = "_f{}".format(frame) if frame is not None else ""

    nt1 = State(DiscreteDistribution(characters={"A":0.0, "C":0.0, "G":0.0, "T":1.0}), name="Stop_nt1{}".format(name))
    
    nt2Probs = {"AA": 0.0,  "AC":0.0,   "AG":0.0,  "AT":0.0,
                "CA": 0.0,  "CC":0.0,   "CG":0.0,  "CT":0.0,
                "TA": 2./3, "TC":0.0,   "TG":1./3,  "TT":0.0,
                "TA": 0.0,  "TC":0.0,   "TG":0.0,  "TT":0.0}

    nt2 = State(HighOrderNucleotide(nucleotides=nt2Probs), name="Stop_nt2{}".format(name))

    nt3Probs = {  "AA": 0.5,  "AC":0.0,   "AG":0.5,  "AT":0.0,
                  "CA": 0.0,  "CC":0.0,   "CG":0.0,  "CT":0.0,
                  "GA": 1.0,  "GC":0.0,   "GG":0.0,  "GT":0.0,
                  "TA": 0.0,  "TC":0.0,   "TG":0.0,  "TT":0.0,
                  "NA": 0.0,  "NC":0.0,   "NG":0.0,  "NT":0.0}
    
    nt3 = State(HighOrderNucleotide(nucleotides=nt3Probs), name="Stop_nt3{}".format(name))
    
    return nt1, nt2, nt3

def Codon(frame):
    """Make codon states and label from what frame they belong to

    Parameters:
    ___________
    frame : int. Frame the codon is in. frame=1,2,3

    Returns:
    ________
    nt1 - first nucleotide
    nt2 - second nucleotide
    nt3 - third nucleotide
    """
    nt1 = State(HighOrderNucleotide(), name="Codon_nt1_f{}".format(frame)) 

    nt2Probs = {"AA": 0.25, "AC":0.25, "AG":0.25, "AT":0.25,
                "CA": 0.25, "CC":0.25, "CG":0.25, "CT":0.25,
                "GA": 0.25, "GC":0.25, "GG":0.25, "GT":0.25,
                "TA": 0.25, "TC":0.25, "TG":0.25, "TT":0.25}

    nt2 = State(HighOrderNucleotide(nucleotides=nt2Probs), name="Codon_nt2_f{}".format(frame))

    nt3Probs = {"AAT": 0.25, "AAG": 0.25, "AAC": 0.25, "AAA": 0.25,
                "ACT": 0.25, "ACG": 0.25, "ACC": 0.25, "ACA": 0.25,
                "AGT": 0.25, "AGG": 0.25, "AGC": 0.25, "AGA": 0.25,
                "ATT": 0.25, "ATG": 0.25, "ATC": 0.25, "ATA": 0.25,
                "CAT": 0.25, "CAG": 0.25, "CAC": 0.25, "CAA": 0.25,
                "CCT": 0.25, "CCG": 0.25, "CCC": 0.25, "CCA": 0.25,
                "CGT": 0.25, "CGG": 0.25, "CGC": 0.25, "CGA": 0.25,
                "CTT": 0.25, "CTG": 0.25, "CTC": 0.25, "CTA": 0.25,
                "GAT": 0.25, "GAG": 0.25, "CAC": 0.25, "CAA": 0.25,
                "GCT": 0.25, "GCG": 0.25, "GCC": 0.25, "GCA": 0.25,
                "GGT": 0.25, "GGG": 0.25, "GGC": 0.25, "GGA": 0.25,
                "GAT": 0.25, "GAG": 0.25, "GAC": 0.25, "GAA": 0.25,
                "GTT": 0.25, "GTG": 0.25, "GTC": 0.25, "GTA": 0.25,
                "TAT": 0.50, "TAG": 0.00, "TAC": 0.50, "TAA": 0.00,
                "TCT": 0.25, "TCG": 0.25, "TCC": 0.25, "TCA": 0.25,
                "TGT": 1./3, "TGG": 1./3, "TGC": 1./3, "TGA": 0.00,
                "TTT": 0.25, "TTG": 0.25, "TTC": 0.25, "TTA": 0.25}
    nt3 = State(HighOrderNucleotide(nucleotides=nt3Probs), name="Codon_nt3_f{}".format(frame))

    return nt1, nt2, nt3

def fill3(frame):
    """Make fill states and label from what frame they belong to

    Parameters:
    ___________
    frame : int. Frame the codon is in. frame=1,2,3

    Returns:
    ________
    nt1 - first nucleotide
    nt2 - second nucleotide
    nt3 - third nucleotide
    """
    nt1 = State(HighOrderNucleotide(), name="Fill_nt1_Frame{}".format(frame))
    nt2 = State(HighOrderNucleotide(), name="Fill_nt2_Frame{}".format(frame))
    nt3 = State(HighOrderNucleotide(), name="Fill_nt3_Frame{}".format(frame))
    return nt1, nt1, nt3

def DraizenModel():
    draizen = Model(name="Draizen")

    #global spacer, random nts before frame stop codons is seen
    global_spacer = State(HighOrderNucleotide(), name="GlobalSpacer")
    draizen.add_state(global_spacer)

    #Fist stop codons seen in all frames
    stop_nt1_f1, stop_nt2_f1, stop_nt3_f1 = Stop(frame=1)
    draizen.add_state(stop_nt1_f1)
    draizen.add_state(stop_nt2_f1)
    draizen.add_state(stop_nt2_f1)

    stop_nt1_f2, stop_nt2_f2, stop_nt3_f2 = Stop(frame=2)
    draizen.add_state(stop_nt1_f2)
    draizen.add_state(stop_nt2_f2)
    draizen.add_state(stop_nt2_f2)

    stop_nt1_f3, stop_nt2_f3, stop_nt3_f3 = Stop(frame=3)
    draizen.add_state(stop_nt1_f3)
    draizen.add_state(stop_nt2_f3)
    draizen.add_state(stop_nt2_f3)

    #Spacers before start codons in each frame
    spacer2_1 = State(HighOrderNucleotide(), name="Spacer1_f2")
    draizen.add_state(spacer2_1)
    spacer3_1 = State(HighOrderNucleotide(), name="Spacer1_f3")
    draizen.add_state(spacer3_1)
    spacer3_2 = State(HighOrderNucleotide(), name="Spacer2_f3")
    draizen.add_state(spacer3_2)

    #Add each codon from each frame
    nt1_f1, nt2_f1, nt3_f1 = Codon(frame=1)
    draizen.add_state(nt1_f1)
    draizen.add_state(nt2_f1)
    draizen.add_state(nt3_f1)

    nt1_f2, nt2_f2, nt3_f2 = Codon(frame=2)
    draizen.add_state(nt1_f2)
    draizen.add_state(nt2_f2)
    draizen.add_state(nt3_f2)

    nt1_f3, nt2_f3, nt3_f3 = Codon(frame=3)
    draizen.add_state(nt1_f3)
    draizen.add_state(nt2_f3)
    draizen.add_state(nt3_f3)

    #Final stop codons seen in all frames
    stop_end_1, stop_end_2, stop_end_3 = Stop()
    draizen.add_state(stop_end_1)
    draizen.add_state(stop_end_2)
    draizen.add_state(stop_end_3)

    #gloabal spacer, random nts before frame stop codons is seen
    global_spacer2 = State(HighOrderNucleotide(), name="GlobalSpacer2")
    draizen.add_state(global_spacer2)

    #Start
    draizen.add_transition(draizen.start, global_spacer, 0.1)
    draizen.add_transition(draizen.start, stop_nt1_f1, 0.3)
    draizen.add_transition(draizen.start, spacer2_1, 0.3)
    draizen.add_transition(draizen.start, spacer3_1, 0.3)
    #draizen.add_transition(draizen.start, nt1_f1, .999)


    #Global Spacer
    draizen.add_transition(global_spacer, global_spacer, 0.01)
    draizen.add_transition(global_spacer, stop_nt1_f1, 0.33)
    draizen.add_transition(global_spacer, spacer2_1, 0.33)
    draizen.add_transition(global_spacer, spacer3_1, 0.33)

    #Frame 1 stop
    draizen.add_transition(stop_nt1_f1, stop_nt2_f1, 1.)
    draizen.add_transition(stop_nt2_f1, stop_nt3_f1, 1.)
    draizen.add_transition(stop_nt3_f1, nt1_f1, 1.)

    #Frame 1
    draizen.add_transition(nt1_f1, nt2_f1, 98./100)
    draizen.add_transition(nt2_f1, nt3_f1, 98./100)
    draizen.add_transition(nt3_f1, nt1_f1, 1-(2./100+1./500))
    draizen.add_transition(nt3_f1, stop_end_1, 1./500)

    #Frame 2 space and fill
    draizen.add_transition(spacer2_1, stop_nt1_f2, 1.0)
    #draizen.add_transition(spacer2_1, nt1_f2, 0.999)

    #Frame 2 stop
    draizen.add_transition(stop_nt1_f2, stop_nt2_f2, 1.0)
    draizen.add_transition(stop_nt2_f2, stop_nt3_f2, 1.0)
    draizen.add_transition(stop_nt3_f2, nt1_f2, 1.0)

    #Frame 2
    draizen.add_transition(nt1_f2, nt2_f2, 98./100)
    draizen.add_transition(nt2_f2, nt3_f2, 98./100)
    draizen.add_transition(nt3_f2, nt1_f2, 1-(2./100+1./500))
    draizen.add_transition(nt3_f2, stop_end_1, 1./500)

    #Frame 3 space and fill
    draizen.add_transition(spacer3_1, spacer3_2, 1.0)
    draizen.add_transition(spacer3_2, stop_nt1_f3, 1.0)
    #draizen.add_transition(spacer3_2, nt1_f3, 0.999)

    #Frame 3 stop
    draizen.add_transition(stop_nt1_f3, stop_nt2_f3, 1.0)
    draizen.add_transition(stop_nt2_f3, stop_nt3_f3, 1.0)
    draizen.add_transition(stop_nt3_f3, nt1_f3, 1.0)

    #Frame 3
    draizen.add_transition(nt1_f3, nt2_f3, 98./100)
    draizen.add_transition(nt2_f3, nt3_f3, 98./100)
    draizen.add_transition(nt3_f3, nt1_f3, 1-(2./100+1./500))
    draizen.add_transition(nt3_f3, stop_end_1, 1./500)

    #End
    draizen.add_transition(stop_end_1, stop_end_2, 1.)
    draizen.add_transition(stop_end_2, stop_end_3, 1.)
    draizen.add_transition(stop_end_3, global_spacer2, 0.5)
    draizen.add_transition(stop_end_3, draizen.end, 0.5)
    draizen.add_transition(global_spacer2, draizen.end, 0.5)
    draizen.add_transition(global_spacer2, global_spacer2, 0.5)

    #Find frameshifts:

    #Insertions: frame i to frame i+2
    draizen.add_transition(nt1_f1, nt1_f3, 1./100)
    draizen.add_transition(nt2_f1, nt2_f3, 1./100)
    draizen.add_transition(nt3_f1, nt3_f3, 1./100)
    draizen.add_transition(nt1_f2, nt2_f1, 1./100)
    draizen.add_transition(nt2_f2, nt2_f1, 1./100)
    draizen.add_transition(nt3_f2, nt3_f1, 1./100)
    draizen.add_transition(nt1_f3, nt1_f2, 1./100)
    draizen.add_transition(nt2_f3, nt2_f2, 1./100)
    draizen.add_transition(nt3_f3, nt3_f2, 1./100)

    #Deletions: frame i to frame i+1
    draizen.add_transition(nt1_f1, nt3_f2, 1./100)
    draizen.add_transition(nt2_f1, nt1_f2, 1./100)
    draizen.add_transition(nt3_f1, nt2_f2, 1./100)
    draizen.add_transition(nt1_f2, nt3_f3, 1./100)
    draizen.add_transition(nt2_f2, nt1_f3, 1./100)
    draizen.add_transition(nt3_f2, nt2_f3, 1./100)
    draizen.add_transition(nt1_f3, nt3_f1, 1./100)
    draizen.add_transition(nt2_f3, nt1_f1, 1./100)
    draizen.add_transition(nt3_f3, nt2_f1, 1./100)

    draizen.bake()

    return draizen

def updateScore(score, state1, state2):
    #return (score + confidence*)
    pass

def update_model_transitions(model, tile):
    """Update the model's transition probabilities"""
    pass

def abInitio(fasta, outfile, correct=True, annotation=None):
    """Find and correct frameshifts is fasta sequences.

    Parameters:
    ___________
    fasta : file-like object contiang FASTA sequences
    outfile : file-like object to write corrected sequences
    correct : Bool. Identify and correct frameshifts, else only identify
              frameshift mutations. Default is True.
    annotation : file-like object to BLASTX output for similairty based 
                 error correction. Optional.
    """
    #Build HMM
    model = DraizenModel()
    best_seq = None
    best_prob = float("-inf")

    for fasta in read_fasta(fasta):
        #Only test with Contig685 from HSP-Tiler v1.0, hacky will be removed
        if not fasta.name == "Contig4721": continue 

        #Get Vitertbi Path for forward stand sequence
        sequence = fasta.sequence
        prob, path = model.viterbi(fasta.sequence)

        #Get Viterbi path got the reverse complement
        rc = revcomp(fasta.sequence)
        prob_rc, path_rc = model.viterbi(rc)


        #Choose which ever has a larger score
        if prob_rc > prob:
            path = path_rc
            prob = prob_rc
            sequence = rc
            print "Using Reverse Complement"

        print prob
        print sequence

        #Updated sequence with out frame shift mutations
        corrected_sequence = Sequence()
        print sequence
        corrected_sequence.name = fasta.name
        
        #Control when sequence begens after first stop codon
        read_sequence = False

        #frame that corrected sequence started in
        frame = 1

        #index of the nt that started the sequence
        start = 0

        #Save preivous state to find framehsifts
        prevState = ""

        for index, state in path[1:]:
            seqChar = sequence[index-1]
            print index, seqChar, state.name, prevState
            if not read_sequence and (state.name.startswith("Stop_nt1") or state.name.startswith("Codon")):
                #Start the sequence
                read_sequence = True
                frame = int(state.name[-1])
                start = index-1

            elif len(corrected_sequence) > 3 and prevState.startswith("Stop_nt3"):
                #Stop the sequence
                read_sequence = False
                break

            if read_sequence:
                print "Reading Sequence"
                if prevState[-1] == state.name[-1] or state.name.startswith("Stop"):
                    #Same frame
                    print "Same Frame"
                    corrected_sequence.sequence += seqChar

                elif ((int(prevState[-1]) == 1 and int(state.name[-1]) == 3) or
                      (int(prevState[-1]) == 2 and int(state.name[-1]) == 1) or
                      (int(prevState[-1]) == 3 and int(state.name[-1]) == 2)):
                    #Insertions jump from frame i to frame i+2
                    #Must delete nt to get back to original frame
                    print "Insertion", state.name[-1], prevState[-1]
                    pass

                elif ((int(prevState[-1]) == 1 and int(state.name[-1]) == 2) or
                      (int(prevState[-1]) == 2 and int(state.name[-1]) == 3) or
                      (int(prevState[-1]) == 3 and int(state.name[-1]) == 1)):
                    #Deletion jump from frame i to frame i+1
                    #Must insert nt to get back to correct frame
                    print "Deletion"
                    corrected_sequence.sequence += random.choice("ACGT") #Fix to account for codon bias?
                    corrected_sequence.sequence += seqChar
                else:
                    raise RuntimeError("Problem with model")

            prevState = state.name

        new_prob, new_path = model.viterbi(corrected_sequence.sequence)
        if new_prob-prob > best_prob:
            best_prob = new_prob-prob
            best_seq = (fasta, corrected_sequence)

        for index, state in path[1:]:
            print index, corrected_sequence.sequence[index-1], state.name


        corrected_sequence.description = "[Old Score={}; New Score={}; start={}; frame={}]".format(math.exp(prob),
                                                                                                  math.exp(new_prob),
                                                                                                  start,
                                                                                                  frame)
        print >> outfile, corrected_sequence

    print best_seq[0]
    print best_seq[1]

        



        #fig, ax = plt.subplots(3)
        #for i in xrange(1,4):
            #for start, end in states[i-1]:
                #print start, end
                #plt.plot((int(start), int(end)), (5,5))
            #print
        #plt.show()

complement_table = string.maketrans("ACGT", "TGCA")
def revcomp(sequence):
    """returns the reverse complement of a sequence. Adapted from Kevin Karplus'
    BME 205 assignment 1 at UCSC.

    Input:
    sequence - string containing nucleotide sequence
    """
    return sequence[::-1].translate(complement_table) 


def parse_args():
    """Parsing command line options
    """
    parser = argparse.ArgumentParser(description="Takes a fasta file of sequences and a BLASTX annotation of that file in xml format.  Attempts to tile Hsps for the highest scoring hit for each sequence, correcting frameshifts in order to improve subsequent annotations.")
    # name of fasta file 
    parser.add_argument("-f", "--fasta", 
                        required=True, 
                        type=argparse.FileType('r'),
                        help="Fasta file containing sequences")
    # name of annotation file (in xml format as code currently stands)
    parser.add_argument("-a", "--annotation", 
                        required=False,
                        default=None, 
                        type=argparse.FileType('r'),
                        help="Blastx xml file containing annotations for sequences")
    # gap limit
    parser.add_argument("-g", "--gap_limit", 
                        type=int, 
                        default=15, 
                        help="Cutoff for distance between hsps. If the gap between hsps in nucleotides is greater than this value, the hsp will not be added to the tile.  Default = 15nt")
    #evlaue cutoff
    parser.add_argument("-e", "--evalue_cutoff",
                        type=float,
                        default=1e-10,
                        help="Only allow Blast Hist less than or equal to cutoff. Default is 1e-10")
    #Use every hit for eqch query instead of just the 1st one
    parser.add_argument("--allHits",
                        default=False,
                        help="Use all hits from annotation. Default is to use only the first, high scoring hit. Optional.",
                        action="store_true")
 
    #Define output
    parser.add_argument("-o", "--outfile",
                        type=argparse.FileType('wt'),
                        default=sys.stdout,
                        help="File to save corrected sequences")
    parser.add_argument("-l", "--logfile",
                        type=argparse.FileType('wt'),
                        default=sys.stderr,
                        help="File to save log")

    #Parse args
    return parser.parse_args()

if __name__ == "__main__":
    #Parse args
    args = parse_args()

    abInitio(args.fasta, args.outfile)


 





