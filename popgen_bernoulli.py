#! /usr/bin/env python3

"""
A module that calculates the probability of finding a number of observed traits given the sample size and frequency  within a population
"""

from Bio import Seq
from Bio import SeqIO
import sys

class Individual:

    """
    Individual is a class that contains attributes for individuals in a FASTA file

    Attributes
    ----------
    determine_phenotype
        Function that determines the phenotype (orange vs blue) based on the fourth nucleotide in the sequence

    translate_seq
        Function that takes a DNA nucleotide sequence and translates it into an amino acid sequence
    """

    def __init__(self, sample_id, seq_date, nuc_seq):

        """
        Initializes Individual object

        Parameters
        ----------
        sample_id: obj
            Identification description for the individual sample

        seq_date: obj
            Date on which sample was collected/uploaded

        nuc_seq: obj
            Sequence of nucleotides from the sample

        Returns
        -------
        None
        """

        self.id = sample_id
        self.seq_date = seq_date
        self.nuc_seq = nuc_seq
        self.aa_seq = self.translate_seq()
        self.phenotype = self.determine_phenotype()

    def translate_seq(self):

        """
        Translates a given DNA nucleotide sequence into an amino acid sequence

        Parameters
        ----------
        None

        Returns
        -------
        aa_seq: obj
            Translated sequence of amino acids for sample
        """

        aa_seq = self.nuc_seq.translate()
        return str(aa_seq)


    def determine_phenotype(self):

        """
        Identifies the phenotype of a seqeunce as either orange or blue based on the the fourth amino acid

        Parameters
        ----------
        None

        Returns
        -------
        phenotype: str
            Either "orange" or "blue" base on the fourth amino acid
        """

        phenotype = ""
        if self.aa_seq[3] == "R":
            phenotype = "orange"
        else:
            phenotype = "blue"
        return phenotype


def read_fasta(filepath):

    """
    Reads and parses through a FASTA file for sequence information into Individual objects

    Parameters
    ----------
    filepath: str
        Filepath for input file that will be read and parsed

    Returns
    -------
    individuals: list
        List of Individual objects, each represnting an individual sample in the FASTA file
    """

    individuals = []
    for record in SeqIO.parse(filepath, "fasta"):
        header = record.description.split("_")
        sample_id = header[0]
        seq_date = header[1]
        nuc_seq = record.seq
        individual = Individual(sample_id, seq_date, nuc_seq)
        individuals.append(individual)
    return individuals

def total_sample_set(individuals):

    """
    Finds the sample size (total number of samples in a FASTA file)

    Parameters
    ----------
    individuals: list
        List of Individual objects

    Returns
    -------
    n: int
        Total sample size
    """

    n = len(individuals)
    return n

def total_orange_phenotype(individuals):

    """
    Determines the number of samples with the "orange" phenotype

    Parameters
    ----------
    individuals: list
        List of Individual objects

    Returns
    -------
    k: int
        Total number of individuals with the orange phenotype
    """

    k = 0
    for individual in individuals:
        if individual.phenotype == "orange":
            k += 1
    return k

def get_factorial(n):

    """
    Calculates a factorial

    Parameters
    ----------
    n: int
        Input value for factorial (n!)

    Returns
    -------
    1: int
        If n is 1, return 1 (indicates end of factorial)

    n * get_factorial(n-1): int
        If n is not 1, multiply it by the factorial of n-1 (indicates continuation of factorial)
    """

    if n == 1:
        return 1
    else:
        return n * get_factorial(n-1)

def calculate_bernoulli(n, k, p):
    """
    Calculates the Bernoulli trial binomial distribution

    Parameters
    ----------
    n: int
        Total sample size

    k: int
        Total number of individuals with the orange phenotype

    p: float
        Frequency of orange within a population

    Returns
    -------
    combination * powers: int
        Probability of finiding a number of observed orange traits given the sample size and frequency of orange within the population
    """

    combination = get_factorial(n) / get_factorial(n-k) / get_factorial(k)
    powers = pow(p, k) * pow(1-p, n-k)
    return combination * powers

def create_output_file(filepath, n, k, p, bernoulli):

    """
    Creates output result file with necessary information regarding the Bernoulli Trial

    Parameters
    ----------
    filepath: str
        Filename for output file

    n: int
        Total sample size

    k: int
        Total number of individuals with the orange phenotype

    p: float
        Frequency of orange within a population

    bernoulli: float
        Probability of finding a number of observed orange traits given the sample size and frequency of orange within the population

    Returns
    -------
    None
    """

    output_string = (
        f"Results\n"
        f"p (the frequency of \"orange\" in the population) = {p}\n"
        f"n (the number of sample individuals) = {n}\n"
        f"k (the number of \"orange\" individuals in the sample set) = {k}\n\n"
        f"Probability of collecting {n} individuals with {k} being \"orange\" (given a population frequency of {p} = {bernoulli}\n"
    )
    with open(filepath, "w") as outfile:
        outfile.write(output_string)
    print("Output file written.")

def main():
    if len(sys.argv) != 4:
        sys.exit(sys.argv[0] + ": Expecting three arguments (input_file, probability, output_file)")
    input_file = sys.argv[1]
    p = float(sys.argv[2])
    output_file = sys.argv[3]
    individuals = read_fasta(input_file)
    n = total_sample_set(individuals)
    k = total_orange_phenotype(individuals)
    bernoulli = calculate_bernoulli(n, k, p)
    create_output_file(output_file, n, k, p, bernoulli)

if __name__=="__main__":
    main()
