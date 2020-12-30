import os, sys
import argparse, logging
from typing import Dict, List
import secrets
import math

from pprint import pprint


logger = logging.getLogger("rnaConstructor")

parser = argparse.ArgumentParser(description="Optimize RNA sequence for target amino acids by mapping RNA codons")

parser.add_argument('-p', "--property-sequence", type=str, help="Biochemical property sequence")
parser.add_argument('-a', "--amino-sequence", type=str, help="Amino acid chain")
parser.add_argument('-s', "--supported-codons", type=str, default=None, help="List of supported codons for target")
parser.add_argument('-e', "--excluded-codons", type=str, default=None, help="List of codons to exclude for target")
parser.add_argument('-i', "--iterations", type=int, default=1, help="Number of iterations to attempt")
parser.add_argument('-d', "--duplicate", type=int,default=1, help="Number of times to duplicate the input sequence")



biochem_short_label =  {
    'NP':'Non-polar',
    'P': 'Polar',
    'B': 'Basic',
    'A': 'Acidic',
    'S': 'Stop'
}

biochem_amino_pairs = {
    'Non-polar':['Phe', 'Leu', 'Ile', 'Met', 'Val', 'Pro', 'Ala', 'Trp', 'Gly'],
    'Polar': ['Ser', 'Thr', 'Tyr', 'Gln', 'Asn', 'Cys'],
    'Basic': ['His', 'Lys', 'Arg'],
    'Acidic': ['Asp', 'Glu'],
    'Stop': ['Stop_Ochre', 'Stop_Amber', 'Stop_Opal']
}

amino_codon_pairs = {
    'Phe': ['UUU','UUC'],
    'Leu': ['UUA','UUG', 'CUU','CUC', 'CUA', 'CUG'],
    'Ile': ['AUU','AUC','AUA'],
    'Met': ['AUG'],
    'Val': ['GUU', 'GUC', 'GUA', 'GUG'],
    'Ser': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'Pro': ['CCU', 'CCC', 'CCA', 'CCG'],
    'Thr': ['ACU', 'ACC', 'ACA', 'ACG'],
    'Ala': ['GCU', 'GCC', 'GCA', 'GCG'],
    'Tyr': ['UAU', 'UAC'],
    'Stop_Ochre': ['UUA'],
    'Stop_Amber': ['UAG'],
    'His': ['CAU', 'CAC'],
    'Gln': ['CAA', 'CAG'],
    'Asn': ['AAU', 'AAC'],
    'Lys': ['AAA', 'AAG'],
    'Asp': ['GAU', 'GAC'],
    'Glu': ['GAA', 'GAG'],
    'Cys': ['UGU', 'UGC'],
    'Stop_Opal': ['UGA'],
    'Trp': ['UGG'],
    'Arg': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'Gly': ['GGU', 'GGC', 'GGA', 'GGG']
}

class RNAChainBuilder(object):
    def __init__(self):
        pass

    def longestRunLength(self, input_string: str) ->int:
        repeatedChar = ""
        longestRun = -1
        currentChar =""
        currentRun = 0
        for char in input_string:
            if char != currentChar:
                if currentRun > longestRun:
                    longestRun = currentRun
                    repeatedChar = currentChar
                currentRun = 0
                currentChar = char
            currentRun+=1
        return longestRun, repeatedChar


    def build_best_RNA_chain_from_property(self, property_list: List[str], useable_biochem_amino_pairs: Dict[str, List[str]],  useable_amino_codon_pairs: Dict[str, List[str]], iterations:int) -> str:
        best_result = self.build_RNA_chain_from_property(property_list, useable_biochem_amino_pairs,useable_amino_codon_pairs, iterations) 
        for _ in range(iterations-1):
            new_result = self.build_RNA_chain_from_property(property_list, useable_biochem_amino_pairs,useable_amino_codon_pairs, iterations)
            if new_result[1][1][0] < best_result[1][1][0]:
                best_result = new_result
        return best_result


    def build_RNA_chain_from_property(self, property_list: List[str], useable_biochem_amino_pairs: Dict[str, List[str]],  useable_amino_codon_pairs: Dict[str, List[str]], iterations:int) -> str:
        amino_chain = []
        for biochem_property in property_list:
            try:
                num_useable_amino = len(useable_biochem_amino_pairs[biochem_property])
            except KeyError:
                logger.error(f"Biochem property cannot be built from supported codons missing: {biochem_property}")
                return None
            if num_useable_amino > 1:
                amino_idx = secrets.randbelow(num_useable_amino)
            else:
                amino_idx = 0
            amino_chain.append(useable_biochem_amino_pairs[biochem_property][amino_idx])

        return ",".join(amino_chain), self.build_best_RNA_chain_from_amino(amino_chain, useable_amino_codon_pairs, iterations)


    def build_best_RNA_chain_from_amino(self, amino_list: List[str], useable_amino_codon_pairs: Dict[str, List[str]], iterations: int)-> str:
        best_result = self.build_RNA_chain_from_amino(amino_list, useable_amino_codon_pairs) 
        for _ in range(iterations-1):
            new_result = self.build_RNA_chain_from_amino(amino_list, useable_amino_codon_pairs)
            if new_result[1] < best_result[1]:
                best_result = new_result
        return best_result


    def build_RNA_chain_from_amino(self, amino_list: List[str], useable_amino_codon_pairs: Dict[str, List[str]]) -> str:
        rna_chain =""
        for amino in amino_list:
            try:
                num_useable_codons = len(useable_amino_codon_pairs[amino])
            except KeyError:
                logger.error(f"Amino acid cannot be built from supported codons missing: {amino}")
                return None
            if num_useable_codons > 1:
                codon_idx = secrets.randbelow(num_useable_codons)
            else:
                codon_idx = 0
            rna_chain+= useable_amino_codon_pairs[amino][codon_idx]

        # return rna_chain, self.longestRepeatedSubstring(rna_chain)
        return rna_chain, self.longestRunLength(rna_chain)


def get_supported_amino_codon_pairs(supported_codons: str, excluded_codons: str) -> Dict[str, List[str]]:
    if supported_codons is not None:
        supported_codons = supported_codons.strip("'").split(',')
    if excluded_codons is not None:
        excluded_codons = excluded_codons.strip("'").split(',') 
    
    useable_codons = {}
    for amino_acid, codon_list in amino_codon_pairs.items():
        for codon in codon_list:
            do_add = True
            if supported_codons is not None:
                do_add = False
                if codon in supported_codons:
                    do_add = True
            if excluded_codons is not None:
                if codon in excluded_codons:
                    do_add = False
            
            if do_add:
                try:
                    useable_codons[amino_acid].append(codon)
                except KeyError:
                    useable_codons[amino_acid] = [codon]

    if len(useable_codons) < len(amino_codon_pairs):
            logger.warn("Not all amino acids are buildable with the specified set of supported codons")

    return useable_codons


def get_supported_biochem_amino_pairs(useable_amino_codon_pairs: Dict[str, List[str]]) -> Dict[str, List[str]]:
    useable_biochem_amino_pairs = {}
    for biochem_property, amino_list in biochem_amino_pairs.items():
        for amino in amino_list:
            if amino in useable_amino_codon_pairs:
                try:
                    useable_biochem_amino_pairs[biochem_property].append(amino)
                except KeyError:
                    useable_biochem_amino_pairs[biochem_property] = [amino]

    if len(useable_biochem_amino_pairs) < len(biochem_amino_pairs):
            logger.warn("Not all biochem properties are buildable with the specified set of supported codons")

    return useable_biochem_amino_pairs



def build_amino_sequence(amino_sequence:str, supported_codons: Dict[str, List[str]], iterations: int):
    amino_list = amino_sequence.strip("'").split(',')
    
    builder = RNAChainBuilder()
    rna_chain = builder.build_best_RNA_chain_from_amino(amino_list, supported_codons, iterations)
    pprint(rna_chain)

def build_property_sequence(property_sequence:str, supported_biochem: Dict[str, List[str]], supported_codons: Dict[str, List[str]], iterations: int):
    property_list = [ biochem_short_label[k] for k in property_sequence.strip("'").split(',') if len(k) >0]
    
    builder = RNAChainBuilder()
    rna_chain = builder.build_best_RNA_chain_from_property(property_list, supported_biochem, supported_codons, iterations)
    pprint(rna_chain)



if __name__ == "__main__":
    try:
        args = parser.parse_args()
    except argparse.ArgumentError:
        logger.error("Cannot parse arguments")

    if args.supported_codons is None and args.excluded_codons is None:
        useable_amino_codon_pairs = amino_codon_pairs
    else:
        useable_amino_codon_pairs = get_supported_amino_codon_pairs(args.supported_codons, args.excluded_codons) 

    if args.amino_sequence is not None:
        if (args.duplicate > 1):
             args.amino_sequence += ("," + args.amino_sequence.strip("'")) * (args.duplicate -1)
        build_amino_sequence(args.amino_sequence, useable_amino_codon_pairs, args.iterations)
    
    if args.property_sequence is not None:
        useable_biochem_amino_pairs = get_supported_biochem_amino_pairs(useable_amino_codon_pairs)
        if (args.duplicate > 1):
            args.property_sequence = args.property_sequence.strip("'")
            args.property_sequence += ("," + args.property_sequence.strip("'")) * (args.duplicate -1)
        build_property_sequence(args.property_sequence, useable_biochem_amino_pairs, useable_amino_codon_pairs, args.iterations)
