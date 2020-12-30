# rnaConstructor

Brute force attempt to minimize the number of repeated RNA proteins.

Based on [RNA Codon tables](https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables).

I personally have no idea about the validity of the generated sequences.

## Example Usage

* An amino acid chain can be specified directly to be built.
* A biochemistry chain can be specified to be built, and the appropriate amino acids will be selected.
* A codon exclusion list can be specified to disable codons from being used
* A codon inclusion list can be specified to use only the specified codons


### Building an amino acid chain
`main.py -a='Phe,Leu,Ile,Ile,Ile,Ile' -s='UUU,CUU,AUG,AUA,AUC'`

This will try to built the specified amino acid chain, using the specified codons. 

### Building a biochemistry property chain
`main.py -i=100 -d=5, -p='A,NP' -e='UCG,CUC,CCG,UTU,GUC,GUG,GCG,GGG'`

This will try to build a sequence of "Acid, Non-Polar"x5, excluding the specified codons.
This also means that any amino acids that cannot be built via the exclusion of these codons, will be excluded.
If a desired biochemistry property cannot be built, a warning will be logged. If a required biochemistry property cannot be fufilled the system will fail.