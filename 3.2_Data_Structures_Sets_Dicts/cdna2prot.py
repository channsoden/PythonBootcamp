'''
This script will convert cDNA sequence to protein
Use as follows:
python <cdna fasta file name> <amino acid fasta file name>

a sample fasta file named sample_nucl.fa can be used
'''
import sys

def codon_table():
    '''
    This function will create the codon table
    Returns a dictionary of codons and their corresponding AA
    '''
    bases      = ['T', 'C', 'A', 'G']
    codons     = [a+b+c for a in bases for b in bases for c in bases]
    aminoAcids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codonTable = dict(zip(codons, aminoAcids))

    return codonTable


def fasta(inFile):
    '''
    This function reads in a fasta file.
    Returns a dictionary of sequence name and sequences.
    If fasta is multiline fasta, sequences will be fragmented.
    This will be fixed later on in the code.
    '''
    fastaDict = {} # initialize dictionary
    
    # loop through file
    with open(inFile) as i:
        for line in i:
            line = line.strip()
            
            if not line: # empty line
                continue
            
            if line.startswith('>'): # name of sequence
                seqName = line[1:]
                if seqName not in fastaDict:
                    fastaDict[seqName] = []
                continue
            
            seq = line
            fastaDict[seqName].append(seq)

    return fastaDict


def convert_AA(sequence):
    # initialize codon table and amino acid sequence variable
    codonTable = codon_table()
    aaSeq      = []

    #loop through sequence
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]

        if codonTable[codon] != '*': # if the codon is not a stop codon, add AA
            aaSeq.append(codonTable[codon])
        else: # if codon is a stop codon, stop the loop--we'll only sequence up to this point
            break

    return ''.join(aaSeq)


def main():
    outFile = open(sys.argv[2], 'w') # open output file
    fastaDict = fasta(sys.argv[1])   # open 
    
    for gene in fastaDict:
        fullSequence = ''.join(fastaDict[gene]) # in case we have multi-line fasta, we combine all lines in a list
        peptide      = convert_AA(fullSequence) # convert sequence to amino acid sequence

        outFile.write('>{}\n{}\n' .format(gene, peptide))

if __name__ == '__main__':
    main()
