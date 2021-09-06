#!/usr/bin/env python

'''

   Author: Sofia Papadimitriou
   
   Script that takes a TAB file with variant information and a fasta file
   with protein sequences and appends amino acid property differences
   between the reference and the alternative protein alleles.

   == INPUT ==: 1. a TAB-delimited variants file with one variant per line and columns with information in the
   following order: translation start position, reference and alternative protein alleles and the gene name
   (or other identifier present in the fasta file, see below).
   More columns can be available, but only the first 4 will be used in the script.

                2. a .fasta file that contain the protein sequences of the genes involved and the identifiers present
   in the variants file.


   == PROCESS ==: For SNPs it uses the difference between the aminoacids in the mutation
   position. For insertions/deletions, it uses the protein sequence FASTA file to get a
   range of 5 aminoacids prior and after the position to calculate the 
   difference, along with the deleted or inserted amino acids.


   == OUTPUT == : aadiff.txt, containing the input file with appended at the end the
   Flexibility and Hydrophobicity amino acid differences as extra columns.

   
   Usage examples 
   ===========================================================================
   
   To run the script with the example variant and fasta files:
   Usage: $ python calc_aa_diff.py
   
   To provide a prefix at the output file, use option -f: 
   Usage: $ python calc_aa_diff.py -f test

 
'''

import argparse
import tarfile
import os.path

def input_parser():
 
    """ Function to construct and parse input parameters """
 
    parser = argparse.ArgumentParser(usage=__doc__)    
    parser.add_argument('-v','--variants_file', 
                        help = 'Input variants file (DEFAULT:./example_data/example_variants.txt',
                        default='./example_data/example_variants.txt')
    parser.add_argument('-s','--fasta_file', 
                        help = 'Fasta file with gene/Uniprot identifiers and protein sequences (DEFAULT: ./example_data/example_sequences.fasta)',
                        default='./example_data/example_sequences.fasta'
                               )
    parser.add_argument('-f','--filename',
                        help = 'Prefix name of the output file',
                        default='')
                        
    return vars(parser.parse_args())                    


def open_file(filename):
    
    '''     
    Function that takes the full name of a file and opens it for reading or 
    prints an error otherwise.    
    '''

 
    if filename.endswith('.tar') or filename.endswith('.gz'):
        files=tarfile.open(filename, 'r') 
        for member in files:       
            input_file=files.extractfile(member)
    else:
        input_file=open(filename)
        
        
    return input_file
    

def parse_fasta(filename):
    
    ''' Function that takes a file, opens it, parses fasta sequences and returns
    a dictionary with the header as key and the sequence as value'''
    
    #initiate dictionary and first key value
    prot_seqs={}
    key=''
    
    fasta_file=open(filename)
    
    #iterate the file
    for line in fasta_file:
        if line.startswith('>'):            
            prot_seqs[line.rstrip('\n').rstrip('\r').lstrip('>')]=''
            key=line.rstrip('\n').rstrip('\r').lstrip('>')
        else:
            prot_seqs[key]+=line.rstrip('\n').rstrip('\r')
    
    fasta_file.close()
    
    return prot_seqs      


def diffwwH (wt,mt):
    
    #  Amino acid hydrophobicity scale from Wimley and White, Nat Struct Biol 3:842 (1996).
    wwH={"D":-1.23, "E":-2.02,"N":-0.42,"Q":-0.58,"K":-0.99,"R":-0.81}
    wwH.update({"H":-0.96,"G":-0.01,"P":-0.45,"S":-0.13,"T":-0.14,"C":0.24})
    wwH.update({"M":0.23,"A":-0.17,"V":-0.07,"I":0.31,	"L":0.56,"F":1.13})
    wwH.update({"W":1.85,"Y":0.94})
      
    wwHwt=0
    wwHmt=0
    
    #calculate value for each aminoacid
    
    for aa in wt:
        try:
            wwHwt+=wwH[aa]
        except KeyError:    
            
            return 'NaN'
            
    for aa in mt:
        try:
            wwHmt+=wwH[aa]
        except KeyError:
            
            return 'NaN'
    
    #calculate means
    mean_wwHwt=wwHwt/len(wt) 
    mean_wwHmt=wwHmt/len(mt)
    
    return mean_wwHwt-mean_wwHmt


#average flexibility index
def diff_flex(wt,mt):
    
    flex={"A":0.360, "C":0.350,"D":0.510,"E":0.500,"F":0.310,"G":0.540,}
    flex.update({"H":0.320,"I":0.460,"K":0.470,"L":0.370,"M":0.300,"N":0.460})
    flex.update({"P":0.510, "Q":0.490,"R":0.530,"S":0.510,"T":0.440,"V":0.390})
    flex.update({"W":0.310,"Y":0.420})
    
    flexwt=0
    flexmt=0
    
    #calculate value for each aminoacid
    for aa in wt:
        try:
            flexwt+=flex[aa]
        except KeyError:
            return 'NaN'
    for aa in mt:
        try:
            flexmt+=flex[aa]
        except KeyError:
            return 'NaN'
    
    #calculate means
    mean_flexwt=flexwt/len(wt) 
    mean_flexmt=flexmt/len(mt)
    
    return mean_flexwt-mean_flexmt
                       
                                                                     
if __name__=='__main__':
    
    #Calling the parser function to start parsing the argument
    parser = input_parser()    
    
    # check if input files exist
    if not os.path.exists(parser['variants_file']):    
        print('##### Input file %s or directory path does not exist.'%\
        (parser['variants_file']))
        print('##### Script is terminated')
        quit()
    
    if not os.path.exists(parser['fasta_file']):    
        print('##### Input file %s or directory path does not exist.'%\
        (parser['fasta_file']))
        print('##### Script is terminated')
        quit()
        
    ######################################## Parse input files 
    print('##### Parsing variants file')
    
    variants_file=open_file(parser['variants_file'])
    input_header=variants_file.readline().rstrip('\n').rstrip('\r') 

    print('##### Parsing fasta file')
    
    prot_dict=parse_fasta(parser['fasta_file']) 
    
    
    ######################################## Prepare the output file 
    
    if parser['filename']=='':
        output_file=open('aadiff.txt','w')
    else:
        output_file=open(parser['filename']+'_aadiff.txt','w')
        
    output_file.write(input_header+'\tFlex_diff\tHydr_diff\n')
    
    
    ######################################## Calculate aa differences 
    print('##### Appending amino acid property differences')

    #iterate over the input file     
    for line in variants_file:     
        
        line_list=line.rstrip('\n').rstrip('\r').split('\t')
        
        #check position,reference,alternative
        pos=line_list[0]
        ref=line_list[1]
        alt=line_list[2]
        uni=line_list[3]
      
        ############# if missense SNV
        if (alt not in ['*','N/A','-','X','*']) and\
        (ref not in ['*','-','X','*']) and len(alt)<2 and len(ref)<2:            
            
            h_diff=diffwwH(ref,alt)           
            fl_diff=diff_flex(ref,alt)            
      
        ############# if nonsense mutation
        elif alt=='*' and ref not in ['N/A','-','X','*']:  
           
            pos=int(pos)-1 #substracting 1 as string starts with zero index in python            
            prot_seq=''
            ref_seq=''
            alt_seq=''
            
            #take the range of the closest aminoacids and the reference
            for key in prot_dict:
                if uni==key:   
                    
                    prot_seq=prot_dict[key]
                    
                    #if the position is inside the protein length                    
                    if pos<len(prot_seq):
                        
                        
                        if (pos-5)>=0 and (pos+6)<=(len(prot_seq)-1):                            
                            ref_seq=prot_seq[pos-5:pos]+prot_seq[pos]+prot_seq[pos+1:pos+6]
                            alt_seq=prot_seq[pos-5:pos] #stop counting after the position
                        elif (pos-5)<0:
                            ref_seq=prot_seq[0:pos]+prot_seq[pos]+prot_seq[pos+1:pos+6]
                            alt_seq=prot_seq[0:pos]
                        elif (pos+6)>(len(prot_seq)-1): #stop counting after the position
                            ref_seq=prot_seq[pos-5:pos]+prot_seq[pos]+prot_seq[pos+1:len(prot_seq)+1]
                            alt_seq=prot_seq[pos-5:pos] #stop counting after the position
                          
                            
                        h_diff=diffwwH(ref_seq,alt_seq)                        
                        fl_diff=diff_flex(ref_seq,alt_seq)     
                    
                    #else, consider it an error and append 'N/A's
                    else:                       
                        h_diff='NaN'
                        fl_diff='NaN'
                        
                   
        ############# if stop loss mutation
        elif ref=='*' and alt not in ['N/A','-','X','*']:   
      
            h_diff='NaN'
            fl_diff='NaN'
    
 
        ############# if non synonymous stop mutation
        elif alt=='*' and ref=='*':        
           
            h_diff=0.0           
            fl_diff=0.0 
  
          
        ############# if insertion only
        elif ref=='-' and alt not in ['N/A','X','*','-'] and 'X' not in alt:  
            pos=int(pos)-1 #substracting 1 as string starts with zero index in python              
            prot_seq=''
            ref_seq=''
            alt_seq=''
                        
            #calculate the scale differences
            for key in prot_dict:
                if uni==key:
                    prot_seq=prot_dict[key]
                    #if the position is inside the protein length
                    if pos<len(prot_seq):  
                        if (pos-5)>=0 and (pos+6)<=(len(prot_seq)-1): 
                            ref_seq=prot_seq[pos-5:pos]+prot_seq[pos]+prot_seq[pos+1:pos+6]
                            alt_seq=prot_seq[pos-5:pos]+prot_seq[pos]+alt+prot_seq[pos+1:pos+6]
                        elif (pos-5)<0:
                            ref_seq=prot_seq[0:pos]+prot_seq[pos]+prot_seq[pos+1:pos+6]
                            alt_seq=prot_seq[0:pos]+prot_seq[pos]+alt+prot_seq[pos+1:pos+6]
                        elif (pos+6)>(len(prot_seq)-1):
                            ref_seq=prot_seq[pos-5:pos]+prot_seq[pos]+prot_seq[pos+1:len(prot_seq)+1]
                            alt_seq=prot_seq[pos-5:pos]+prot_seq[pos]+alt+prot_seq[pos+1:len(prot_seq)+1]
                                                   
                        h_diff=diffwwH(ref_seq,alt_seq)                        
                        fl_diff=diff_flex(ref_seq,alt_seq)
                        
                        
                    #else consider it an error and append 'N/A's
                    else:                     
                        h_diff='NaN'
                        fl_diff='NaN'
                        
        ############# if codon change and insertion
        elif ref not in ['N/A','X','*','-'] and alt not in ['N/A','X','*','-'] and\
        'X' not in ref and 'X' not in alt and len(alt)>len(ref):  
            
            pos=int(pos)-1 #substracting 1 as string starts with zero index in python             
            prot_seq=''
            ref_seq=''
            alt_seq=''
               
            #calculate the scale differences
            for key in prot_dict:
                if uni==key:
                    prot_seq=prot_dict[key]
                    #if the position is inside the protein length
                    if pos<len(prot_seq):  
                        if (pos-5)>=0 and (pos+6)<=(len(prot_seq)-1): 
                            ref_seq=prot_seq[pos-5:pos]+prot_seq[pos]+prot_seq[pos+1:pos+6]
                            alt_seq=prot_seq[pos-5:pos]+alt+prot_seq[pos+1:pos+6]
                        elif (pos-5)<0:
                            ref_seq=prot_seq[0:pos]+prot_seq[pos]+prot_seq[pos+1:pos+6]
                            alt_seq=prot_seq[0:pos]+alt+prot_seq[pos+1:pos+6]
                        elif (pos+6)>(len(prot_seq)-1):
                            ref_seq=prot_seq[pos-5:pos]+prot_seq[pos]+prot_seq[pos+1:len(prot_seq)+1]
                            alt_seq=prot_seq[pos-5:pos]+alt+prot_seq[pos+1:len(prot_seq)+1]
                            
                        
                        h_diff=diffwwH(ref_seq,alt_seq)                        
                        fl_diff=diff_flex(ref_seq,alt_seq)
                        
                    #else consider it an error and append 'N/A's
                    else:                        
                        h_diff='NaN'
                        fl_diff='NaN'
                        
       
        ############## if deletion only
        elif alt=='-' and ref not in ['N/A','X','*','-'] and 'X' not in ref:  
            
            pos=int(pos)-1 #substracting 1 as string starts with zero index in python 
            prot_seq=''
            ref_seq=''
            alt_seq=''
            
            #calculate the scale differences
            for key in prot_dict:
                if uni==key:
                    
                    prot_seq=prot_dict[key]
                    
                    #if the position is inside the protein length
                    if pos<len(prot_seq): 
                        
                        if (pos-5)>=0 and (pos+len(ref)+5)<=(len(prot_seq)-1): 
                            alt_seq=prot_seq[pos-5:pos]+prot_seq[pos+len(ref):pos+len(ref)+5]
                            ref_seq=prot_seq[pos-5:pos]+prot_seq[pos:pos+len(ref)+5]
                        elif (pos-5)<0:
                            
                            alt_seq=prot_seq[0:pos]+prot_seq[pos+len(ref):pos+len(ref)+5]
                            ref_seq=prot_seq[0:pos]+prot_seq[pos:pos+len(ref)+5]
                        elif (pos+len(ref)+5)>len(prot_seq)-1:
                            alt_seq=prot_seq[pos-5:pos]+prot_seq[pos+len(ref):len(prot_seq)+1]
                            ref_seq=prot_seq[pos-5:pos]+prot_seq[pos:len(prot_seq)+1]
                        
                        
                        h_diff=diffwwH(ref_seq,alt_seq)                        
                        fl_diff=diff_flex(ref_seq,alt_seq)                      
                    
                    #else consider it an error and append 'N/A's
                    else:
                        vol_diff='NaN'
                        h_diff='NaN'
                        fl_diff='NaN'
                        pol_diff='NaN'
        
        ############## if codon change and deletion
        elif alt not in ['N/A','X','*','-'] and ref not in ['N/A','X','*','-'] and\
        'X' not in ref and 'X' not in alt and len(ref)>len(alt):             
            
            pos=int(pos)-1 #substracting 1 as string starts with zero index in python 
            prot_seq=''
            ref_seq=''
            alt_seq=''
            
            #calculate the scale differences
            for key in prot_dict:
                if uni==key:
                    prot_seq=prot_dict[key]
                    #if the position is inside the protein length
                    if pos<len(prot_seq): 
                        if (pos-5)>=0 and (pos+len(ref)+5)<=(len(prot_seq)-1): 
                            alt_seq=prot_seq[pos-5:pos]+alt+prot_seq[pos+len(ref)-len(alt):pos+len(ref)+5]
                            ref_seq=prot_seq[pos-5:pos]+prot_seq[pos:pos+len(ref)+5]
                        elif (pos-5)<0:
                            alt_seq=prot_seq[0:pos]+alt+prot_seq[pos+len(ref)-len(alt):pos+len(ref)+5]
                            ref_seq=prot_seq[0:pos]+prot_seq[pos:pos+len(ref)+5]
                        elif (pos+len(ref)+5)>len(prot_seq)-1:
                            alt_seq=prot_seq[pos-5:pos]+alt+prot_seq[pos+len(ref)-len(alt):len(prot_seq)+1]
                            ref_seq=prot_seq[pos-5:pos]+prot_seq[pos:len(prot_seq)+1]                        
                       
                        h_diff=diffwwH(ref_seq,alt_seq)                        
                        fl_diff=diff_flex(ref_seq,alt_seq)
                        
                    
                    #else consider it an error and append 'N/A's
                    else:                        
                        h_diff='NaN'
                        fl_diff='NaN'
                        
             
        ############## if intronic or splicing variant
        elif alt=='-' and ref=='-':
            h_diff='0.0'
            fl_diff='0.0'

        else:
            h_diff = 'NaN'
            fl_diff = 'NaN'
           
        
        #write aa differences in the output file
        output_file.write(line.rstrip('\n').rstrip('\r')+'\t'+str(fl_diff)+'\t'+str(h_diff)+'\n')

 
    variants_file.close()
    output_file.close()
            
        

