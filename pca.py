OA#11.25.19 extract n-mer from mutated protein 

import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.Seq import Seq
from Bio import Alphabet
from Bio.Alphabet import IUPAC, ProteinAlphabet

def get_prot(swissprot):

    uni_prot = swissprot

    import httplib2 as http
    import json

    try:
        from urlparse import urlparse
    except ImportError:
        from urllib.parse import urlparse

    headers = {
            'Accept': 'application/rdf+xml',
    }

    uri = 'https://www.uniprot.org'
    path = '/uniprot/%s.fasta'%(uni_prot)
    target = urlparse(uri+path)
    method = 'GET'
    body = ''

    h = http.Http()

    response, content = h.request(
            target.geturl(),
            method,
            body,
            headers)

    record = "".join(content.decode("utf-8").split('\n')[1:])
    return(record)

def get_pca(test_id,pca_size):

    flank = int(pca_size / 2)

    mut_test = mut_parsed[mut_parsed.patient_id == test_id]
    mut_test = mut_test.fillna(0)
    mut_filt = mut_test[~np.array(mut_test.SWISSPROT == 0)]
    mut_filt = mut_filt[~np.array(mut_filt.HGVSp == 0)]
    ctr = 0
    failed = 0 

    psuedo_cancer_antigen = []
    

    for index,row in mut_filt.iterrows():
        mutation_index = int(''.join(list(filter(str.isdigit, row.HGVSp))))
        try:
            prot_seq = get_prot(row.SWISSPROT)
            mut_nucleotide = row.Amino_acids.split("/")[1]
            new_seq = prot_seq[:mutation_index-1] + mut_nucleotide + prot_seq[mutation_index:]      
            native_antigen = prot_seq[mutation_index-1 - flank : mutation_index-1 + flank]
            tumor_antigen = new_seq[mutation_index-1 -flank : mutation_index-1 + flank].upper()

            if Alphabet._verify_alphabet(Seq(tumor_antigen,Alphabet.IUPAC.protein)):
                psuedo_cancer_antigen += [tumor_antigen]
            else:
                print(test_id,row.Hugo_Symbol,tumor_antigen,'Invalid pca')

        except:
            print(test_id,row.Hugo_Symbol,"Failed to get antigen")
            
        
    len_bool = [len(a) >= 8 and len(a) <= 15 for a in psuedo_cancer_antigen]
    filt_pca = list(np.array(psuedo_cancer_antigen)[len_bool])

    return filt_pca

'''
Load mini_challenge_data

'''

DATA = "mini_challenge_data/"

#clin_train = pd.read_csv(DATA + "clinical_test.txt",sep = "\t")
mut_test = pd.read_csv(DATA + "mut_train.txt", sep="\t")

cols = ['patient_id','Amino_acids','Variant_Type','Variant_Classification','IMPACT','Hugo_Symbol','HGVSc','HGVSp','SWISSPROT']
mut_parsed = mut_test[cols]

#mut_parsed = mut_parsed[~(mut_parsed.Variant_Classification.isin(
#    ["RNA","Silent", "5'Flank", "Intron", "IGR", "3'UTR", "5'UTR", "3'Flank", "Splice_Site"]))]

#mut_parsed = mut_parsed[mut_parsed.Variant_Type == 'SNP']

mut_parsed = mut_parsed[mut_parsed.Variant_Classification == 'Missense_Mutation']

from mhcflurry import Class1AffinityPredictor
import operator
predictor = Class1AffinityPredictor.load()

import sys
test_id = str(sys.argv[1])                                                                                                                       
filt_pca = get_pca(test_id,10)

print(test_id,len(filt_pca),filt_pca,"\n")
prediction = predictor.predict(allele="HLA-A0101", peptides= filt_pca)
zipped = zip(filt_pca,prediction)

prediction_sorted = sorted(zipped, key = operator.itemgetter(1)) 
print(prediction_sorted,"\n")
print(prediction_sorted[0])
    
