import sys

barcode_oligo_txt_file = "/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/associations/MH_association.count_ready.tsv"
fasta_file = "/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/associations/oligos_barcodes.fa"


with open(barcode_oligo_txt_file) as f:
    
    with open(fasta_file, 'w') as outfile:
        count = 0
        for line in f:
            count += 1
            if count == 1:  # skip the header line
                continue
            info = [x.strip() for x in line.split(',')]
            oligo = info[0]
            bc_seq = info[1]
            outfile.write('>'+oligo+'_'+bc_seq)
            outfile.write('\n')
            outfile.write(bc_seq)
            outfile.write('\n')
            
