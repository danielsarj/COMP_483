# loading libraries
import os
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline

# checking if results folder exists. creates one if it doesnt. 
current_dir = os.getcwd()
if os.path.isdir(current_dir+'/results') == False:
    os.makedirs(current_dir+'/results')
results_dir = current_dir+'/results'

# creates log file
log_output = open(results_dir+'/miniproject.log', 'w')

# retrieving SRR8185310 K-12 SRA Illimuna reads
print('INFO: Attempting to fetch SRR8185310 SRA Illumina reads.')
os.system('prefetch SRR8185310')
if os.path.isdir(current_dir+'/SRR8185310'):
    print('INFO: SRA Illumina reads successfully downloaded.')
    os.system('fasterq-dump SRR8185310')
    if os.path.isfile(current_dir+'/SRR8185310.fastq'):
        print('INFO: SRA format file successfully translated into fastq format.')
    else:
        print('INFO: Could not translate SRA format to fastq format. Check if SRA toolkit has been correctly installed and try again.')
        exit()
else:
    print('INFO: Could not download SRA Illumina reads. Check if SRA toolkit has been correctly installed or if the given accession number is correct, and try again.')
    exit()
os.system('mv '+current_dir+'/SRR8185310 '+results_dir)
os.system('mv '+current_dir+'/SRR8185310.fastq '+results_dir)

# assembling genome using SPAdes
print('INFO: Attempting to assemble genome using SPAdes.')
spades_command = 'python3 '+current_dir+'/SPAdes-3.15.4-Linux/bin/spades.py -t 2 -k 33,55,77,99,127 --only-assembler -s '+results_dir+'/SRR8185310.fastq -o '+results_dir+'/SRR8185310_assembly/'
log_output.write(spades_command+'\n')
os.system(spades_command)
if os.path.isfile(results_dir+'/SRR8185310_assembly/scaffolds.fasta'):
    print('INFO: Genome successfully assembled.')
else:
    print('INFO: Could not assemble genome. Check input files and try again.')
    exit()

# assessing genome assembly info
records = list(SeqIO.parse(results_dir+'/SRR8185310_assembly/scaffolds.fasta', 'fasta'))
final_contigs = 0
assembly_size = 0
with open(results_dir+'/SRR8185310_assembly_1000bp.fasta', 'w') as outfile:
    for item in records:
        if len(str(item.seq)) > 1000:
            outfile.write('>'+str(item.id)+'\n'+str(item.seq)+'\n')
            final_contigs += 1
            assembly_size += len(str(item.seq))

log_output.write('There are '+str(final_contigs)+' contigs > 1000 in the assembly.\n')
log_output.write('There are '+str(assembly_size)+' bp in the assembly.\n')
print('There are '+str(final_contigs)+' contigs > 1000 in the assembly.')
print('There are '+str(assembly_size)+' bp in the assembly.')

# annotating genome using GeneMarkS-2
print('INFO: Attempting to annotate genome using GeneMarkS-2.')
gms2_command = 'perl '+current_dir+'/gms2_linux_64/gms2.pl --seq '+results_dir+'/SRR8185310_assembly_1000bp.fasta --genome-type bacteria --output '+results_dir+'/SRR8185310_anno.gtf --format gtf --faa '+results_dir+'/SRR8185310_anno_AA.fasta'
os.system(gms2_command)
if os.path.isfile(results_dir+'/SRR8185310_anno.gtf'):
    print('INFO: Genome successfully annotated.')
else:
    print('INFO: Could not assemble genome. Check input files and try again.')
    exit()

# predicting protein function using BLAST
print('INFO: Attempting to predict function of annotated proteins using BLAST.')
if os.path.isfile(current_dir+'/Ecoli_db.phr') == False and os.path.isfile(current_dir+'/Ecoli_db.pin') == False and os.path.isfile(current_dir+'/Ecoli_db.psq') == False:
    print('INFO: BLAST database not found. Attempting to make one.')
    os.system('makeblastdb -in Ecoli.fasta -out Ecoli_db -title Ecoli_db -dbtype prot')
    if os.path.isfile(current_dir+'/Ecoli_db.phr') == False and os.path.isfile(current_dir+'/Ecoli_db.pin') == False and os.path.isfile(current_dir+'/Ecoli_db.psq') == False:
        print('INFO: BLAST database could not be created. Check if the reference FASTA file exists and is in the working folder, or if makeblastdb is installed.\n')
        exit()
    else: 
        print('INFO: BLAST database successfully created. Will run BLAST next.')
blast_command = NcbiblastpCommandline(query = results_dir+'/SRR8185310_anno_AA.fasta', db = 'Ecoli_db', num_threads = 2, max_target_seqs = 1, outfmt = '10 qseqid sseqid pident qcovs', out = results_dir+'/SRR8185310_blastp_results.csv')
stdout, stderr = blast_command()
if os.path.isfile(results_dir+'/SRR8185310_blastp_results.csv'):
    print('INFO: Successfully ran BLAST.')
else:
    print('INFO: Could not run BLAST. Check input files and try again.')
    exit()
os.system('mv '+current_dir+'/Ecoli_db.phr '+results_dir)
os.system('mv '+current_dir+'/Ecoli_db.pin '+results_dir)
os.system('mv '+current_dir+'/Ecoli_db.psq '+results_dir)

with open(results_dir+'/SRR8185310_blastp_results.csv', 'r') as blast_file:
    last_line = blast_file.readlines()[-1]
    cds_n = int(last_line.split(',')[0])
    if cds_n > 4140:
        difference = cds_n - 4140
        log_output.write('GeneMarkS-2 found '+str(difference)+' more CDS than the RefSeq.')
        print('GeneMarkS-2 found '+str(difference)+' more CDS than the RefSeq.')
    elif cds_n < 4140:
        difference = 4140 - cds_n
        log_output.write('GeneMarkS-2 found '+str(difference)+' less CDS than the RefSeq.')
        print('GeneMarkS-2 found '+str(difference)+' less CDS than the RefSeq.')
    else: 
        log_output.write('GeneMarkS-2 found the same number of CDS than the RefSeq.')
        print('GeneMarkS-2 found the same number of less CDS than the RefSeq.')
    
# retrieving SRR1411276 K-12 derivative BW38028 SRA transcriptomic project reads 
print('INFO: Attempting to fetch SRR1411276 SRA Illumina reads.')
os.system('prefetch SRR1411276')
if os.path.isdir(current_dir+'/SRR1411276'):
    print('INFO: SRA Illumina reads successfully downloaded.')
    os.system('fasterq-dump SRR1411276')
    if os.path.isfile(current_dir+'/SRR1411276.fastq'):
        print('INFO: SRA format file successfully translated into fastq format.')
    else:
        print('INFO: Could not translate SRA format to fastq format. Check if SRA toolkit has been correctly installed and try again.')
        exit()
else:
    print('INFO: Could not download SRA Illumina reads. Check if SRA toolkit has been correctly installed or if the given accession number is correct, and try again.')
    exit()
os.system('mv '+current_dir+'/SRR1411276 '+results_dir)
os.system('mv '+current_dir+'/SRR1411276.fastq '+results_dir)

print('INFO: Attempting to index E. coli K-12 assembled genome using Bowtie.')
os.system('bowtie-build '+results_dir+'/SRR8185310_assembly_1000bp.fasta NC_000913')
if os.path.isfile(current_dir+'/NC_000913.1.ebwt') and os.path.isfile(current_dir+'/NC_000913.2.ebwt') and os.path.isfile(current_dir+'/NC_000913.3.ebwt') and os.path.isfile(current_dir+'/NC_000913.4.ebwt') and os.path.isfile(current_dir+'/NC_000913.rev.1.ebwt') and os.path.isfile(current_dir+'/NC_000913.rev.2.ebwt'):
    print('INFO: E. coli K-12 assembled genome successfully indexed.')
else:
    print('INFO: Could not index E. coli K-12 assembled genome. Check if input file exists or Bowtie is correctly installed.')
    exit()

# mapping transcriptomic reads to reference genome & quantifying expression
print('INFO: Attempting to map SRR1411276 K-12 derivative BW38028 SRA transcriptomic project reads to the E. coli K-12 assembled genome.')
os.system('tophat2 --no-novel-juncs -o '+results_dir+'/SRR1411276_mapped_NC_000913 NC_000913 '+results_dir+'/SRR1411276.fastq')
if os.path.isfile(results_dir+'/SRR1411276_mapped_NC_000913/accepted_hits.bam'):
    print('INFO: Mapping was successful.\n')
else:
    print('INFO: Could not map reads. Check if input file exists or Tophat2 is correctly installed.')
    exit()
os.system('mv '+current_dir+'/NC_000913.fa '+results_dir)
os.system('mv '+current_dir+'/NC_000913.1.ebwt '+results_dir)
os.system('mv '+current_dir+'/NC_000913.2.ebwt '+results_dir)
os.system('mv '+current_dir+'/NC_000913.3.ebwt '+results_dir)
os.system('mv '+current_dir+'/NC_000913.4.ebwt '+results_dir)
os.system('mv '+current_dir+'/NC_000913.rev.1.ebwt '+results_dir)
os.system('mv '+current_dir+'/NC_000913.rev.2.ebwt '+results_dir)

print('INFO: Attempting to quantify expression of mapped reads using Cufflinks.')
os.system('cufflinks -o '+results_dir+'/SRR1411276_cufflinks -p 2 -G '+results_dir+'/SRR8185310_anno.gtf '+results_dir+'/SRR1411276_mapped_NC_000913/accepted_hits.bam')
if os.path.isfile(results_dir+'/SRR1411276_cufflinks/transcripts.gtf'):
    print('INFO: Expression successfully quantified.')
else:
    print('INFO: Could not quantify expression of mapped reads. Check if input file exists or Cufflink is correctly installed.')
    exit()

print('INFO: Summarizing Cufflinks output.')
with open(results_dir+'/SRR1411276_cufflinks/transcripts.gtf','r') as cuff_infile:
    with open(results_dir+'/transcriptome_data.fpkm', 'w') as cuff_outfile:
        for line in cuff_infile:
            seqname = ''
            start = ''
            end = '' 
            strand = ''
            FPKM = ''
            line = line.split('\t')
            if line[2] == 'transcript':
                seqname = line[0]
                start = line[3]
                end = line[4]
                strand = line[6]
                att = line[8].split(' ')
                FPKM = att[7].replace('"','').replace(';','')
                cuff_outfile.write(seqname+','+start+','+end+','+strand+','+FPKM+'\n')

if os.path.isfile(results_dir+'/transcriptome_data.fpkm'):
    print('INFO: Cufflinks output successfully summarized.')
    print('Wrapper successfully finished running. Have a nice day. :)')
else:
    print('INFO: Could not summarize Cufflinks output. Check if input file exists.')
    exit()

# closes log file
log_output.close()