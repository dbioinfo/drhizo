#mothur uses its own console and commands, so use these commands after the mothur command

#locate fastq
make.file(inputdir=raw_data/fastq/, type=fastq, prefix=raw_data/fastq/stability)
#match into contigs
make.contigs(file=raw_data/fastq/stability.files)
#check lengths
summary.seqs(fasta=current, count=current)
#filter bad seqs
screen.seqs(fasta=current, count=current, maxambig=0, maxlength=350, maxhomop=8)
#check em again
summary.seqs(fasta=current, count=current)
#reduce dataset by unique seq
unique.seqs(fasta=raw_data/fastq/stability.trim.contigs.good.fasta, count=raw_data/fastq/stability.trim.contigs.good.count_table)
#check em again
summary.seqs(count=raw_data/fastq/stability.trim.contigs.good.count_table)
#from here, need to download aligned silva ref
pcr.seqs(fasta=raw_data/silva/

