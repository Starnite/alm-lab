# -*- coding: utf-8 -*-
"""
author: Vivian Zhong
date: August 7, 2017
"""

import click
@click.command()
@click.option('--fasta', help="Fasta file of quality filtered/trimmed sequences identified as 16S sequences.")
@click.option('--full_fastq', help="The corresponding, original fastq file with metagenomic data.")
@click.option('--filt_fastq', help="The new fastq file to be written.")
@click.option('--trunclen', help="The length of the sequences after quality filtering/trimming.")

def quality_matcher(fasta, full_fastq, filt_fastq, trunclen):
    """
    Get the quality scores for sequences in a fasta file and write to a new fastq file.
    """
    with open(fasta, "r") as fasta, open(full_fastq, "r") as fastq, open(filt_fastq, "w") as new_fastq:
        #make lists of the fasta and fastq files, where every successive value is a successive line
        #purpose of -1: to avoid the "\n" newline character at the end of the lines
        fastq_list = [line[:-1] for line in fastq]
        fasta_list = [line[:-1] for line in fasta]
        #iterate through the sequence ids in the fasta file
        for fasta_index, fasta_id in enumerate(fasta_list):
            if fasta_id[0] == ">":
                #get the list index of the matching sequence id in the metagenomic fastq file
                fastq_index = fastq_list.index("@{}".format(fasta_id[1:]))
                #print and write a new fastq entry with the quality scores string truncated to the same length as the sequence from the fasta file
                print(str("@{}".format(fasta_id[1:])) + "\n" + str(fasta_list[fasta_index+1]) + "\n" + str("+{}".format(fasta_id[1:])) + "\n" + str(fastq_list[fastq_index+3][:int(trunclen)]))
                new_fastq.write(str("@{}".format(fasta_id[1:])) + "\n" + str(fasta_list[fasta_index+1]) + "\n" + str("+{}".format(fasta_id[1:])) + "\n" + str(fastq_list[fastq_index+3][:int(trunclen)]))

if __name__ == '__main__':
    quality_matcher()
