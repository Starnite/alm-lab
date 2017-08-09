#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 09:29:32 2017

@author: vivzhong
"""

import click
import pandas as pd
import numpy as np
@click.command()
@click.option('--fasta', help="Fasta file of quality filtered/trimmed sequences identified as 16S sequences.")
@click.option('--full_fastq', help="The corresponding, original fastq file with metagenomic data.")
@click.option('--filt_fastq', help="The new fastq file to be written.")
@click.option('--trunclen', help="The length of the sequences after quality filtering/trimming.")

def quality_matcher(fasta, full_fastq, filt_fastq, trunclen):
    """
    Get the quality scores for sequences in a fasta file and write to a new fastq file.
    """
    fasta_seq = pd.read_csv(fasta, header=None, engine="c", skiprows=lambda x:x%2==0, names=["Seq"])
    fasta_id = pd.read_csv(fasta, header=None, engine="c", skiprows=lambda x:x%2!=0, names=["ID"])
    fastq_id = pd.read_csv(full_fastq, header=None, engine="c", skiprows=lambda x:x%4!=0, names=["ID"])
    fastq_qual = pd.read_csv(full_fastq, header=None, engine="c", skiprows=lambda x:(x+1)%4!=0, names=["Qual"])
#    fasta_seq = pd.read_csv(fasta, header=None, engine="c", names=["Seq"]).iloc[1::2,:].reset_index(drop=True)
#    fasta_id = pd.read_csv(fasta, header=None, engine="c", names=["ID"]).iloc[::2].reset_index(drop=True)
#    fastq_id = pd.read_csv(full_fastq, header=None, engine="c", names=["ID"]).iloc[::4].reset_index(drop=True)
#    fastq_qual = pd.read_csv(full_fastq, header=None, engine="c", names=["Qual"]).iloc[3::4].reset_index(drop=True)
    fasta_id["ID"] = fasta_id["ID"].map(lambda x: x.lstrip('>'))
    fastq_id["ID"] = fastq_id["ID"].map(lambda x: x.lstrip('@'))
    fasta = pd.concat([fasta_id,fasta_seq], axis=1).set_index("ID")
    fastq = pd.concat([fastq_id, fastq_qual], axis=1).set_index("ID")
    fastq["Qual"] = fastq["Qual"].map(lambda x:str(x)[:int(trunclen)])
    new_fastq = fasta
    new_fastq["IDQ"] = new_fastq.index.map(lambda x: "+{}".format(x))
    new_fastq["Qual"] = fastq["Qual"]
    new_fastq.index = new_fastq.index.map(lambda x: "@{}".format(x))
    new_fastq = new_fastq.reset_index()
    np.savetxt(filt_fastq, new_fastq.values.flatten(), fmt='%s', newline="\n")


if __name__ == '__main__':
    quality_matcher()

