#!/usr/bin/env python3

import os
import pytest
import shutil
import subprocess

import pyfastaq


COVID_REF_FASTA = "MN908947.3.fa"
READ_IT_AND_KEEP = os.path.join(os.pardir, "src", "readItAndKeep")
ART_ILLUMINA = "art_illumina"
assert shutil.which(READ_IT_AND_KEEP) is not None
assert shutil.which(ART_ILLUMINA) is not None


def riak_debug_check_badreads_failed_reads_ok(debug_file, identity_cutoff=86):
    with open(debug_file) as f:
        for line in f:
            if not line.startswith("REJECTED_READ\t"):
                continue

            description = line.split("\t")[2]
            if description.startswith("junk_seq") or description.startswith("random_seq"):
                continue

            # description example:
            # MN908947.3,-strand,26462-27502 length=971 error-free_length=1060 read_identity=70.47%
            read_identity = description.split()[-1]
            assert read_identity.startswith("read_identity=")
            percent_identity = float(read_identity.split("=")[-1].strip("%"))
            if percent_identity > identity_cutoff:
                return False

    return True


def run_read_it_and_keep(reads_files, outprefix, tech):
    assert 1 <= len(reads_files) <= 2
    reads = " ".join([f"--reads{i+1} " + f for i, f in enumerate(reads_files)])
    command = f"{READ_IT_AND_KEEP} --debug --tech {tech} --ref_fasta {COVID_REF_FASTA} {reads} -o {outprefix} > {outprefix}.out"
    subprocess.check_output(command, shell=True)
    results = {}
    with open(f"{outprefix}.out") as f:
        for line in f:
            key, val = line.rstrip().split("\t")
            assert key not in results
            results[key] = int(val)
    return results


def shred_ref_genome(outfile, k):
    seqs = {}
    pyfastaq.tasks.file_to_dict(COVID_REF_FASTA, seqs)
    assert len(seqs) == 1
    ref = list(seqs.values())[0]
    k = 75
    with open(outfile, "w") as f:
        for i in range(len(ref)):
            end = i + k - 1
            if end > len(ref) - 1:
                break
            print(f">{i+1}.{end+1}", file=f)
            print(ref[i: end+1], file=f)


def simulate_nanopore(outfile, coverage=10, read_length_stdev="15000,13000"):
    """Uses 'badread' to simulate nanopore reads. The default read_length_stdev
    for this function is the same as the default for badread."""
    coverage = 10
    command = " ".join([
        "badread simulate",
        "--seed 42",
        "--reference", COVID_REF_FASTA,
        "--quantity", f"{coverage}x",
        "--length", read_length_stdev,
        ">", outfile
    ])
    subprocess.check_output(command, shell=True)


def simulate_illumina(outprefix):
    machine = "HS25"  # This is HiSeq 2500 (125bp, 150bp)
    read_length = 150
    read_depth = 10
    mean_fragment_length = 250
    fragment_length_sd = 10
    command = " ".join(
        [
            ART_ILLUMINA,
            "--in",
            COVID_REF_FASTA,
            "--out",
            outprefix,
            "--noALN",  # do not output alignment file
            "--seqSys",
            machine,
            "--len",
            str(read_length),
            "--fcov",
            str(read_depth),
            "--mflen",
            str(mean_fragment_length),
            "--sdev",
            str(fragment_length_sd),
            "--rndSeed 42",
        ]
    )
    subprocess.check_output(command, shell=True)


def test_shredded_ref_genome():
    outprefix = "out.shredded_ref_genome"
    subprocess.check_output(f"rm -rf {outprefix}*", shell=True)
    shred_length = 75
    shredded_ref = f"{outprefix}.k-{shred_length}.fa"
    shred_ref_genome(shredded_ref, shred_length)
    riak_out = f"{outprefix}.riak"
    riak_counts = run_read_it_and_keep([shredded_ref], riak_out, "illumina")
    number_of_reads = 29829
    assert riak_counts["Input reads file 1"] == number_of_reads
    assert riak_counts["Input reads file 2"] == 0
    assert riak_counts["Kept reads 1"] == number_of_reads
    assert riak_counts["Kept reads 2"] == 0
    reads_out = f"{riak_out}.reads.fasta.gz"
    assert os.path.exists(reads_out)
    assert pyfastaq.tasks.count_sequences(reads_out) == number_of_reads
    # check it didn't make output reads files that would only be made when
    # there were two input files
    reads_out_1 = f"{riak_out}.reads_1.fasta.gz"
    assert not os.path.exists(reads_out_1)
    reads_out_2 = f"{riak_out}.reads_2.fasta.gz"
    assert not os.path.exists(reads_out_2)
    subprocess.check_output(f"rm -rf {outprefix}*", shell=True)


def test_illumina_hiseq_sim():
    outprefix = "out.illumina_hiseq_sim"
    subprocess.check_output(f"rm -rf {outprefix}*", shell=True)
    illumina_read_prefix = f"{outprefix}.reads"
    simulate_illumina(illumina_read_prefix)
    illumina_1 = f"{illumina_read_prefix}1.fq"
    illumina_2 = f"{illumina_read_prefix}2.fq"
    riak_out = f"{outprefix}.riak"
    riak_counts = run_read_it_and_keep([illumina_1, illumina_2], riak_out, "illumina")
    number_of_reads = 995
    assert riak_counts["Input reads file 1"] == number_of_reads
    assert riak_counts["Input reads file 2"] == number_of_reads
    assert riak_counts["Kept reads 1"] == number_of_reads
    assert riak_counts["Kept reads 2"] == number_of_reads
    reads_out_1 = f"{riak_out}.reads_1.fastq.gz"
    assert os.path.exists(reads_out_1)
    assert pyfastaq.tasks.count_sequences(reads_out_1) == number_of_reads
    reads_out_2 = f"{riak_out}.reads_2.fastq.gz"
    assert os.path.exists(reads_out_2)
    assert pyfastaq.tasks.count_sequences(reads_out_2) == number_of_reads
    # check it didn't make output read file that would be made if input
    # was only one file
    reads_out = f"{riak_out}.reads.fasta.gz"
    assert not os.path.exists(reads_out)
    subprocess.check_output(f"rm -rf {outprefix}*", shell=True)


def test_nanopore_sim_default():
    outprefix = "out.nanopore_sim_default"
    subprocess.check_output(f"rm -rf {outprefix}*", shell=True)
    nano_reads = f"{outprefix}.reads.fq"
    simulate_nanopore(nano_reads)
    number_of_input_reads = 35
    assert pyfastaq.tasks.count_sequences(nano_reads) == number_of_input_reads
    riak_out = f"{outprefix}.riak"
    riak_counts = run_read_it_and_keep([nano_reads], riak_out, "ont")
    reads_out = f"{riak_out}.reads.fastq.gz"
    number_of_output_reads = 32
    assert os.path.exists(reads_out)
    assert riak_debug_check_badreads_failed_reads_ok(f"{riak_out}.reads.debug")
    assert pyfastaq.tasks.count_sequences(reads_out) == number_of_output_reads
    assert riak_counts["Input reads file 1"] == number_of_input_reads
    assert riak_counts["Input reads file 2"] == 0
    assert riak_counts["Kept reads 1"] == number_of_output_reads
    assert riak_counts["Kept reads 2"] == 0
    subprocess.check_output(f"rm -rf {outprefix}*", shell=True)


def test_nanopore_sim_1k_reads():
    outprefix = "out.nanopore_sim_1k_reads"
    subprocess.check_output(f"rm -rf {outprefix}*", shell=True)
    nano_reads = f"{outprefix}.reads.fq"
    simulate_nanopore(nano_reads, read_length_stdev="1000,50")
    number_of_input_reads = 301
    assert pyfastaq.tasks.count_sequences(nano_reads) == number_of_input_reads
    riak_out = f"{outprefix}.riak"
    riak_counts = run_read_it_and_keep([nano_reads], riak_out, "ont")
    reads_out = f"{riak_out}.reads.fastq.gz"
    number_of_output_reads = 282
    assert os.path.exists(reads_out)
    assert riak_debug_check_badreads_failed_reads_ok(f"{riak_out}.reads.debug")
    assert pyfastaq.tasks.count_sequences(reads_out) == number_of_output_reads
    assert riak_counts["Input reads file 1"] == number_of_input_reads
    assert riak_counts["Input reads file 2"] == 0
    assert riak_counts["Kept reads 1"] == number_of_output_reads
    assert riak_counts["Kept reads 2"] == 0
    subprocess.check_output(f"rm -rf {outprefix}*", shell=True)

