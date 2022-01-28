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
seqs = {}
pyfastaq.tasks.file_to_dict(COVID_REF_FASTA, seqs)
assert len(seqs) == 1
COVID_REF_SEQ = list(seqs.values())[0]


def file_to_lines(filename):
    f = pyfastaq.utils.open_file_read(filename)
    lines = [x.rstrip() for x in f]
    f.close()
    return lines


def riak_debug_check_badreads_failed_reads_ok(
    debug_file, identity_cutoff=87, min_length=50
):
    with open(debug_file) as f:
        for line in f:
            if not line.startswith("REJECTED_READ\t"):
                continue

            _, _, description, seq, _ = line.rstrip().split("\t")
            # description = line.split("\t")[2]
            if description.startswith("junk_seq") or description.startswith(
                "random_seq"
            ):
                continue

            # description example:
            # MN908947.3,-strand,26462-27502 length=971 error-free_length=1060 read_identity=70.47%
            read_identity = description.split()[-1]
            assert read_identity.startswith("read_identity=")
            percent_identity = float(read_identity.split("=")[-1].strip("%"))
            if percent_identity > identity_cutoff and len(seq) >= min_length:
                print("REJECT FAIL:", line)
                return False

    return True


def run_read_it_and_keep(reads_files, outprefix, tech, enumerate_names=False):
    enumerate_opt = "--enumerate_names" if enumerate_names else ""
    assert 1 <= len(reads_files) <= 2
    reads = " ".join([f"--reads{i+1} " + f for i, f in enumerate(reads_files)])
    command = f"{READ_IT_AND_KEEP} --debug {enumerate_opt} --tech {tech} --ref_fasta {COVID_REF_FASTA} {reads} -o {outprefix} > {outprefix}.out"
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
    k = 75
    with open(outfile, "w") as f:
        for i in range(len(COVID_REF_SEQ)):
            end = i + k - 1
            if end > len(COVID_REF_SEQ) - 1:
                break
            print(f">{i+1}.{end+1}", file=f)
            print(COVID_REF_SEQ[i : end + 1], file=f)


def simulate_nanopore(
    outfile, coverage=10, read_length_stdev="15000,13000", ref=COVID_REF_FASTA
):
    """Uses 'badread' to simulate nanopore reads. The default read_length_stdev
    for this function is the same as the default for badread."""
    coverage = 10
    command = " ".join(
        [
            "badread simulate",
            "--seed 42",
            "--reference",
            ref,
            "--quantity",
            f"{coverage}x",
            "--length",
            read_length_stdev,
            ">",
            outfile,
        ]
    )
    subprocess.check_output(command, shell=True)


def simulate_illumina(outprefix, ref=COVID_REF_FASTA):
    machine = "HS25"  # This is HiSeq 2500 (125bp, 150bp)
    read_length = 150
    read_depth = 10
    mean_fragment_length = 250
    fragment_length_sd = 10
    command = " ".join(
        [
            ART_ILLUMINA,
            "--in",
            ref,
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


def make_ref_genome_with_deletion(outfile, del_start, del_end, keep_start, keep_end):
    """Writes a FASTA file of the covid reference genome, but with positions
    del_start to del_end inclusive deleted. Only writes out the part of the
    genome from keep_start to keep_end - we only really care about testing if
    reads that contain the deletion are kept by readItAndKeep.
    Coordinates are 1-based."""
    seqs = {}
    pyfastaq.tasks.file_to_dict(COVID_REF_FASTA, seqs)
    assert len(seqs) == 1
    ref = list(seqs.values())[0]
    ref.seq = ref[keep_start - 1 : del_start - 1] + ref[del_end:keep_end]
    with open(outfile, "w") as f:
        print(ref, file=f)


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


# The next tests is a deletion that is seen in real samples. This is a
# good summary of known deletions, and has the coordinates to use:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7577707/
# We will simulate reads from the covid reference with the deletion added,
# and check that readItAndKeep keeps the reads.

# "Delta 382" deletion is deletion from 27848 to 28229. It's the longest
# deletion in that paper.
def test_delta_382_deletion():
    outprefix = "out.delta_382_deletion"
    subprocess.check_output(f"rm -rf {outprefix}*", shell=True)
    sim_reads_fasta = f"{outprefix}.for_reads.fa"
    make_ref_genome_with_deletion(sim_reads_fasta, 27848, 28229, 26848, 29229)

    nano_reads = f"{outprefix}.reads.fq"
    simulate_nanopore(nano_reads, read_length_stdev="1000,50", ref=sim_reads_fasta)
    number_of_input_reads = 32
    assert pyfastaq.tasks.count_sequences(nano_reads) == number_of_input_reads
    riak_out = f"{outprefix}.riak"
    riak_counts = run_read_it_and_keep([nano_reads], riak_out, "ont")
    reads_out = f"{riak_out}.reads.fastq.gz"
    # Some short reads/low % identity reads are not kept. This is ok. Their descriptions are:
    # MN908947.3,+strand,1807-2703 length=219 error-free_length=220 read_identity=86.83%
    # MN908947.3,-strand,1946-3027 length=78 error-free_length=80 read_identity=86.46%
    # MN908947.3,+strand,1721-2636 length=286 error-free_length=314 read_identity=73.39%
    # MN908947.3,+strand,1985-2995 length=36 error-free_length=38 read_identity=90.57%
    # junk_seq length=997 error-free_length=1022 read_identity=83.73%
    # MN908947.3,+strand,5-982 length=922 error-free_length=980 read_identity=76.17%

    number_of_output_reads = 26
    assert os.path.exists(reads_out)
    assert riak_debug_check_badreads_failed_reads_ok(f"{riak_out}.reads.debug")
    assert pyfastaq.tasks.count_sequences(reads_out) == number_of_output_reads
    assert riak_counts["Input reads file 1"] == number_of_input_reads
    assert riak_counts["Input reads file 2"] == 0
    assert riak_counts["Kept reads 1"] == number_of_output_reads
    assert riak_counts["Kept reads 2"] == 0
    subprocess.check_output(f"rm -rf {nano_reads} {reads_out} {riak_out}*", shell=True)

    illumina_read_prefix = f"{outprefix}.reads"
    simulate_illumina(illumina_read_prefix, ref=sim_reads_fasta)
    illumina_1 = f"{illumina_read_prefix}1.fq"
    illumina_2 = f"{illumina_read_prefix}2.fq"
    riak_counts = run_read_it_and_keep([illumina_1, illumina_2], riak_out, "illumina")
    number_of_reads = 65
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


def test_enumerate_names_unpaired():
    outprefix = "out.enumerate_names_unpaired"
    subprocess.check_output(f"rm -rf {outprefix}*", shell=True)
    read1 = COVID_REF_SEQ[0:100]
    read2 = COVID_REF_SEQ[200:300]
    test_reads_fa = f"{outprefix}.reads.fq"
    with open(test_reads_fa, "w") as f:
        print(">read1", read1, sep="\n", file=f)
        print(">read2", read2, sep="\n", file=f)
    riak_out = f"{outprefix}.riak"

    run_read_it_and_keep([test_reads_fa], outprefix, "illumina", enumerate_names=False)
    expect_lines = [">read1", read1, ">read2", read2]
    assert expect_lines == file_to_lines(f"{outprefix}.reads.fasta.gz")
    run_read_it_and_keep([test_reads_fa], outprefix, "illumina", enumerate_names=True)
    expect_lines = [">0", read1, ">1", read2]
    assert expect_lines == file_to_lines(f"{outprefix}.reads.fasta.gz")
    subprocess.check_output(f"rm -rf {outprefix}*", shell=True)


def test_enumerate_names_paired():
    outprefix = "out.enumerate_names_paired"
    subprocess.check_output(f"rm -rf {outprefix}*", shell=True)
    read1 = COVID_REF_SEQ[0:100]
    read2 = COVID_REF_SEQ[200:300]
    read3 = COVID_REF_SEQ[300:400]
    read4 = COVID_REF_SEQ[400:500]
    test_reads_fa1 = f"{outprefix}.reads.1.fq"
    test_reads_fa2 = f"{outprefix}.reads.2.fq"
    with open(test_reads_fa1, "w") as f:
        print(">read1", read1, sep="\n", file=f)
        print(">read2", read2, sep="\n", file=f)
    with open(test_reads_fa2, "w") as f:
        print(">read1", read3, sep="\n", file=f)
        print(">read2", read4, sep="\n", file=f)
    riak_out = f"{outprefix}.riak"

    fq_list = [test_reads_fa1, test_reads_fa2]
    run_read_it_and_keep(fq_list, outprefix, "illumina", enumerate_names=False)
    expect_lines1 = [">read1", read1, ">read2", read2]
    expect_lines2 = [">read1", read3, ">read2", read4]
    assert expect_lines1 == file_to_lines(f"{outprefix}.reads_1.fasta.gz")
    assert expect_lines2 == file_to_lines(f"{outprefix}.reads_2.fasta.gz")

    run_read_it_and_keep(fq_list, outprefix, "illumina", enumerate_names=True)
    expect_lines1 = [">0/1", read1, ">1/1", read2]
    expect_lines2 = [">0/2", read3, ">1/2", read4]
    assert expect_lines1 == file_to_lines(f"{outprefix}.reads_1.fasta.gz")
    assert expect_lines2 == file_to_lines(f"{outprefix}.reads_2.fasta.gz")
    subprocess.check_output(f"rm -rf {outprefix}*", shell=True)
