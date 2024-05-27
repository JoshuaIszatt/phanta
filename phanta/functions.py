import os
import subprocess
from Bio import SeqIO
from deprecated import deprecated
from .decorators import experimental

def interleave_reads(read_1, read_2, output_file, ram_mb=20000):
    command = [
        "reformat.sh",
        f"-in={read_1}",
        f"-in2={read_2}",
        f"-Xmx{ram_mb}m",
        f"-out={output_file}"
    ]
    try:
        subprocess.run(command, check=True)
        print("Reads interleaved successfully")
    except Exception as e:
        print(f"Read trimming failed {e}")
        raise e


@deprecated(reason="Reads should be interleaved using interleave_reads", version="pre release")
def trim_pe_reads(read_1, read_2, output_file, ram_mb, read_length, trim_length, read_quality, minimum_length):
    tl = 0 + trim_length
    tr = read_length - trim_length
    command = [
        "bbduk.sh",
        f"-Xmx{ram_mb}m",
        "tpe",
        "tbo",
        f"in1={read_1}",
        f"in2={read_2}",
        f"out={output_file}",
        f"ftl={tl}",
        f"ftr={tr}",
        f"minavgquality={read_quality}",
        f"minlength={minimum_length}"
    ]
    try:
        subprocess.run(command, check=True)
        print(f"Reads trimmed successfully, Q:{read_quality}")
    except Exception as e:
        print(f"Read trimming failed {e}")
        raise


def trim_interleaved_reads(reads, output_file,
                           read_length=150, ram_mb=20000, trim_length=10,
                           read_quality=30, minimum_length=0):
    tl = 0 + trim_length
    tr = read_length - trim_length
    command = [
        "bbduk.sh",
        f"-Xmx{ram_mb}m",
        "tpe",
        "tbo",
        f"in={reads}",
        f"out={output_file}",
        f"ftl={tl}",
        f"ftr={tr}",
        f"minavgquality={read_quality}",
        f"minlength={minimum_length}"
    ]
    try:
        subprocess.run(command, check=True)
        print(f"Reads trimmed successfully, Q:{read_quality}")
    except Exception as e:
        print(f"Read trimming failed {e}")
        raise


def convert_bam_to_fasta(input_reads_bam, output_reads_fasta, ram_mb=20000):
    command = [
        'reformat.sh',
        f'in={input_reads_bam}',
        f'out={output_reads_fasta}',
        f'-Xmx{ram_mb}m'
    ]
    try:
        subprocess.run(command, check=True)
    except Exception as e:
        print(f"BAM to FA conversion failed: {e}")
        raise
    return os.path.abspath(output_reads_fasta)


def deduplicate_reads(input_reads, output_reads, ram_mb=20000):
    command = [
        "dedupe.sh",
        f"-Xmx{ram_mb}m",
        f"in={input_reads}",
        f"out={output_reads}"
    ]
    try:
        subprocess.run(command, check=True)
        print(f"Reads deduplicated successfully")
    except Exception as e:
        print(f"Read deduplication failed {e}")
        raise


def normalise_reads(input_reads, output_reads, ram_mb=20000, target_coverage=200):
    command = [
        "bbnorm.sh",
        f"-Xmx{ram_mb}m",
        "min=5",
        f"target={target_coverage}",
        f"in={input_reads}",
        f"out={output_reads}"
    ]
    try:
        subprocess.run(command, check=True)
        print(f"Reads normalised to {target_coverage}x")
    except Exception as e:
        print(f"Read normalisation failed {e}")
        raise


@experimental(description="Unsure if this should be part of the pipeline, needs testing")
def merge_short_reads(input_reads, output_reads, ram_mb=20000, error_correction=True):
    command = [
        "bbmerge.sh",
        f"-Xmx{ram_mb}m",
        f"in={input_reads}",
        f"out={output_reads}",
        f"vstrict"
    ]
    if error_correction:
        command.extend(["ecco", "mix"])
    try:
        subprocess.run(command, check=True)
        print(f"bbmerge ran successfully")
        return output_reads
    except Exception as e:
        print(f"bbmerge failed {e}")
        raise


def fastqc(reads, output_directory):
    """
    @param reads:
    @param output_directory:
    @return:
    """
    command = [
        "fastqc",
        f"{reads}",
        "-o",
        f"{output_directory}"
    ]
    try:
        subprocess.run(command, check=True)
        print("Reads QC success")
    except subprocess.CalledProcessError:
        print("Reads QC failed")


def spades_assembly(input_reads, output_directory, ram_mb=20000, threads=8):
    input_reads_path = os.path.abspath(input_reads)
    output_path = os.path.abspath(output_directory)
    ram_gb = int(ram_mb / 1000)
    command = [
        "spades.py",
        "-t", f"{threads}",
        "-m", f"{ram_gb}",
        "--only-assembler",
        "--careful",
        "-k", "55,77,99,127",
        "-o", f"{output_path}",
        "--12", f"{input_reads_path}"
    ]
    try:
        subprocess.run(command, check=True)
        print(f"SPAdes genome assembly finished successfully")
    except Exception as e:
        print(f"SPAdes genome assembly failed {e}")
        raise


def sam_sort_index(input_reads, output_directory):
    """
    @param input_reads: Must be SAM formatted reads
    @param output_directory: Output location
    @return:
    """
    output_bam = os.path.join(output_directory, "mapped_sorted.bam")
    view_command = [
        "samtools", "view", "-bS", "-F4", f"{input_reads}"
    ]
    sort_command = [
        "samtools", "sort", "-", "-o", output_bam
    ]
    # Execute the first command and pipe the output to the second command
    with subprocess.Popen(view_command, stdout=subprocess.PIPE) as proc1:
        with subprocess.Popen(sort_command, stdin=proc1.stdout) as proc2:
            proc1.stdout.close()
            proc2.communicate()
    print(f"Sorted BAM file created at: {output_bam}")
    index_command = [
        "samtools", "index", output_bam
    ]
    subprocess.run(index_command, check=True)
    return output_bam


def pilon_polish(genome_fasta, reads_bam, output_directory):
    command = [
        "pilon",
        "--genome", genome_fasta,
        "--frags", reads_bam,
        "--outdir", output_directory,
        "--changes"
    ]
    try:
        subprocess.run(command, check=True)
    except Exception as e:
        print(f"Pilon failed {e}")
        raise


def read_mapping(contigs_fasta, reads, output_directory, ram_mb=20000, mapped_sam=False):
    covstats = os.path.join(output_directory, "covstats.tsv")
    basecov = os.path.join(output_directory, "basecov.tsv")
    scafstats = os.path.join(output_directory, "scafstats.tsv")
    qhist = os.path.join(output_directory, "qhist.tsv")
    aqhist = os.path.join(output_directory, "aqhist.tsv")
    lhist = os.path.join(output_directory, "lhist.tsv")
    gchist = os.path.join(output_directory, "gchist.tsv")
    command = [
        "bbmap.sh",
        f"-Xmx{ram_mb}m",
        f"ref={contigs_fasta}",
        f"in={reads}",
        f"covstats={covstats}",
        f"basecov={basecov}",
        f"scafstats={scafstats}",
        f"qhist={qhist}",
        f"aqhist={aqhist}",
        f"lhist={lhist}",
        f"gchist={gchist}",
        "nodisk"
    ]
    if mapped_sam:
        mapped = os.path.join(output_directory, "mapped.sam")
        command.append(f"out={mapped}")
    else:
        mapped = False
    try:
        subprocess.run(command, check=True)
        print(f"Reads mapped")
        return basecov, covstats, scafstats, mapped
    except Exception as e:
        print(f"Read mapping failed: {e}")
        raise Exception(f"Read mapping failed {e}")


def split_mapped_reads(putative_genome, reads, output_directory, ram_mb=20000):
    mapped = os.path.join(output_directory, "mapped.fastq.gz")
    unmapped = os.path.join(output_directory, "unmapped.fastq.gz")
    command = [
        "bbmap.sh",
        f"-Xmx{ram_mb}m",
        f"ref={putative_genome}",
        f"in={reads}",
        f"outm={mapped}",
        f"outu={unmapped}",
        "nodisk"
    ]
    try:
        subprocess.run(command, check=True)
        print(f"Reads mapped and split")
        return mapped, unmapped
    except Exception as e:
        print(f"Read mapping failed: {e}")
        raise


def extract_contig(contigs_fasta, header, output_file, rename=None):
    """
    @param contigs_fasta:
    @param header:
    @param output_file:
    @param rename:
    @return:
    """
    with open(contigs_fasta, 'r') as handle:
        entries = SeqIO.parse(handle, 'fasta')
        with open(output_file, 'w') as textfile:
            for entry in entries:
                if header in entry.id:
                    if rename:
                        entry.id = rename
                    try:
                        SeqIO.write(entry, textfile, 'fasta')
                        break
                    except Exception as e:
                        print(f"Could not extract {output_file}")
                        raise Exception(f"Could not extract {output_file}: {e}")


@experimental(description="Needs improvements")
def split_multifasta(input_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    with open(input_file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            output_file = os.path.join(output_dir, f"{record.id}.fasta")
            with open(output_file, "w") as output_handle:
                SeqIO.write(record, output_handle, "fasta")


def checkv(contigs, output_directory):
    command = [
        "checkv", "end_to_end",
        f"{contigs}",
        f"{output_directory}"
    ]
    try:
        subprocess.run(command, check=True)
        print("CheckV successful")
    except subprocess.CalledProcessError:
        print("CheckV failed")


def snippy(reference_genome, reads, output_directory):
    command = [
        "snippy",
        "--ref", reference_genome,
        "--outdir", output_directory,
        "--peil", reads,
        "--report",
        "--cleanup"
    ]
    try:
        subprocess.run(command, check=True)
        print("snippy successful")
    except subprocess.CalledProcessError:
        print("snippy failed")
    pass
