import pandas as pd


# Config error
class ConfigFileError(Exception):
    """Raised when the configuration file is invalid or corrupt."""
    pass


# Pipeline error
class PipelineError(Exception):
    """General exception for pipeline errors."""

    def __init__(self, message="An error occurred in the pipeline."):
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'PipelineError: {self.message}'


# Reads QC class
class Reads(object):
    def __init__(self, name, read_type, read_1, read_2):
        self.name = name
        self.type = read_type
        self.read_1 = read_1
        self.read_2 = read_2


# Genome assembly class
class Phage(object):
    def __init__(self, name, tag, batch, fasta_hash, reads_hash, hostname):
        self.name = name
        self.tag = tag
        self.passage = batch
        self.fasta_hash = fasta_hash
        self.reads_hash = reads_hash
        self.hostname = hostname


# Mapping files class (bbmap)
class BBpath:
    def __init__(self, dir_path):
        self.aqhist = f"{dir_path}/aqhist.tsv"
        self.basecov = f"{dir_path}/basecov.tsv"
        self.covstats = f"{dir_path}/covstats.tsv"
        self.gchist = f"{dir_path}/gchist.tsv"
        self.lhist = f"{dir_path}/lhist.tsv"
        self.qhist = f"{dir_path}/qhist.tsv"
        self.scafstats = f"{dir_path}/scafstats.tsv"
        self.mapped = f"{dir_path}/mapped.fastq.gz"
        self.unmapped = f"{dir_path}/unmapped.fastq.gz"

    def find_genomes(self, mincov=90, minlen=4000):
        # Reading data
        df = pd.read_csv(self.scafstats, sep='\t')
        df2 = pd.read_csv(self.covstats, sep='\t')
        # Minlength
        parse = list((df2[df2['Length'] >= minlen])['#ID'])
        df = df[df['#name'].isin(parse)]
        # Coverage
        df = df[df['%unambiguousReads'] >= mincov]
        # Return header
        try:
            contig_header = df.loc[df['%unambiguousReads'].idxmax(), '#name']
            return contig_header
        except ValueError:
            return None


# Mapping files class (bbmap)
class CheckvPath:
    def __init__(self, dir_path):
        self.complete_genomes = f"{dir_path}/complete_genomes.tsv"
        self.completeness = f"{dir_path}/completeness.tsv"
        self.contamination = f"{dir_path}/contamination.tsv"
        self.proviruses = f"{dir_path}/proviruses.fna"
        self.quality_summary = f"{dir_path}/quality_summary.tsv"
        self.viruses = f"{dir_path}/viruses.tsv"

    def find_genomes(self, mincomplete=90, minlen=4000):
        df = pd.read_csv(self.quality_summary, sep='\t')
        # Minlength
        df = df[df['contig_length'] >= minlen]
        # Mincomplete
        df = df[df['completeness'] >= mincomplete]
        # Return headers
        contig_headers = list(df['contig_id'])
        return contig_headers

    def return_data(self):
        pass  # todo Concatenate and return all checkv data
