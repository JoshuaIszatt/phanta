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
