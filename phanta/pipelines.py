import os
import pandas as pd
from .classes import PipelineError
from .functions import (
    configure_defaults,
    configure_log,
    find_paired_reads,
    find_interleaved_reads,
    interleave_reads,
    trim_interleaved_reads,
    deduplicate_reads,
    merge_short_reads,
    normalise_reads,
    fastqc,
    spades_assembly,
    read_mapping,
    extract_contig,
    generate_coverage_graph,
    checkv
)


# _____________________________________________________INPUT


def detect_reads(input_directory, config_file=None):
    # Set configuration
    if config_file is None:
        config = configure_defaults()
    else:
        config = configure_defaults(config_file)

    # Obtain read file extensions
    r1_ext = config["input"]["r1_ext"]
    r2_ext = config["input"]["r2_ext"]
    int_ext = config["input"]["interleaved_ext"]

    # Finding PE read pairs
    print(f"Using {r1_ext} and {r2_ext} to find read pairs")
    pe_reads = find_paired_reads(
        input_directory=input_directory,
        file_extension_1=r1_ext,
        file_extension_2=r2_ext,
        read_type="paired_end",
        exclude=[int_ext]
    )

    # Find interleaved PE reads
    print(f"Using {int_ext} to find interleaved reads")
    se_reads = find_interleaved_reads(
        input_directory=input_directory,
        file_extension=int_ext,
        read_type="interleaved",
        exclude=[r1_ext, r2_ext]
    )

    # Combining reads
    reads_list = pe_reads + se_reads

    # Return reads list
    if len(reads_list) == 0:
        print("No read pairs found")
        raise ValueError("No read pairs found")
    return reads_list


# _____________________________________________________PIPELINES


def assembly_pipeline(reads_class, output_dir, config_file=None, qc_only=False, stringent=False,
                      logger=None, logfile_location=None, logfile_configuration=None):
    # Setting logger if None has been passed
    if logger is None:
        if logfile_location is None:
            logfile_location = output_dir
        logger = configure_log(
            location=logfile_location,
            configuration=logfile_configuration
        )

    # Check booleans
    assert isinstance(qc_only, bool), 'Error: qc_only should be a boolean value.'
    assert isinstance(stringent, bool), 'Error: stringent should be a boolean value.'
    logger.info(f"Assembly pipeline settings: stringent={stringent}, qc_only={qc_only}")

    # Set pipeline configuration
    if config_file is None:
        config = configure_defaults()
        logger.info(f"Setting default pipeline configuration:")
        logger.info(config)
    else:
        config = configure_defaults(config_file)
        logger.info(f"Setting custom pipeline configuration: {config}")
        logger.info(config)

    # Check single element
    if isinstance(reads_class, (list, tuple)):
        if len(reads_class) == 1:
            reads_class = reads_class[0]
        else:
            e = "Batch sets of reads classes should be run using batch assembly pipeline"
            logger.error(e)
            raise PipelineError(e)

    # Attempt to use logger
    try:
        logger.info(f'Beginning assembly pipeline: {reads_class.name}')
    except Exception as e:
        raise Exception(f"Logging error: {e}")

    # Creating output directory
    out_dir = os.path.join(output_dir, f"{reads_class.name}")
    if os.path.exists(out_dir):
        e = "Output directory already exists"
        logger.error(e)
        raise PipelineError(e)
    else:
        os.makedirs(out_dir)

    # Interleave reads
    raw_reads = os.path.join(out_dir, "raw_reads.fastq.gz")
    if reads_class.type == "paired_end":
        logger.info("Interleaving paired end reads")
        interleave_reads(
            read_1=reads_class.read_1,
            read_2=reads_class.read_2,
            output_file=raw_reads,
            ram_mb=config['system']['RAM']
        )
    else:
        raw_reads = reads_class.read_1

    # Trimming reads
    trim_reads = os.path.join(out_dir, "trim_reads.fastq.gz")
    trim_interleaved_reads(
        reads=raw_reads,
        output_file=trim_reads,
        read_length=config['reads']['read_length'],
        ram_mb=config['system']['RAM'],
        trim_length=config['reads']['trim_length'],
        read_quality=config['reads']['read_quality'],
        minimum_length=config['reads']['minimum_length']
    )

    # Removing duplicates
    deduplicated_reads = os.path.join(out_dir, "deduplicated_reads.fastq.gz")
    deduplicate_reads(
        input_reads=trim_reads,
        output_reads=deduplicated_reads,
        ram_mb=config['system']['RAM']
    )

    # Merging reads
    merged_reads = os.path.join(out_dir, "merged_reads.fastq.gz")
    merge_short_reads(
        input_reads=deduplicated_reads,
        output_reads=merged_reads,
        ram_mb=config['system']['RAM'],
        error_correction=config['reads']['error_correction']
    )

    # Read normalisation
    normalised_reads = os.path.join(out_dir, "normalised_reads.fastq.gz")
    normalise_reads(
        input_reads=merged_reads,
        output_reads=normalised_reads,
        ram_mb=config['system']['RAM'],
        target_coverage=config['reads']['target_coverage']
    )

    # Reads QC
    out = os.path.join(out_dir, "Reads_QC")
    os.makedirs(out, exist_ok=False)
    reads = [raw_reads, trim_reads, deduplicated_reads, merged_reads, normalised_reads]
    for read in reads:
        fastqc(
            reads=read,
            output_directory=out
        )
    # todo Add multiqc

    # Checkpoint
    if qc_only:
        print(f"Finished reads QC for {reads_class.name}")
        return

    # Assembly 1
    out = os.path.join(out_dir, "SPAdes_initial")
    contigs = spades_assembly(
        input_reads=normalised_reads,
        output_directory=out,
        ram_mb=config['system']['RAM'],
        threads=config['system']['threads'],
        kmers=config['assembly']['kmers'],
    )

    # Read mapping (normalised - QC)
    out = os.path.join(out_dir, 'QC_read_mapping')
    basecov, covstats, scafstats, mapped, unmapped = read_mapping(
        contigs_fasta=contigs,
        reads=merged_reads,
        output_directory=out,
        ram_mb=config['system']['RAM'],
        keep_reads=True
    )

    # Extracting potential genomes
    df = pd.read_csv(scafstats, sep='\t')
    df = df[df['%unambiguousReads'] >= 90]
    if len(df) == 0:
        print("No contig with >90% reads mapped")
        raise Exception(f"No contig with >90% reads in {contigs}")
    contig_header = list(df['#name'])[0]

    put_genome = os.path.join(out_dir, 'putative_initial_genome.fasta')
    extract_contig(
        contigs_fasta=contigs,
        header=contig_header,
        output_file=put_genome,
        rename=f"{reads_class.name}_initial"
    )

    # Generating coverage graph
    generate_coverage_graph(
        header=contig_header,
        basecov=basecov,
        output_directory=out_dir
    )

    # CheckV
    if os.getenv('CHECKVDB'):
        logger.info(f"CHECKVDB variable detected, running analysis: {os.getenv('CHECKVDB')}")
        if not os.path.isdir(os.getenv('CHECKVDB')):
            logger.error(f"CheckV database does not exist: {os.getenv('CHECKVDB')}")
        else:
            try:
                outdir = os.path.join(out, 'CheckV')
                checkv(
                    contigs=contigs,
                    output_directory=out
                )
            except Exception as e:
                logger.warning(e)
    else:
        logger.info(f"No Checkv database detected, skipping...")

    # _________ SECOND SET

    if stringent:
        # Assembly 2
        out = os.path.join(out_dir, "SPAdes_mapped")
        contigs = spades_assembly(
            input_reads=mapped,
            output_directory=out,
            ram_mb=config['system']['RAM'],
            threads=config['system']['threads'],
            kmers=config['assembly']['kmers'],
        )

        # Read mapping
        out = os.path.join(out_dir, 'QC_read_remapping')
        basecov, covstats, scafstats, mapped, unmapped = read_mapping(
            contigs_fasta=contigs,
            reads=merged_reads,
            output_directory=out,
            ram_mb=config['system']['RAM'],
            keep_reads=False
        )

        # Extracting potential genomes
        df = pd.read_csv(scafstats, sep='\t')
        df = df[df['%unambiguousReads'] >= 90]
        if len(df) == 0:
            print("No contig with >90% reads mapped")
            raise Exception(f"No contig with >90% reads in {contigs}")
        contig_header = list(df['#name'])[0]

        put_genome = os.path.join(out_dir, 'putative_mapped_genome.fasta')
        extract_contig(
            contigs_fasta=contigs,
            header=contig_header,
            output_file=put_genome,
            rename=f"{reads_class.name}_mapped"
        )

        # Generating coverage graph
        generate_coverage_graph(
            header=contig_header,
            basecov=basecov,
            output_directory=out_dir
        )


# _____________________________________________________BATCHES


def batch_assembly_pipeline(input_dir, output_dir, config_file=None, qc_only=False, stringent=False,
                            logger=None, logfile_location=None, logfile_configuration=None):
    # Check booleans
    assert isinstance(qc_only, bool), 'Error: qc_only should be a boolean value.'
    assert isinstance(stringent, bool), 'Error: stringent should be a boolean value.'

    # Batch pipeline
    input_dir = os.path.expanduser(input_dir)
    output_dir = os.path.expanduser(output_dir)

    # Setting logger if None has been passed
    if logger is None:
        if logfile_location is None:
            logfile_location = output_dir
        logger = configure_log(
            location=logfile_location,
            configuration=logfile_configuration
        )

    # Setting in/out to log
    logger.info(f"input_dir: {input_dir}")
    logger.info(f"output_dir: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)
    count = 0

    # Check inputs
    if not os.path.isdir(input_dir):
        e = f'Input directory {input_dir} is not a directory'
        logger.error(e)
        raise ValueError(e)
    else:
        e = f'Batch assembly pipeline:'
        logger.info(e)

    # Detecting reads
    reads = detect_reads(
        input_directory=input_dir,
        config_file=config_file
    )
    logger.info(f"Found {len(reads)} reads")
    for read in reads:
        count = count+1
        logger.info(f"Running input file #{count}")
        try:
            assembly_pipeline(
                reads_class=read,
                output_dir=output_dir,
                config_file=config_file,
                qc_only=qc_only,
                stringent=stringent,
                logger=logger
            )
        except Exception as e:
            logger.error(f'PipelineError: Reads: {read.name},  {e}')
            continue

    # Finish batch pipeline
    logger.info("Batch pipeline complete")
