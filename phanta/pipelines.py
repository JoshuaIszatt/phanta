import os
from .functions import (
    interleave_reads
)


# _____________________________________________________PIPELINES


def assembly_pipeline(input_file, output_dir, isolate='unknown', log=False):
    pass

# _____________________________________________________BATCHES


def batch_assembly_pipeline(input_dir, output_dir, isolates="unknown"):
    os.makedirs(output_dir, exist_ok=True)
    for file in os.listdir(input_dir):
        path = os.path.join(input_dir, file)
        try:
            assembly_pipeline(path, output_dir, isolates)
        except Exception as e:
            print(f"ERROR {e}")
            continue
