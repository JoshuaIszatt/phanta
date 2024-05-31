# Phanta
Configurable short read assembly pipeline for phages:
![Phage pipeline](pipeline.png)
**Figure 1:** Rough phage assembly pipeline, dotted line indicates where processing begins.

## Install dependencies
1. Download the yml environment file:
```
wget https://anaconda.org/JoshIszatt/phanta/2024.05.31.130738/download/phanta.yml
```

2. Create the environment
```
conda env create --file phanta.yml
```

3. Activate the environment:
```
conda activate phanta
```

4. Install phanta pipeline:
```
pip install phanta
```

5. Optional: Setup the checkv database 
```
checkv download_database /path/to/checkv-db
```

```
export CHECKVDB=/path/to/checkv-db
```

## Usage
Open a python terminal and enter:
```py
import phanta
dir(phanta)
```

Detect reads:
This will return a list of Reads class objects that can be passed individually or as a batch to assembly pipelines
```py
import phanta
reads = phanta.detect_reads('path_to_input_directory/')
```

Single phage assembly:
```py
import phanta
reads = phanta.detect_reads('path_to_input_directory/')
phanta.assembly_pipeline(reads[0], 'output_directory/')
```

Batch phage assembly:
```py
import phanta
reads = phanta.detect_reads('path_to_input_directory/')
phanta.batch_assembly_pipeline(reads, 'output_directory/')
```

## Dependencies:
  - python>=3
  - checkv==1.0.3
  - biopython==1.83
  - bbmap==39.06
  - pandas==2.2.1
  - matplotlib==3.8.4
  - spades==3.15.5
  - fastqc==0.12.1
  - multiqc==1.22.1
  - pilon==1.24
  - packaging
  - deprecated

### Configuration
These are the defaults, amend as required and pass to the pipeline functions
```json
{
 "system": {
  "RAM": 24000,
  "threads": 8
 },
 "input": {
  "interleaved_ext": "_interleaved.fastq",
  "r1_ext": "_R1.fastq.gz",
  "r2_ext": "_R2.fastq.gz"
 },
 "reads": {
  "read_length": 150,
  "trim_length": 10,
  "minimum_length": 100,
  "read_quality": 30,
  "error_correction": true,
  "target_coverage": 200
 },
 "assembly": {
  "kmers": "55,77,99,127"
 }
}
```