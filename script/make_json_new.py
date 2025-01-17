import json
import argparse
import os
import sys
import pdb
# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("--fastq_dir", help="Please input the full path to the folder containing fastq.gz files")
args = parser.parse_args()

# Ensure the provided path exists
assert args.fastq_dir is not None, "Please provide the path to the folder containing fastq.gz files"
if not os.path.exists(args.fastq_dir):
    print("The input folder does not exist")
    sys.exit()

# Variables
fastq_dir = args.fastq_dir
fastq_files = []
files_by_sample = {}

# Traverse subdirectories to collect fastq files
for dirpath, _, filenames in os.walk(fastq_dir):
    for filename in filenames:
        if filename.endswith(".fastq.gz"):
            fastq_files.append(os.path.abspath(os.path.join(dirpath, filename)))
#pdb.set_trace()
# Extract sample names and organize files
#samples = set(f.split('_')[0] for f in map(os.path.basename, fastq_files))
samples = set(f.rsplit('_1.fastq.gz', 1)[0].rsplit('_2.fastq.gz', 1)[0] for f in map(os.path.basename, fastq_files))

for sample in samples:
    read1_files = sorted(f for f in fastq_files if sample in f and "_1.fastq.gz" in f)
    read2_files = sorted(f for f in fastq_files if sample in f and "_2.fastq.gz" in f)
    files_by_sample[sample] = {"1": read1_files, "2": read2_files}

#pdb.set_trace()
# Output the JSON file
output_json = "../samples.json"
with open(output_json, "w") as f:
    json.dump(files_by_sample, f, indent=4, sort_keys=True)

print(f"Sample details have been written to {output_json}")

