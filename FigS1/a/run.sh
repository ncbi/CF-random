#!/bin/bash
module load alphafold2/2.2.0

run_singularity \
--model_preset=monomer \
--fasta_paths=$PWD/${1}.fasta \
--model_config=$PWD/config.py \
--max_template_date=1960-01-01 \
--use_precomputed_msas \
--output_dir=$PWD
