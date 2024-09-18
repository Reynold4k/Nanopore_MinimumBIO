# Pipeline 2

Once you have executed the initial analysis with your gene expression counts and identified potential hits, it’s time to further investigate them using AlphaFold for protein structure prediction. Below are detailed instructions for running the script while ensuring the relevant paths are properly set.

### Setting Up Paths
# Create a new conda environment
```bash
conda create -n alphafold_env python=3.9
conda activate alphafold_env
```

# Install canu, if it is not installed:
```bash
wget https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.bz2
tar -xvjf canu-2.2.Linux-amd64.tar.bz2
export PATH=$PATH:path/to/your/canu-2.2/bin
```
# Step A: Extracting potential hits into a new folder

#Edit file path in the script:

```bash
#Path to your bam files
PARENT_DIR="/path/to/your/experiments/R3/step2"
#Path to the new folder: potential hits
OUTPUT_DIR="/path/to/your/experiments/potential_hit"
#Annotation files
ANN_FILE="path/to/your/reference/Homo_sapiens.GRCh38.110.gtf"

Hit_LOCATIONS=$(awk -v gene_name="your_interested potential hits name" '$3 == "exon" && $0 ~ gene_name {print "chr"$1":"$4"-"$5}' "$ANN_FILE")

```

```bash

bash ./a_extract_reads.pbs

```

If running without error, you will see:

![image](https://github.com/user-attachments/assets/6114dda2-d07d-456c-9ac2-c96fe67e32f4)


# Step B: Error Correction, trimming, unitigging, Assembly

#Edit file path in the script:

OUTPUT_DIR="/path/to/your/experiments/potential_hit"

```bash

bash ./b_Correction_Assembly.pbs

```

It won't take a long time, but you need to wait till the contigs files are created!!!

You must wait until seeing the contigs file:
![image](https://github.com/user-attachments/assets/9c6e2bd7-c371-4376-84dd-e84f345fe0f4)


Then execute the translation script, but still need to change the path first:

contig_fasta_file = "/path/to/your/corrected.contigs.fasta"
output_protein_file = "/path/to/your/protein_sequences_long.fasta"

```bash
python3 translation.py
```

If succeed, you will see:
![image](https://github.com/user-attachments/assets/e30b575e-6021-4724-9a16-80a87febc39c)

Make a blast online to check whether it is still the potential hits that you extracted in the Step A via: 

https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_SPEC=GeoBlast&PAGE_TYPE=BlastSearch


# Step C: AlphaFold

It is not easy to run AlphaFold on Katana, please debug and solve every error comes with the execution of this command, I've listed several known errors below, good luck.

```bash

cd path/to/your/alphafold
pip install -r requirements.txt

```

```bash

python3 run_alphafold.py \
    --fasta_paths=/srv/scratch/z3546698/true/Small_Molecule/JQ1/T7MB-1/231104/potential_hit/FASTQ/canu_out/protein_sequences_long.fasta \
    --output_dir=/srv/scratch/z3546698/true/Small_Molecule/JQ1/T7MB-1/231104/potential_hit/alphafold_output \
    --max_template_date=2024-09-16 \
    --model_preset=monomer \
    --data_dir=/data/bio/alphafold \
    --uniref90_database_path=/data/bio/alphafold/uniref90/uniref90.fasta \
    --uniref30_database_path=/data/bio/alphafold/uniref30/UniRef30_2021_03 \
    --mgnify_database_path=/data/bio/alphafold/mgnify/mgy_clusters_2022_05.fa \
    --template_mmcif_dir=/data/bio/alphafold/pdb_mmcif/mmcif_files \
    --obsolete_pdbs_path=/data/bio/alphafold/pdb_mmcif/obsolete.dat \
    --use_gpu_relax=True \
    --bfd_database_path=/data/bio/alphafold/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --pdb70_database_path=/data/bio/alphafold/pdb70/pdb70

```

If you beat the final boss and win, you'll see:
![image](https://github.com/user-attachments/assets/9b0579bd-2c10-4f23-885c-ffaf2e07e874)

These pdb files are what you may want to check as the result of AlphaFold, you can view them via Linux:

```bash
conda install -c conda-forge pymol

pymol your_protein.pdb

```

If you want to run it faster, please add GPU into calculation of AlphaFold:

```bash
pip install --upgrade "jax[cuda]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
```

### when showing error:

ImportError: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.26' not found (required by /srv/scratch/z3546698/true/ls/envs/alphafold_env/lib/python3.9/site-packages/scipy/linalg/_matfuncs_sqrtm_triu.cpython-39-x86_64-linux-gnu.so)

```bash
conda install -c conda-forge libgcc-ng

ls $CONDA_PREFIX/lib | grep libstdc++.so

export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

strings $CONDA_PREFIX/lib/libstdc++.so.6 | grep GLIBCXX

```

Check the output for GLIBCXX_3.4.26, confirming that the symbol is part of the new library version.





ImportError: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.26' not found (required by /srv/scratch/z3546698/true/ls/envs/alphafold_env/lib/python3.9/site-packages/scipy/linalg/_matfuncs_sqrtm_triu.cpython-39-x86_64-linux-gnu.so)

Solution:
```bash

cd $CONDA_PREFIX/lib
ls -l | grep libstdc++
ln -s libstdc++.so.6.0.33 libstdc++.so.6
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

# Once the symlink has been successfully created, run the following command to verify that the GLIBCXX 
# version is included:
# This should list your support for all GLIBCXX versions of libstdc++.
strings $CONDA_PREFIX/lib/libstdc++.so.6 | grep GLIBCXX
```

# For Kalign, please refer:
https://github.com/TimoLassmann/kalign?tab=readme-ov-file

```bash

mkdir build
cd /srv/scratch/z3546698/true/alphafold/kalign-3.4.0/build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/local/kalign
make
make install
export PATH="$HOME/local/kalign/bin:$PATH"
```
