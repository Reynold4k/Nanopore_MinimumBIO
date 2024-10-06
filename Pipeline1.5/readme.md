
# Bioinformatics Data Processing Pipeline and Gene Coverage Analysis

## Overview of the Bioinformatics Data Processing Pipeline

This pipeline script is designed to process high-throughput sequencing data by executing a series of automated data processing steps on specified folders containing FASTQ files. The main functionalities include:

1. **Merging and Trimming FASTQ Files**: It merges all FASTQ files within each sample directory into a single file and trims the merged file using Porechop to remove low-quality sequences and adapter sequences.

2. **Quality Control**: Utilizes NanoPlot to assess the quality of the trimmed sequences, ensuring the reliability of the data.

3. **Alignment and Coverage Analysis**: Leverages Minimap2 to align the trimmed sequences to a reference genome, generating BAM files and calculating read coverage, which is then output to `coverage.txt`.

4. **Position Extraction and BED File Generation**: Extracts relevant gene positions from the coverage data and generates a BED file for easier visualization and downstream analysis.

5. **Sequence Extraction**: Based on the generated BED file, extracts specified regions from the reference sequences and saves them in FASTA format for further analysis and applications.

These scripts provides an automated workflow for bioinformatics research, significantly increasing data processing efficiency while ensuring accuracy and consistency in the results.

### Overall of all the scripts you would go through here is shown in order below:

pipeline1.5.pbs -> Compare.pbs -> Analysis.R

### Step 1: Modifying and Executing `pipeline1.5.pbs`

Function introduction: Pipeline1.5.pbs is the script which goes through the working path you specified in the script, you need to edit and modify the working path in the script before running, please note, if you're using Katana, you would need to load several modules coz Katana is not able to remember the modules you used before——not like your local PC, but don't worry the stable versions of these modules have been comprehensive definied in the very begining of each script you would go through.

In conclusion, what is going on in the pipeline1.5.pbs is to go through every rounds in your specified experimental folder and load all the fastq first, then it would be merging all the fastq files into one to accelerate the analysing process, more details are available in the start of the script.

The software I recommend you to use is VS Code, you could download it via:https://code.visualstudio.com/.

So let's start!

To run the `pipeline1.5.pbs` script, follow these steps:

1. **Open the Script**: Open the file `pipeline1.5.pbs` in a text editor, if you're using Katana, you could also open the script in the VS Code.

2. **Modify Paths**:
    - Change **Line 31**: Locate the line that defines `FASTQ_FOLDER` and replace it with the path to your folder containing the FASTQ files. For example:
      ```bash
      #You can add as many folder as you want as long as they all have the same reference file:
        FASTQ_FOLDER=(
            "/mnt/d/Bait_Glue/CRBN/glue/TON/240427"
            "/mnt/d/Bait_Glue/CRBN/MB014/TON/230827")
       ```
      This will be the folder path like below:
      
      ![image](https://github.com/user-attachments/assets/9cfb1d87-747c-4e3d-8f47-ec78c661b40b)

    - Change **Line 36**: Locate the line that defines `REFERENCE` and update the path to your reference file. For example:
      ```bash
      REFERENCE="/path/to/your/reference/file"
      ```
    - If you haven't index your reference fasta file, please index it first using the command below:
      ```bash
      module load bwa
      bwa index your_reference_fasta
      
      ```
![image](https://github.com/user-attachments/assets/aeb9ee35-ca29-4b82-abae-15250df71707)


3. **Set Up the Environment**: At the beginning of the script, add commands to load the required modules and set up a Conda environment. Follow these steps:

    - **Load Python Module**: Load the desired version of Python. (This is only applied when using Katana)
    ```bash
    module load python/3.10.8
    ```

    - **Create the Conda Environment**: Create a Conda environment in your home directory.
    ```bash
    conda create -n my_conda python=3.10
    ```
    
    If you're unsure about your current directory path, using the below command:

    ```bash
    pwd
    ```
    
    Use the environment name from the directory where you want it:

    ```bash
    conda create --prefix /your/current/directory/my_conda python=3.10
    ```

    - **Activate the Conda Environment**: Activate the created Conda environment.
    ```bash
    conda activate my_conda
    ```
    Or, if you used a specific path to create the environment, you can activate it with:
    ```bash
    conda activate /your/current/directory/my_conda
    ```
    
    - **Install Required Packages**: Install the necessary Python packages with pip.
    ```bash
    pip3 install pandas matplotlib

    #Check if the installation is successful:
    which pandas
    which matplotlib
    ```

5. **Execute the Script**: Once you have modified the paths and set up the environment, save the file and submit it to your job scheduler. For example, use the following command:

    If you didn't book any CPU, then you need to create another pbs script having the contents below:

```bash
   #!/bin/bash
    #how many CPU and Memory you would like to use when running:
    #PBS -l select=1:ncpus=1:mem=4gb
    #how much time you would like to use when running:
    #PBS -l walltime=12:00:00
    #how to let you know when to finish:
    #PBS -M your.name.here@unsw.edu.au
    #PBS -m ae
    #PBS -j oe

    #don't change it
    cd $PBS_O_WORKDIR
    #Make sure the name of the script you would to run is correct:
    ./pipeline1.5.pbs
```

    ```bash
    #If you didn't book any CPU:
    qsub run.pbs
    ```
### Highly recommended!!!
    ```bash
    #If you did book any CPU, directly run:
    bash ./pipeline1.5.pbs
    ```

### Notes:
- **Virtual Environment Activation**: After activating the virtual environment, the shell prompt will change, indicating that the environment is active. You can now use Python and pip to install and manage packages within this environment.
- Replace `your_username` in the paths with your actual username.
- This process ensures that you are using an isolated environment for your Python scripts, making dependency management easier.

### Expected Output

The output of the `pipeline1.5.pbs` script will include:
- A merged and trimmed FASTQ file.
- Quality assessment report from NanoPlot.
- BAM files generated through the alignment process.
- Coverage data saved in `coverage.txt`.
- A BED file containing relevant gene positions.
- FASTA format files for specified regions.

These outputs will be stored in the specified output directories.

---


## Overview of the Gene Coverage Analysis Script `Compare.pbs`

The main purpose of the `Compare.pbs` script is to analyze the gene coverage between experimental and control groups, calculate differential coverage, and generate a corresponding horizontal bar plot for visualization. The specific steps are as follows:

1. **Load Required Modules**: Use the `module load` command to load the necessary Python and computational tools.

2. **Define Paths**: Specify the paths for the experimental group, control group, and reference files.

3. **Find the Latest Experimental Directory**: Retrieve the latest experimental and control round directories based on the folder structure.

4. **Check Coverage Files**: Confirm that the coverage files exist.

5. **Read Coverage Data**: Read the coverage data for the control and experimental groups into associative arrays.

6. **Calculate Differences**: Compare the coverage between the control and experimental groups, calculate the difference, and write the results to an output file.

7. **Generate Visualization Charts**: Use Python to create a horizontal bar plot that displays the coverage for the experimental group.

### Step 2: Modifying and Executing `Compare.pbs`

Function introduction: As you may find out yourself, if you don't really have an annotation gtf file rather than your personalized reference fasta, what you could do is to do the similar calculation on those data you could have at this stage. So Compare.pbs compare the difference of the coverage files you generated during step2 between experimental and control groups when you specify them into the correct directory.

So let's continue!

To run the `Compare.pbs` script, follow these steps:

1. **Open the Script**: Open the file `Compare.pbs` in a text editor, or whatever you can to edit the text, if you're a person of super natural, please use the command:

Before running the script, you would be required to fully understand the functions of Nano command through, including how to edit and save it:
https://ioflood.com/blog/nano-linux-command/#:~:text=To%20use%20the%20nano%20command,will%20create%20it%20for%20you.

   ```bash
    Nano Compare.pbs
   ```

2. **Modify Paths**:
    - Change **Line 25**: Locate the line that defines `EXPERIMENTAL_FOLDER` and replace it with the path to your experimental group's FASTQ parent folder. For example:
      ```bash
      EXPERIMENTAL_FOLDER="/path/to/your/exp/fastq"
      ```
    - Change **Line 26**: Locate the line that defines `CONTROL_FOLDER` and replace it with the path to your control group's FASTQ parent folder. For example:
      ```bash
      CONTROL_FOLDER="/path/to/your/control/fastq"
      ```

![image](https://github.com/user-attachments/assets/5a2f7537-9bba-4d75-9802-6b22a9ee7991)


3. **Execute the Script**: Once you have modified the paths, save the file and submit it to your job scheduler. For example, use the following command:
```bash
   #!/bin/bash
    #how many CPU and Memory you would like to use when running:
    #PBS -l select=1:ncpus=1:mem=4gb
    #how much time you would like to use when running:
    #PBS -l walltime=12:00:00
    #how to let you know when to finish:
    #PBS -M your.name.here@unsw.edu.au
    #PBS -m ae
    #PBS -j oe

    #don't change it
    cd $PBS_O_WORKDIR
    #Make sure the name of the script you would to run is correct:
    ./compare.pbs
```

    ```bash
    #If you didn't book any CPU:
    qsub run.pbs
    ```
### Highly recommended!!!
    ```bash
    #If you did book any CPU, directly run:
    bash ./compare.pbs
    ```
   
### Step3 Analysis: Once you have executed the step3 above successfully, analyse it using the R script (Analysis.R) following command:
   
The following command illustrates how to intiate the environment and make a new R library, this is important!!!

## What you what to modify:

 1.Change the line 12,13,14 to your/newly_created_library_for_R
```bash
#Please replace the z3546698 with your own zid when running, for example:
mkdir /srv/scratch/z3546698/R_library

```

And then change the folder path of the Pipeline1.R in the lines 12,13,14 with a text editor of Analysis.R

# 2.Change the line 17 EXPERIMENTAL_FOLDER to your/exp/fastq parent path

```r
# Set the experimental folder path
EXPERIMENTAL_FOLDER <- "/mnt/d/Bait_Glue/VHL/MB012/TON/230827"
```

# 3.Change the line 19 id_mapping to your/id_mapping_file, you could save the id_mapping file in your previous reference folder if you want to:

```r
id_mapping <- read.table("path/to/your/reference/idmapping_2024_10_01.tsv", 
                         header = TRUE, 
                         sep = "\t", 
                         stringsAsFactors = FALSE, 
                         fill = TRUE, 
                         quote = "",  
                         comment.char = "") 
```


Then run through the whole R scripts and check the result plots.

```bash

module load r/4.4.0

Rscript script.R
```
![image](https://github.com/user-attachments/assets/1b616d91-ed8a-4f04-a898-c80ba191a47a)

### After you modified the path directory in the Analysis.R, run the R script through:
   
   ```bash
    Rscript Analysis.R
   ``` 
### Expected Output

The output of the `Compare.pbs` script will include:
- A `differential_coverage.txt` file, which contains the gene IDs, control coverage, experimental coverage, and the calculated differences.
- A horizontal bar plot image saved in a specified output directory, visually representing the differences in coverage between the experimental and control groups.

### Interpreting the Output Results

- **differential_coverage.txt**: This file contains the following columns:
  - **Gene**: The identifier for each gene analyzed.
  - **Control_Coverage**: The coverage value for the control group.
  - **Experimental_Coverage**: The coverage value for the experimental group.
  - **Difference**: The difference in coverage values between the experimental and control groups.
![image](https://github.com/user-attachments/assets/0179be74-56a4-4577-8207-5dbf1be378d2)

![volcano_plot](https://github.com/user-attachments/assets/d70d1959-62a4-426c-a840-0de4fa5e4c46)


