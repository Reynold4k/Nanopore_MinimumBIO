
# Bioinformatics Data Processing Pipeline and Gene Coverage Analysis

## Overview of the Bioinformatics Data Processing Pipeline

This pipeline script is designed to process high-throughput sequencing data by executing a series of automated data processing steps on specified folders containing FASTQ files. The main functionalities include:

1. **Merging and Trimming FASTQ Files**: It merges all FASTQ files within each sample directory into a single file and trims the merged file using Porechop to remove low-quality sequences and adapter sequences.

2. **Quality Control**: Utilizes NanoPlot to assess the quality of the trimmed sequences, ensuring the reliability of the data.

3. **Alignment and Coverage Analysis**: Leverages Minimap2 to align the trimmed sequences to a reference genome, generating BAM files and calculating read coverage, which is then output to `coverage.txt`.

4. **Position Extraction and BED File Generation**: Extracts relevant gene positions from the coverage data and generates a BED file for easier visualization and downstream analysis.

5. **Sequence Extraction**: Based on the generated BED file, extracts specified regions from the reference sequences and saves them in FASTA format for further analysis and applications.

This script provides an automated workflow for bioinformatics research, significantly increasing data processing efficiency while ensuring accuracy and consistency in the results.

### Step 1: Modifying and Executing `pipeline1.5.pbs`

To run the `pipeline1.5.pbs` script, follow these steps:

1. **Open the Script**: Open the file `pipeline1.5.pbs` in a text editor.

2. **Modify Paths**:
    - Change **Line 31**: Locate the line that defines `FASTQ_FOLDER` and replace it with the path to your folder containing the FASTQ files. For example:
      ```bash
      FASTQ_FOLDER="/path/to/your/fastq/files"
      ```
      This will be the folder path like below:
      
      ![image](https://github.com/user-attachments/assets/9cfb1d87-747c-4e3d-8f47-ec78c661b40b)

    - Change **Line 36**: Locate the line that defines `REFERENCE` and update the path to your reference file. For example:
      ```bash
      REFERENCE="/path/to/your/reference/file"
      ```
    - If you haven't index your reference fasta file, please index it first using the command below:
      ```bash
      
      bwa index your_reference_fasta
      
      ```

3. **Set Up the Environment**: At the beginning of the script, add commands to load the required modules and set up a Python virtual environment. Follow these steps:

    - **Load Python Module**: Load the desired version of Python.
    ```bash
    module load python/3.10.8
    ```

    - **Create the Virtual Environment**: Create a virtual environment in your home directory.
    ```bash
    python3 -m venv /home/your_username/environments/pipeline1.5_env
    ```

    - **Activate the Virtual Environment**: Activate the created virtual environment.
    ```bash
    source /home/your_username/environments/pipeline1.5_env/bin/activate
    ```

    - **Install Required Packages**: Install the necessary Python packages with pip.
    ```bash
    pip3 install pandas matplotlib

    #Check if the installation is successful:
    which pandas
    which matplotlib
    ```

4. **Execute the Script**: Once you have modified the paths and set up the environment, save the file and submit it to your job scheduler. For example, use the following command:
    ```bash
    qsub pipeline1.5.pbs
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

To run the `Compare.pbs` script, follow these steps:

1. **Open the Script**: Open the file `Compare.pbs` in a text editor.

2. **Modify Paths**:
    - Change **Line 25**: Locate the line that defines `EXPERIMENTAL_FOLDER` and replace it with the path to your experimental group's FASTQ parent folder. For example:
      ```bash
      EXPERIMENTAL_FOLDER="/path/to/your/exp/fastq"
      ```
    - Change **Line 26**: Locate the line that defines `CONTROL_FOLDER` and replace it with the path to your control group's FASTQ parent folder. For example:
      ```bash
      CONTROL_FOLDER="/path/to/your/control/fastq"
      ```
    - Change **Line 27**: Locate the line that defines `REFERENCE` and update the path to your reference file. For example:
      ```bash
      REFERENCE="/path/to/your/path/to/reference"
      ```

3. **Execute the Script**: Once you have modified the paths, save the file and submit it to your job scheduler. For example, use the following command:
    ```bash
    qsub Compare.pbs
    ```

4. **Analysis**: Once you have executed the step3 above successfully, analyse it using the R script (Analysis.R) following command:
   ```bash
    conda activate bioenv
    mkdir R_library
    module load r/4.4.0
   ```
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


