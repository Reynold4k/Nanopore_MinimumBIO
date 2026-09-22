# Beginner setup guide

This guide covers basic Linux/macOS terminal usage, Windows Subsystem for Linux (WSL), and installing R/RStudio. It is migrated from the original project READMEs for users who are new to command-line bioinformatics.

## Windows users: Ubuntu on WSL

If you are running Windows, install Ubuntu from the Microsoft Store or from <https://ubuntu.com/>.

### Understanding file paths

- Ubuntu uses forward slashes (`/`) instead of backslashes (`\`).
  - Example: `C:\Users\YourName\Documents` becomes `/home/YourName/Documents`.
- File and directory names are case-sensitive.
  - `Documents`, `documents`, and `DOCUMENTS` are different directories.
- Your personal files live in the home directory, `/home/YourName` (similar to `C:\Users\YourName`).
- Windows drives are mounted under `/mnt`:
  - `C:` is `/mnt/c`
  - `D:` is `/mnt/d`
- Files and folders starting with a dot (`.`) are hidden. Press `Ctrl + H` in the file manager to show them.
- Permissions are stricter than on Windows. You may need `chmod` or `chown` for some tasks.

### Installing basic tools

Open the Ubuntu terminal and run:

```bash
sudo apt update
sudo apt install python3 r-base wget build-essential
```

If an installation command fails, try prefixing it with `sudo`.

### Useful Ubuntu tips

- Backup your data before making major changes. Tools like `rsync` are useful.
- Software is usually installed via the `apt` package manager or the Ubuntu Software Center.
- Ubuntu does not rely on file extensions as strictly as Windows; it often identifies files by their content.

## macOS users: using Terminal

Terminal is built into macOS.

- Open Terminal from `Applications > Utilities > Terminal`.
- The command line interface lets you type commands to perform tasks.

## Basic terminal commands

These commands work on both Linux and macOS.

### Navigating the file system

```bash
pwd              # print current directory
ls               # list files and directories
ls -la           # list with hidden files and details
cd directory_name    # change directory
cd ..            # go up one level
```

### Managing files and directories

```bash
mkdir directory_name    # create a directory
touch file_name         # create an empty file
cp source destination   # copy a file or directory
mv source destination   # move or rename a file
rm file_name            # delete a file
rm -r directory_name    # delete a directory and its contents
```

### Editing files

```bash
nano file_name    # edit a file in the terminal
```

Use `Ctrl + O` to save, `Ctrl + X` to exit.

### Permissions and ownership

```bash
chmod permissions file_name    # change file permissions
chown user file_name           # change file owner
```

### Searching and finding files

```bash
grep search_term file_name          # search inside a file
find directory -name file_name      # find a file by name
```

### System information

```bash
df -h        # disk usage
free -h      # memory usage (Linux)
top          # running processes
```

### Useful shortcuts

- Press `Tab` to auto-complete file and directory names.
- Press `Up` / `Down` to scroll through previous commands.
- Press `Ctrl + C` to cancel a running command.

## Installing R and RStudio

### Step 1: Install R

1. Visit the CRAN website: <https://cran.r-project.org/>.
2. Choose your operating system (Windows, macOS, or Linux).
3. Download and run the installer for the latest R version.

On Ubuntu you can also install R from the terminal:

```bash
sudo apt install r-base
```

### Step 2: Install RStudio

1. Visit the RStudio download page: <https://posit.co/download/rstudio-desktop/>.
2. Download RStudio Desktop (Open Source License) for your operating system.
3. Run the installer.

### Step 3: Getting started with R and RStudio

1. Open RStudio.
2. The interface has several panes: console, script editor, environment/history, and files/plots/packages/help.
3. Create a new script: `File > New File > R Script`.
4. Type `print("Hello, World!")` and press `Ctrl + Enter` (macOS: `Cmd + Enter`) to run it.

### Installing R packages

In the R console or RStudio:

```r
install.packages("ggplot2")
```

For Bioconductor packages you first need `BiocManager`:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("edgeR", "rtracklayer"))
```

### Running an R script from the terminal

```bash
Rscript my_analysis.R
```

### Running an R script in RStudio

1. Open the script in RStudio.
2. Highlight the lines you want to run and press `Ctrl + Enter` (`Cmd + Enter` on macOS).
3. To run the whole script, click `Source` or press `Ctrl + Shift + S` (`Cmd + Shift + S` on macOS).

Outputs appear in the Console pane.

## Reference files

Typical reference inputs for this workflow are:

- Reference genome in FASTA format, e.g. `hg38.fa`.
- Gene annotation in GTF format, e.g. `hg38.ensGene.gtf`.

### Downloading the human reference genome

```bash
sudo apt install wget bwa

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

bwa index hg38.fa
```

### Downloading a GTF annotation

```bash
wget http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz
```

You can download GTF files for other species from Ensembl (<ftp://ftp.ensembl.org/pub/release-110/gtf/>) or UCSC (<http://hgdownload.soe.ucsc.edu/downloads.html>).

### Downloading a human protein database for BLAST

```bash
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# Optional: keep only human entries
zcat uniprot_sprot.fasta.gz | grep -A 1 '^>.*OS=Homo sapiens' > human_proteins.fasta

# Install BLAST+ and build the database
sudo apt-get install ncbi-blast+
makeblastdb -in human_proteins.fasta -dbtype prot -out human_proteome_db
```

### BLAST+ open-files limit

On some systems BLAST+ may fail with a memory-mapping error because the default open-files limit is too low. Increase it with:

```bash
ulimit -n unlimited
# or, if that is rejected:
ulimit -n 65536
```

## Next steps

Once your environment is set up, follow the quick-start instructions in `README.md` and edit your workflow config file.
