### Comprehensive Pipeline for Nanopore Biopanning data

**Author**: Chen Zhu  
**Email**: [z3546698@ad.unsw.edu.au]  
**Date**: 04/Sep/2024


### If you don't have a personalized or targeted reference genome file, please go through:

## <span style="color: red;">Pipeline1 & 2</span>

### If you actually have a personalized or targeted reference genome file, please go through:

## <span style="color: green;">Pipeline1.5 & 2</span>

## Please note, it would be different procedures between <span style="color: red;">Pipeline1</span> and <span style="color: green;">Pipeline1.5</span>, please carefully check your settings.

## Prerequisites

Ensure your system meets the following requirements before running the pipeline:

- **Operating System**: Linux or compatible environment.

# For Windows Users Ubuntu from microsoft store is highly recommended:
- You could also obtain this software via: https://ubuntu.com/

![image](https://github.com/user-attachments/assets/9b400296-fc30-4dbd-929c-7c9bdafaf438)

## If you didn't install python ever on your PC, install NOW after you read the Linux instructions below!!
Type the command below in your linux portal
```bash
sudo apt install python3
```

## Guide for Windows Users Transitioning to Ubuntu

## Understanding File Paths
- **Path Delimiters:**  
  Unlike Windows, which uses backslashes (`\`) for file paths, Ubuntu utilizes forward slashes (`/`).  
  Example: `C:\Users\YourName\Documents` becomes `/home/YourName/Documents` in Ubuntu.

- **Case Sensitivity:**  
  File and directory names in Ubuntu are case-sensitive.  
  `Documents`, `documents`, and `DOCUMENTS` are considered different directories.

## Home Directory
- On Ubuntu, your personal files are stored in the home directory, typically located at `/home/YourName`. 
  This is similar to `C:\Users\YourName` on Windows.

## Accessing Drives
- Windows drives (C:, D:, etc.) are mounted in the `/mnt` or `/media` directory in Ubuntu.  
  For example, your C: drive might be accessible under `/mnt/c`.

## Hidden Files and Directories
- Files and folders prefixed with a dot (`.`) in Ubuntu are hidden by default.  
  To view them in the file manager, press `Ctrl + H`.

## Permissions
- Ubuntu enforces file permissions more strictly than Windows.  
  You might need to modify permissions using commands like `chmod` or change the owner with `chown` for certain tasks.

## Using the Terminal
- The terminal is a powerful tool in Ubuntu, used for a variety of tasks.  
  Get familiar with basic commands like `ls` (list), `cd` (change directory), and `cp` (copy).

## Installing Software
- Unlike Windows, software in Ubuntu is often installed via package managers like `apt` (Advanced Package Tool).  
  You might also use software repositories like the Ubuntu Software Center.

## File Extensions
- Ubuntu does not rely on file extensions to identify file types as strictly as Windows.  
  It often determines the file type by its content.

## Backup Your Data
- Before making any major changes, always ensure your data is backed up.  
  Tools like `rsync` can be invaluable for maintaining backups in Ubuntu.


# For Mac Users Click the Launchpad icon in the Dock, type Terminal in the search field, then click Terminal.
# In the Finder , open the /Applications/Utilities folder, then double-click Terminal.

    
## Beginner's Guide to Using Terminal and Command Line on macOS

## Introduction to Terminal
- **Accessing Terminal:**  
  Terminal is a built-in application in macOS, found under `Applications > Utilities > Terminal`.

- **Basic Interface:**  
  Terminal provides a command line interface where you can type commands to perform various tasks.

## Navigating the File System
- **Current Directory:**
  - Use `pwd` to print the current working directory.

- **Listing Files and Directories:**  
  - Use `ls` to list files and directories.  
  - Use `ls -la` to include hidden files and detailed information.

- **Changing Directory:**  
  - Use `cd [directory_name]` to navigate to a different directory.  
  - Use `cd ..` to go up one directory level.

## Managing Files and Directories
- **Creating Directories:**  
  - Use `mkdir [directory_name]` to create a new directory.

- **Creating Files:**  
  - Use `touch [file_name]` to create a new, empty file.

- **Copying Files:**  
  - Use `cp [source] [destination]` to copy files or directories.

- **Moving/Renaming Files:**  
  - Use `mv [source] [destination]` to move or rename files.

- **Deleting Files and Directories:**  
  - Use `rm [file_name]` to delete files.  
  - Use `rm -r [directory_name]` to delete directories and their contents.

## Editing Files
- **Using Nano Editor:**  
  - Use `nano [file_name]` to edit files directly in Terminal.

## Permissions
- **Changing Permissions:**  
  - Use `chmod [permissions] [file_name]` to change file or directory permissions.

- **Changing Ownership:**  
  - Use `chown [user] [file_name]` to change file ownership.

## Searching and Finding Files
- **Search with grep:**  
  - Use `grep [search_term] [file_name]` to search for a term within a file.

- **Finding Files:**  
  - Use `find [directory] -name [file_name]` to search for files by name.

## Useful Tips
- **Auto-Completion:**  
  - Use the `Tab` key to auto-complete file and directory names.

- **Command History:**  
  - Use the `Up` and `Down` arrow keys to cycle through command history.

- **Canceling Commands:**  
  - Use `Ctrl + C` to cancel an ongoing command or process.


**Reference and Annotation Files**:
- Reference genome in FASTA format (e.g., `hg38.fa`).
- Gene annotation file in GTF format (e.g., `hg38.ensGene.gtf`).

# Beginner's Guide to Essential Linux Commands


## Introduction to the Terminal
- **Accessing the Terminal:**  
  The terminal is a command line interface to interact with your Linux system. You can open it from your system's applications menu or by using a keyboard shortcut (usually `Ctrl + Alt + T`).

## Navigating the File System
- **Current Directory:**
  - Use `pwd` (print working directory) to display your current directory path.

- **Listing Files and Directories:**  
  - Use `ls` to list directory contents.  
  - Use `ls -la` for detailed information, including hidden files.

- **Changing Directory:**  
  - Use `cd [directory_name]` to navigate between directories.  
  - Use `cd ..` to move up one level.

## Managing Files and Directories
- **Creating Directories:**  
  - Use `mkdir [directory_name]` to create a new directory.

- **Creating Files:**  
  - Use `touch [file_name]` to create an empty file.

- **Copying Files:**  
  - Use `cp [source] [destination]` to copy files or directories.

- **Moving/Renaming Files:**  
  - Use `mv [source] [destination]` to move or rename files.

- **Deleting Files and Directories:**  
  - Use `rm [file_name]` to delete files.  
  - Use `rm -r [directory_name]` to delete directories and their contents.

## Editing Files
- **Using Nano Editor:**  
  - Use `nano [file_name]` to edit files within the terminal.

## Permissions and Ownership
- **Changing Permissions:**  
  - Use `chmod [permissions] [file_name]` to modify file permissions.

- **Changing Ownership:**  
  - Use `chown [user] [file_name]` to change file ownership.

## Searching and Finding Files
- **Search with grep:**  
  - Use `grep [search_term] [file_name]` to search for a term within files.

- **Finding Files:**  
  - Use `find [directory] -name [file_name]` to locate files by name.

## System Information
- **Check Disk Usage:**  
  - Use `df -h` to display disk space usage.

- **Check Memory Usage:**  
  - Use `free -h` to display memory usage.

- **View Running Processes:**  
  - Use `top` or `htop` to view active processes.

## Useful Tips
- **Auto-Completion:**  
  - Use the `Tab` key for auto-completing commands and file/directory names.

- **Command History:**  
  - Use the `Up` and `Down` arrow keys to scroll through previously used commands.

- **Canceling Commands:**  
  - Use `Ctrl + C` to stop an ongoing command or process.



### Please follow the instructions and check the results after each step.
## Good Luck
![image](https://github.com/user-attachments/assets/33e84a7b-7fde-481b-a0ee-999fbe9a18d2)

