# cabsflex_restraints
This repository contains a Python script designed to generate new Category restraints for CABS-flex.

This project is written in Python 3, and we recommend using version 3.9. It requires the following libraries to run: `numpy`, `pandas`, and `biopandas`. We recommend setting up the environment using either **Conda** or **pip**.

## Installation

### Option 1: Using Conda

1. **Clone the repository** (if applicable):
```
git clone https://github.com/kwroblewski7/cabsflex_restraints.git
cd cabsflex_restraints
```

2. **Create the Conda environment** using the `environment.yml` file:
```
conda env create -f environment.yml
```

3. **Activate the environment**:
```
conda activate cabsflex_restraints
```

This will create a Conda environment named `cabsflex_restraints` with Python 3.9 and install `numpy`, `pandas`, and `biopandas`.

### Option 2: Using pip

1. **Clone the repository** (if applicable):
```
git clone https://github.com/kwroblewski7/cabsflex_restraints.git
cd cabsflex_restraints
```

2. **Create a virtual environment** with Python 3.9:
```
python3.9 -m venv cabsflex_restraints
```

3. **Activate the virtual environment**:

```
source cabsflex_restraints/bin/activate
```

4. **Install the required packages** using `requirements.txt`:
```
pip install -r requirements.txt
```

This will create a Python 3.9 virtual environment named `cabsflex_restraints` with `numpy`, `pandas`, and `biopandas` installed.

## Usage

After setting up the environment with either method, you can run the scripts or Jupyter notebooks in this project.

### What You'll Need to Provide:
- **PDB File (`-i` or `--pdb_file`)**  
  Path to your PDB file.<br>
  Example: `-i 6avj.pdb`

- **pLDDT File (`-p` or `--plddt_file`)**  
  This can be one of two things:
    - A path to a .json file with pLDDT scores.
    - `'bf'` if you want to use B-factor values as pLDDT scores directly.  
    Example: `-p 6avj_scores.json` or `-p bf`

- **Secondary Structure File (`-s` or `--ss_file`)**  
  Path to the file containing the protein's secondary structure information in old fixed column DSSP format.<br>
  For generating DSSP output, you can use the tool available at https://github.com/PDB-REDO/dssp. <br>
  Example: `-s 6avj_dssp.txt`

### Optional Settings (You Can Skip These):
- **Working Directory (`-w` or `--work_dir`)**  
  Directory where the script looks for files and saves results. Default is the current folder.<br>
  Example: `-w example`

- **Output File Name (`-o` or `--output_file`)**  
  Name for the output file with the generated restraints. Default is `ca_restraints_file.txt`.<br>
  Example: `-o 6avj_ca_restraints_file.txt`

### Example Commands:
``` 
python3 create_restraints.py -i 6avj.pdb -p 6avj_scores.json -s 6avj_dssp.txt -w example -o 6avj_ca_restraints_file.txt
```
``` 
python3 create_restraints.py -i 6avj.pdb -p bf -s 6avj_dssp.txt -w ./example -o 6avj_ca_restraints_file.txt
```
## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.