# model_validator
model_validator is a tool to validate profile HMMs and adjust their cutoff score

##   Instalation
model_validator does not need to be installed. The user should only download the model_validator.py file.

## Requirements

## Usage
```
python model_validator.py -conf <configuration_file>
python model_validator.py -i <hmm_file> -model_type <'l'|'long'|'s'|'short'> -db_type 'short' -cell_org_db_frag <fasta_file> <optional parameters>
python model_validator.py -i <hmm_file> -model_type <'l'|'long'|'s'|'short'> -db_type 'long' -cell_org_db <fasta_file> -vir_db <directory> <optional parameters>
```
### Mandatory parameters:
```
-cell_org_db <fasta_file>	 	Cellular organisms database, mandatory if -db_type = 'long' or 'l'.
-cell_org_db_frag <fasta_file>	 	Cellular organisms fragmented database, mandatory if -db_type = 'short' or 's'.
-conf <file>	 	Configuration file (parameters in format 'parameter=value', one parameter per line and parameter names equal to those on the command line).
-db_type <'l'|'long'|'s'|'short'>	 	Type of database. If -db_type = 'short' or 's', -cell_org_db_frag is mandatory. Else if -db_type = 'long' or 'l', -cell_org_db and -vir_db are mandatory
-i <hmm_file>	 	Profile HMMs to be validated.
-model_type <'l'|'long'|'s'|'short'> Type of model.
-vir_db <directory>	 	Directory containing viral sequences, mandatory if -db_type = 'long' or 'l'.
```

### Optional parameters:
```
-lh_lr <string>	 	For -model_type = 'long' or 'l' and -db_type = 'long' or 'l', hmm-prospector cell_org run parameters.
-lh_sr <string>	 	For -model_type = 'long' or 'l', hmm-prospector virus run parameters.
-out <name>	 	Name of output directory (default: hmm_validated).
-pd <decimal>	 	Minimum accepted percentage detection rate of viral sequences (default = 80) after applying new cutoff scores to the HMMs profile, must be between 0 and 100.
-pt <decimal>	 	Maximum percentage ratio (default = 80), must be between 0 and 100. Given similarity results against cell organisms (OrgCellDB) and viral sequences (VirusDB), highest score obtained for OrgCellDB must be lower than pt% of the lowest score obtained against Virus, otherwise the model is not accepted.
-sh_lr <string>	 	For -model_type = 'short' or 's' and -db_type = 'long' or 'l', hmm-prospector cell_org run parameters.
-sh_rv <string>	 	For -model_type = 'short' or 's', hmm-prospector virus run parameters.
-sh_sr <string>	 	For -model_type = 'short' or 's' and -db_type = 'short' or 's', hmm-prospector cell_org run parameters.
``` 

## Contact
To report bugs, to ask for help and to give any feedback, please contact Arthur Gruber (argruber@usp.br) or Giuliana L. Pola (giulianapola@usp.br).

## Versions
### 1.0.9B
- Creation of invalidated.csv CSV file that contains invalidated model, protein, taxon and family
- Add line break at the end of generated templates

### 1.0.8B
- Change in the rules for creating the new cutoff score, generated less invalidation of the models

### 1.0.8
- Fixed the error when getting the username to put in the log

### 1.0.7
- Delete validHMMs file
- Invalid selected.hmm becomes invalid.hmm
- Valid selected.hmm becomes valid.hmm
- Log file ends with a summary

### 1.0.6
- Addition of fasta file path in viral database (-virdb) when taxon is subfamily
- Fixed hmm-prospector default parameters for short sequences (db_type=short)
- Changed the rank of those who were not matched from “Mismatch” to “No match”

### 1.0.5
- Embedding the get_rank function in the get_vir_db function
- Using the Bio.Entrez.efetch (db="Taxonomy") package to get the taxon classification and family

### 1.0.4
- Replacement of the requests module with the Bio.Entrez package when getting the taxonomic classification of the taxon (Entrez.efetch)
- Adding more details about each analyzed model

### 1.0.3
- Correction when getting the lowest viral score (V) 