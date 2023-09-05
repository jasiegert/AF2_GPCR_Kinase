# A Python interface to model user-defined funtional and structural features of Kinases and GPCRs with AlphaFold2.

This repository is an expansion of our previous [colabfold](https://github.com/sokrypton/ColabFold)-based work ["Sampling alternative conformational states of transporters and receptors with AlphaFold2"](https://elifesciences.org/articles/75751) by Diego del Alamo, Davide Sala, Hassane S. Mchaourab, and Jens Meiler. All the functionalities have been kept. To have a general overview, please read the README.md at https://github.com/delalamo/af2_conformations. Here, our previous workflow is extended with the aim of predicting a user-defined conformational state of GPCR and Kinase with minimal effort. We also introduced known features like using a custom pdb template or print the predicted pTM score for each model. 

### How to use the code in this repository

Before importing the code contained in the `scripts/` folder, the user needs to install the AlphaFold source code and download the parameters to a directory named `params/`. Additional Python modules that must be installed include [Numpy](https://numpy.org/), [Requests](https://docs.python-requests.org/en/latest/), and [Logging](https://abseil.io/docs/python/guides/logging).

The scripts can be imported and used out-of-the-box to fetch multiple sequence alignments and/or templates of interest. Note that the `max_msa_clusters` and `max_extra_msa` options can be provided to reduce the size of the multiple sequence alignment. If these are not provided, the networks default values will be used. Additional options allow the number of recycles, as well as the number of loops through the recurrent Structure Module, to be specified. In addition, ptm can be enabled to print pTM score as a suffix of model name. 

## Predicting a user-defined GPCR functional state

To predict a specific activation state of a GPCR target, the pdbs list must contain one of the following string in the first position ("Inactive", "Active", "Intermediate", "G protein", "Arrestin"). The script will retrieve templates in the annotated functional state from GPCRdb.org. Template PDBs can be excluded simply by adding PDB IDs without chain the ID specified to the list. Example: ["G protein", "7FII"] to predict the active state of your target by using G protein bound templates but excluding 7FII. Which PDB ids have been used to bias the prediction can be retrieved from the log file (example by typing 'grep PDBS example.log').
For GPCR targets, PDBs of proteins from a given GPCR subfamily can be excluded. The subfamily can be fetched from the GPCRdb given a protein name. For example for LSHR (GPCRdb "family" entry: 001_003_003_002), the subfamily would be "001_003_003" and passing it along as exclude_gpcr_subfamily would exclude the closely related Glycoprotein hormone receptors FSHR (001_003_003_001), LSHR (001_003_003_002) and TSHR (001_003_003_001).
Templates can also be randomized for each model. 
Below, an example of outputting info into example.log and predict 50 models of LSHR by using the best 4 templates determined in active state, but excluding all LSHR PDBs released and all PDBS of proteins from the same subfamily (LSHR, FSHR, TSHR).  

```python
from AF2_GPCR_Kinase.scripts import mmseqs2
import multiprocessing
import logging
logging.basicConfig(filename='example.log', level=logging.DEBUG) # print log with debug level

# Jobname for reference
jobname = 'lshr_gprot_4t'

# Amino acid sequence. Whitespace and inappropriate characters are automatically removed
sequence = ("YDFLRVLIWLINILAIMGNMTVLFVLLTSRYKLTVPRFLMCNLSFADFCMGLYLLLIASVDSQTKGQYYNHAIDWQTGSGCSTAGFFTVFASELSVYTLTVITLERWHTITYAIHLDQKLRLRHAILIMLGGWLFSSLIAMLPLVGVSNYMKVSICFPMDVETTLSQVYILTILILNVVAFFIICACYIKIYFAVRNPELMATNKDTKIAKKMAILIFTDFTCMAPISFFAISAAFKVPLITVTNSKVLLVLFYPINSCANPFLYAIFTKTFQRDFFLLLSKFGCC")

# State annotation followed by PDB IDs to be excluded.
pdbs = ["Active", "7FII", "7FIG", "7FIH", "7FIJ"]
# Fetch subfamily from GPCRdb (e.g. 001_003_003 for lshr -> Class A, Protein receptors, Glycoprotein hormone receptors)
target_gpcr_subfamily = mmseqs2.get_subfamily("lshr_human")

#parameters
max_msa_clusters = 32 # Number of sequence clusters
max_extra_msa = 64   # Number of extra sequences not clustered
max_recycles = 1     # Number of neural network iterations
n_struct_module_repeats = 8 # Number of the structure module iterations
n_models = 50   # Number of models to be predicted
model_id = -1   # Which AF neural network. -1 = Randomize
model_params = -1 # Which AF neural network parameters. -1 = Randomize
ptm = True    # Print pTM value before .pdb of each model
_rank = 1   # Number assigned to the first predicted model
remove_msa_for_template_aligned = False  # Remove the genetic information for regions already covered by templates. Copied from Heo L. et al., DOI: 10.1002/prot.26382.
n_templates = 4 # Number of templates to be used

# Initializes the Runner object that queries the MMSeqs2 server
mmseqs2_runner = mmseqs2.MMSeqs2Runner( jobname, sequence, n_templates = n_templates )

# Fetches the data and saves to the appropriate directory
# Pass along activation state in pdbs[0], structures to exclude in pdbs[1:] and subfamily to exclude
a3m_lines, template_path = mmseqs2_runner.run_job(templates = pdbs , exclude_gpcr_subfamily = target_gpcr_subfamily )

from AF2_GPCR_Kinase.scripts import predict

for i in range( n_models ):
  model_name = str(jobname + "_" + str(_rank) + ".pdb")
  
  # Optionally, templates can be shuffled within the list of PDBs passing filters. 
  # Uncomment line below to enable templates randomization.
  # template_path = mmseqs2_runner.shuffle_templates()
  
  # Run a prediction with templates
  predict.predict_structure_from_templates( sequence, model_name, a3m_lines, template_path=template_path,  model_id=model_id, max_msa_clusters=max_msa_clusters, max_extra_msa=max_extra_msa, max_recycles=max_recycles, n_struct_module_repeats=n_struct_module_repeats, ptm=ptm, remove_msa_for_template_aligned=remove_msa_for_template_aligned )
  
  #Two alternatives to predict 1. without templates or 2. with local pdb as a template
  
  # 1. Run a prediction without templates 
  predict.predict_structure_no_templates( sequence, model_name, a3m_lines, model_id=model_id, max_msa_clusters=max_msa_clusters, max_extra_msa=max_extra_msa, max_recycles=max_recycles, n_struct_module_repeats=n_struct_module_repeats, ptm=ptm)
         
# 2. Run a prediction with a local pdb template. 
  predict.predict_structure_from_custom_template( sequence, model_name, a3m_lines, template_pdb="pdb_file",  model_id=model_id, max_msa_clusters=max_msa_clusters, max_extra_msa=max_extra_msa, max_recycles=max_recycles, n_struct_module_repeats=n_struct_module_repeats, ptm=ptm,      remove_msa_for_template_aligned=remove_msa_for_template_aligned)
  
  _rank += 1
```


## Predicting user-defined structural features of kinases

Similar to predicting user-defined GPCRs functional states, users can force AF2 to retrieve kinase templates matching three structural feature values from KLIFS. Allowed structural properties are: 1. DFG 2. aC_helix 3. Salt-bridge (KIII.17 and EαC.24). Allowed values for the corresponding property 1. DFG: out, in, out-like, all 2. aC_helix: out, in, all 3. Salt-bridge: yes, no, all. Following the example above, the format must be a three members list in the first position of the pdbs list. Besides sequence and jobname, what to change in the script above to predict models with templates in DFG=out is reported below.

```python
# structural properties of kinase templates to be used
kinase_temps = ["out", "all", "all" ] # Format is [DFG, aC_helix, salt_bridge]

# kinases annotation list followed by PDB IDs to be excluded
pdbs = [kinase_temps, "6N3O", "6N3L", "6N3N", "7QQ6", "7QWK"]
```
### Introducing mutations into MSA

There is also functionality to introduce mutations (e.g. alanines) across the entire MSA to remove the evolutionary evidence for specific interactions (see [here](https://www.biorxiv.org/content/10.1101/2021.11.29.470469v1) and [here](https://twitter.com/sokrypton/status/1464748132852547591) on why you would want to do this). This can be achieved as follows:

```python
# Define the mutations and introduce into the sequence and MSA
residues = [ 41,42,45,46,56,59,60,63,281,282,285,286,403,407 ]
muts = { r: "A" for r in residues }
mutated_msa = util.mutate_msa( a3m_lines, muts )
```
### Known issues

Here is a shortlist of known problems that we are currently working on:
* The MMSeqs2 server queries the PDB70, rather than the full PDB. This can cause some structures to be missed if their sequences are nearly identical to those of other PDB files.
* Multimer prediction is not currently supported.
* Custom MSAs are not currently supported.
* Additional annotations of both GPCRs and kinases are not currently supported.

If you find any other issues please let us know in the "issues" tab above.

### Citations

If the code in this repository has helped your scientific project, please consider citing our papers:

```bibtex
@article {Sala2022.12.11.519936,
author = {Sala, Davide and Meiler, Jens},
title = {Biasing AlphaFold2 to predict GPCRs and Kinases with user-defined functional or structural properties},
elocation-id = {2022.12.11.519936},
year = {2022},
doi = {10.1101/2022.12.11.519936},
publisher = {Cold Spring Harbor Laboratory},
URL = {https://www.biorxiv.org/content/early/2022/12/11/2022.12.11.519936},
eprint = {https://www.biorxiv.org/content/early/2022/12/11/2022.12.11.519936.full.pdf},
journal = {bioRxiv}
}
@article {10.7554/eLife.75751,
article_type = {journal},
title = {Sampling alternative conformational states of transporters and receptors with AlphaFold2},
author = {del Alamo, Diego and Sala, Davide and Mchaourab, Hassane S and Meiler, Jens},
editor = {Robertson, Janice L and Swartz, Kenton J and Robertson, Janice L},
volume = 11,
year = 2022,
month = {mar},
pub_date = {2022-03-03},
pages = {e75751},
citation = {eLife 2022;11:e75751},
doi = {10.7554/eLife.75751},
url = {https://doi.org/10.7554/eLife.75751},
journal = {eLife},
issn = {2050-084X},
publisher = {eLife Sciences Publications, Ltd},
}
```

