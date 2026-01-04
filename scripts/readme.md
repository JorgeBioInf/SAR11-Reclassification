# What does each script perform?
The following scripts are adapted to a specific workflow. They may need to be adapted for related purposes.

Requirements are specified in each script's code. 

- 1) `genomes_download.py`: retrieves and filters SAR11 genomes according to completeness and contamination thresholds.
- 2) `fastANI`
  3) `ANI_grouping.py`: clusters genomes based on ANI.
  4) `ANI_byClade.ipynb`: computes inter-clade and intra-clade ANI, and performs visual representations.
  5) `GTDB_processer.py`: retrieves GTDB classifications for each genome.
  6) `GTDB-tk` for the remaining genomes.
  7) `prodigal_analysis.sh`: executes `prodigal` for each source group.
  8) `clade_dividing.sh`: divides test genomes according to their clade.
  9) `prodigal_group_anaysis.sh`: executes `prodigal` for each clade group.
  10) `CSF_sources_analysis.sh`: performs `ConSpeciFix` analysis among source groups.
  11) `CSF_clades_analysis.sh`: performs `ConSpeciFix` analysis between clade groups and all source groups.
  12) `PopCOGenT`
  13) `summary_table.py`: summarizes genome-wise classification in a single table.
  14) `ByClade_Table_Grouping.ipynb`: groups genome-wise classification by clade.
