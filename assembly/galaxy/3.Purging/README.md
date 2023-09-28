# Purge Haplotigs/Overlaps and QC
Galaxy Workflows for removing haplotigs and contig overlaps from assemblies based on read depth using Purge_Dups.

## Hap1/Hap2 assemblies
The workflow takes a trimmed HiFi reads collection, Hap1/Hap2 contigs, and the values for transition parameter and max coverage depth (calculated from WF1) to run Purge_Dups. It produces purged Hap1 and Hap2 contigs assemblies, and runs all the QC analysis (gfastats, BUSCO, and Merqury).
![PurgeHap1Hap2](pics/PurgeHap1Hap21908.png)

## Pri/Alt assemblies
\[in preparation]
