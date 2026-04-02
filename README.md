# BIGapp-PanelHub

A central metadata repository of marker-panel definitions for downstream workflows in [BIGapp](https://github.com/Breeding-Insight/BIGapp) and [BIGr](https://github.com/Breeding-Insight/BIGr).

## `.botloci` Files

- **Purpose:** Enumerate all marker IDs for which bottom-strand probes were designed.  
- **Workflow Role:** During MADC-to-VCF conversion, these files tell BIGapp/BIGr exactly which sequences to extract as reverse complements.  
- **Format:** Plain-text lists—one marker ID per line (optionally namespaced by panel or species).

## `_lut.csv` Files

- **Purpose:** Provide mappable CloneID with the botloci files in case MADC CloneID doesn't match; provide REF and ALT bases for speedup madc2vcf_targets; provide indels position and length. 
- **Workflow Role:** Used by all `madc2vcf`. If absent, results are limited depending on the format of the input.  
- **Format:** CSV containing columns: Panel_markerID, BI_markerID, Chr, Pos, Ref, Alt, Type, Indel_pos, Note

## `_v001.fa` Files

- **Purpose:** first version of the microhaplotype DB. Serves to recover possible missing Ref and Alt tags 
- **Workflow Role:** Used by all `madc2vcf`. If absent, results are limited depending on the format of the input.  
- **Format:** fasta file containing header with CloneID followed by the AlleleSequence


---

For full usage details, file format specifications, and version history, see the [BIGapp ](https://github.com/Breeding-Insight/BIGapp) and [BIGr](https://github.com/Breeding-Insight/BIGr) functions `madc2vcf_all`, `madc2vcf_targets` and `madc2vcf_multi` documentations.
