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

## MADC test files:

### Alfalfa

* `alfalfa_IUPAC.csv`: 
    * 300 CloneIDs and 30 samples simulated
    * AlleleSequences have IUPAC codes
* `alfalfa_lowercase.csv`: 
    * 300 CloneIDs and 30 samples simulated
    * AlleleSequences have mix of lower and upper case
    * has 3 Ref and 1 Alt tag missing
* `alfalfa_madc_raw.csv`: 
    * 300 CloneIDs and 30 samples simulated
    * has the raw MADC header
    * AlleleID is not fixed (doesn't have _0001 or _0002)
* `alfalfa_madc_wrongID.csv`: 
    * 300 CloneIDs and 30 samples simulated
    * CloneID doesn't match botloci (CloneID = chr01_[position]; botloci = chr01.1_[position])
* `alfalfa_madc.csv`: 
    * 300 CloneIDs and 30 samples simulated
    * all clean
* `alfalfa_marker_info_ChromPos.csv`: 
    * version of _lut.csv file without the REF and ALT

### Potato

* `potato_indel_IUPAC.csv`: 
    * 300 CloneIDs and 30 samples simulated
    * contains indels and AlleleSequences have IUPAC codes
* `potato_indel_lowercase.csv`: 
    * 300 CloneIDs and 30 samples simulated
    * contains indels 
    * AlleleSequences have mix of lower and upper case 
* `potato_indel_madc.csv`: 
    * 300 CloneIDs and 30 samples simulated 
    * contains indels
* `potato_marker_info_chrompos.csv`: 
    * version of _lut.csv file without the REF and ALT
* `potato_more_indels_madc_ChromPosFALSE.csv`: 
    * 300 CloneIDs and 30 samples simulated 
    * contains indels 
    * not all CloneID have the Chr_Pos format

### Simulation script
* `generate_simu_madc.R`: script used to generate the simulated MADC files

---

For full usage details, file format specifications, and version history, see the [BIGapp ](https://github.com/Breeding-Insight/BIGapp) and [BIGr](https://github.com/Breeding-Insight/BIGr) functions `madc2vcf_all`, `madc2vcf_targets` and `madc2vcf_multi` documentations.
