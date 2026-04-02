# Simulated files for tests

## Alfalfa

* alfalfa_IUPAC.csv: 300 CloneIDs and 30 samples simulated - AlleleSequences have IUPAC codes
* alfalfa_lowercase.csv: 300 CloneIDs and 30 samples simulated - AlleleSequences have mix of lower and upper case and it has 3 Ref and 1 Alt tag missing
* alfalfa_madc_raw.csv: 300 CloneIDs and 30 samples simulated - It has the raw MADC header and AlleleID is not fixed (doesn't have _0001 or _0002)
* alfalfa_madc_wrongID.csv: 300 CloneIDs and 30 samples simulated - CloneID doesn't match botloci (CloneID = chr01_[position]; botloci = chr01.1_[position])
* alfalfa_madc.csv: 300 CloneIDs and 30 samples simulated - all clean
* alfalfa_marker_info_ChromPos.csv: version of _lut.csv file without the REF and ALT

## Potato

* potato_indel_IUPAC.csv: 300 CloneIDs and 30 samples simulated - contains indels and AlleleSequences have IUPAC codes
* potato_indel_lowercase.csv: 300 CloneIDs and 30 samples simulated - contains indels and AlleleSequences have mix of lower and upper case 
* potato_indel_madc.csv: 300 CloneIDs and 30 samples simulated - all clean
* potato_marker_info_chrompos.csv: version of _lut.csv file without the REF and ALT
* potato_more_indels_madc_ChromPosFALSE.csv: 300 CloneIDs and 30 samples simulated - contains indels and not all CloneID have the Chr_Pos format

* generate_simu_madc.R: script used to generate the simulated MADC files