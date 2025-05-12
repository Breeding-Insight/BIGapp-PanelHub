# BIGapp-PanelHub

A central metadata repository of marker-panel definitions for downstream workflows in [BIGapp](https://github.com/Breeding-Insight/BIGapp) and [BIGr](https://github.com/Breeding-Insight/BIGr).

## `.botloci` Files

- **Purpose:** Enumerate all marker IDs for which bottom-strand probes were designed.  
- **Workflow Role:** During MADC-to-VCF conversion, these files tell BIGapp/BIGr exactly which sequences to extract as reverse complements.  
- **Format:** Plain-text lists—one marker ID per line (optionally namespaced by panel or species).

---

For full usage details, file format specifications, and version history, see the [BIGapp ](https://github.com/Breeding-Insight/BIGapp) and [BIGr](https://github.com/Breeding-Insight/BIGr) functions `madc2vcf_all` and `madc2vcf_targets` documentations.
