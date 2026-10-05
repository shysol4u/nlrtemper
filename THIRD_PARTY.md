# Third-party materials and review status

The MIT license applies to original project code, not third-party code, sequences, annotations, models or publication material. Original website text and graphics are additionally licensed by Shyam Solanki under CC BY 4.0; see CONTENT-LICENSE.md. This does not relicense third-party material.

- 3Dmol.js 2.5.5 is in `dist/vendor/`, with its original license and bundled-component notices in `3Dmol-LICENSE`.
- RefPlantNLR sequences/annotations derive from [Kourelis et al., 2021](https://doi.org/10.1371/journal.pbio.3001124). Preserve attribution and applicable source-data terms. The article is CC BY; sequence origins are recorded in the phenotype ledger.
- MorexV3 and Arabidopsis candidates derive from Ensembl Plants release 62. URLs and checksums are in `validation/sources.json`.
- Sr6 derives from PP949235.1 / XCA47734.1; [Hewitt et al., 2025](https://doi.org/10.1038/s41467-025-60030-x).
- Optional AlphaFold DB models are retrieved from their original service, not redistributed here. Retain AlphaFold/UniProt attribution and source terms when exporting models.
- Reconstruction requires Biopython, pyhmmer and openpyxl (versions in `validation/requirements.txt`); they are not browser dependencies and are not vendored.
- NLRtracker is not copied or integrated; its inspected MIT license does not cover every dependency/database in a future integration.

Dependency and data licensing review is preliminary; public release requires completion. Relevant patent-family review remains pending, including WO2025019283A1. No freedom-to-operate conclusion is made.
