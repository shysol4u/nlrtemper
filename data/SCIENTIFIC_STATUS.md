# Scientific status — rules v0.3.1

This is an exploratory evidence viewer, not a validated temperature-phenotype predictor.

NLEEL-like hits retain their one-based coordinates, identity, conservative similarity and fourth-residue (second-E) status. Hits fully inside annotated LRR 1 take priority, ranked by earliest start (v0.4.0). Otherwise the earliest qualifying hit in the annotated LRR region is preferred, with the ±50-residue proximal window as an unresolved fallback; the most conserved and nearest hits remain separate. The window is an exploratory choice, not a validated biological cutoff. No annotated boundary means no verified first-LRR localization.

First-LRR proximity is a preference, not a universal property: in the RefPlantNLR computational annotations, the user-supplied SNC1 reference has NLEEL637–641 in a later projected repeat, whereas the first projected LRR begins at 543. Sequence conservation does not override conflicting biological evidence. Motif absence is not evidence of temperature tolerance.

The user-supplied `axxaxxa/LxxL/IxaxxCxaxxaxx` shorthand is mapped as an exploratory hypothesis, with `a=AILMV` and `x=any standard amino acid`. Its exact grammar, substitutions and ordering have not been independently verified as the published Grech-Baran classifier. The reported interdomain stabilization mechanism motivates inspection, but does not establish these regex rules or a phenotype prediction.

WHD–LRR contact counts require a uniquely matched structure chain and explicit WHD/LRR annotations. Counts are unique residue pairs with any heavy atoms within 4.5 Å (sequence separation at least 5). They are geometric descriptors, not hydrogen-bond counts, binding energies or proof of thermal stabilization. Missing atoms, conformational state and uncertain domain placement affect interpretation. Generic crystallographic B factors are not pLDDT; PAE is not assessed.

Exact-sequence curated phenotypes can inform the integrated display, while the sequence-only call remains visible. Phenotype-informed calls must never be scored against their own input labels. Unknown phenotype is not a negative control. Counts alone do not trigger TT-like calls.

Sources: [RefPlantNLR](https://doi.org/10.1371/journal.pbio.3001124), [Zhu et al.](https://doi.org/10.1371/journal.ppat.1000844), [Ry_sto](https://doi.org/10.1111/pbi.13230), [Sr6](https://doi.org/10.1038/s41467-025-60030-x), [Grech-Baran preprint](https://doi.org/10.64898/2025.12.17.694812).

IP review remains incomplete. [WO2025019283A1](https://patents.google.com/patent/WO2025019283A1/en) is relevant to engineered immune-receptor thermostability. This code provides observational analysis; publication and open-source licensing do not establish freedom to operate.
