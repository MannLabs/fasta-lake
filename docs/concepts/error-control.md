# Error control and entrapment

SAGE's reported q-values address the search performed against the supplied candidates under its model. FastaLake selected those candidates using evidence from the same acquisitions. This adaptivity needs validation beyond a nominal 1% search threshold. The acceptance experiment tests implementation behavior and empirical repeatability. It does not resolve the calibration question.

Entrapment requires sequences or peptides that satisfy the chosen control's absence and exchangeability assumptions. A random half of the real reference is not automatically absent from the sample. FastaLake therefore rejects the legacy random-split mode. It also rejects an empty entrapment input, because zero opportunities to observe an entrapment match cannot certify a method.

Construction-stage matches to a shuffled control are reported as construction diagnostics. They are not labelled as search FDR certification. Whole-protein shuffling with cleavage residues fixed does not automatically create the one-to-one, unique peptide pairing required by a paired entrapment design. Interpret a scalar point estimate under its stated ratio and counting level, and consider sampling uncertainty separately. [Wen et al., Nature Methods (2025)](https://doi.org/10.1038/s41592-025-02719-x) provides the relevant methodological treatment.

!!! warning

    Do not apply Benjamini-Hochberg to an existing vector of q-values and call the result global FDR control. FastaLake rejects that legacy option. Combining independently searched acquisitions or regrouping proteins after searching requires an explicit statistical procedure if a study-wide error claim is intended.
