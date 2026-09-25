"""Keep null denominators paired while measuring increased matching opportunities."""

import gzip
import json

from fasta_lake.bags import select_databases
from tools.run_bag_nulls import NULL_ARMS, retrieval


def test_null_bag_inflation_is_visible_without_extra_evidence_votes(tmp_path):
    evidence = tmp_path / "evidence.jsonl.gz"
    candidates = tmp_path / "candidates.jsonl.gz"
    with gzip.open(evidence, "wt") as handle:
        for q, control in ((0, 0), (1, 1), (2, 2), (3, 1)):
            handle.write(
                json.dumps(dict(query=q, control=control, lookup=True, sequence="PEPTIDERK")) + "\n"
            )
    records = [
        dict(header="bag_only_1", sequence="PETPIDERK", lookup_routes=16, evidence=[[1, 16]]),
        dict(header="bag_only_2", sequence="PETPIDERK", lookup_routes=16, evidence=[[1, 16]]),
    ]
    with gzip.open(candidates, "wt") as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")
    result = retrieval(candidates, evidence, tmp_path)
    paired = {r["arm"]: r for r in result["paired_comparisons"]}
    assert paired["shuffle_1_exact"]["queried"] == paired["shuffle_1_bags"]["queried"] == 2
    assert paired["shuffle_1_exact"]["matched"] == 0
    assert paired["shuffle_1_bags"]["matched"] == 1
    assert paired["shuffle_1_bags"]["gained"] == 1
    assert paired["shuffle_1_bags"]["ratio_to_exact"] is None
    assert paired["shuffle_1_bags"]["candidate_proteins"] == 2
    selected = select_databases(candidates, tmp_path, 4, arms=NULL_ARMS)
    assert selected["arms"]["shuffle_1_exact"]["selected_proteins"] == 0
    assert selected["arms"]["shuffle_1_bags"]["assigned_original_queries"] == 1
    assert selected["arms"]["shuffle_1_bags"]["selected_proteins"] == 1
