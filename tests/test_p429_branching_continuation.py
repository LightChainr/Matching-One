import csv
import sys


sys.path.insert(0, "scripts")
from score_p429_branching_continuation import environment_estimate, load_rows  # noqa: E402


def test_branching_estimands_are_separate(tmp_path):
    path = tmp_path / "branch.csv"
    fields = ["n", "orientation", "batch", "replica",
              "checkpoint_b1_safe_count", "branch_common_safe",
              "branch_clone1_survives", "branch_clone2_survives",
              "branch_both_survive"]
    values = [(1, 1, 1), (1, 0, 0), (0, 1, 0), (0, 0, 0)]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for replica, (y1, y2, both) in enumerate(values):
            writer.writerow({"n": 325, "orientation": "first", "batch": 0,
                             "replica": replica, "checkpoint_b1_safe_count": 7,
                             "branch_common_safe": 1, "branch_clone1_survives": y1,
                             "branch_clone2_survives": y2,
                             "branch_both_survive": both})
    estimate = environment_estimate(load_rows(path))
    assert estimate["b2_survival_estimate"] == 0.5
    assert estimate["branch_success"] == 0.25
    assert estimate["clone_dependence_gap"] == 0.0
