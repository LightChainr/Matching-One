# Prequential evidence ledger

Only `primary` rows enter cumulative evidence.

| Block | Role | Status | Model | chi2 | NLPD |
|---|---|---|---|---:|---:|
| issue43_n185_n265_deltaM | primary | SCORED | H4_x17_over_4 | 30.2461 | -3.468488 |
| issue43_n185_n265_deltaM | primary | SCORED | H4_x21_over_4 | 3.04598 | -17.141523 |
| issue43_n185_n265_deltaM | primary | SCORED | zero_effect | 29.4094 | -4.009083 |
| issue43_n185_n265_sprime_training_only | sensitivity | SCORED | analytic_inverse_N | 0.862214 | -0.891629 |
| issue43_n185_n265_sprime_training_only | sensitivity | SCORED | rank2_Jordan_log | 1.2036 | -0.376709 |
| issue43_deltaS_literal_channel_mismatch | protocol_history | PROTOCOL_FAILURE_CHANNEL_MISMATCH | — | — | — |
| p37_two_lineage_fixed_coordinate | primary | PENDING_ARTIFACT_EXTRACTION | — | — | — |
| p50_third_fixed_coordinate | primary | PENDING_ARTIFACT_EXTRACTION | — | — | — |
| p45_angular_normalized_root | primary | PENDING_ARTIFACT_EXTRACTION | — | — | — |
| issue57_norm5 | primary | PENDING_REVEAL | — | — | — |

## Pairwise comparisons on matched primary endpoints

Negative delta favors the left model.

| Left | Right | Blocks | Delta NLPD | Preferred |
|---|---|---|---:|---|
| H4_x17_over_4 | H4_x21_over_4 | issue43_n185_n265_deltaM | 13.673035 | H4_x21_over_4 |
| H4_x17_over_4 | zero_effect | issue43_n185_n265_deltaM | 0.540595 | zero_effect |
| H4_x21_over_4 | zero_effect | issue43_n185_n265_deltaM | -13.132441 | H4_x21_over_4 |
