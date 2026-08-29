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
| p37_two_lineage_fixed_coordinate | primary | SCORED | H4_x21_over_4 | 0.0344532 | -16.554960 |
| p50_third_fixed_coordinate | primary | SCORED | H4_x21_over_4 | 0.233028 | -8.916631 |
| p50_third_fixed_coordinate | primary | SCORED | zero_effect | 15.7018 | -1.343361 |
| p45_angular_normalized_root | sensitivity | SCORED | H4_x21_over_4 | 2.42667 | -3.917646 |
| p45_angular_normalized_root | sensitivity | SCORED | zero_effect | 461.275 | 225.136793 |
| issue57_norm5 | primary | PENDING_REVEAL | — | — | — |
| p49_clean_full_curve | primary | SCORED | H4_pure_N_minus_5_over_4 | 37.8868 | 4.273907 |
| p49_clean_full_curve | primary | SCORED | q2_even_scalar | 1.78993 | -13.272009 |
| p49_clean_full_curve | primary | SCORED | rank2_Jordan_log | 0.676706 | -13.758148 |
| p49_clean_full_curve | primary | SCORED | zero_effect | 1044.57 | 507.290005 |
| p91_n145_to_n290_full_curve | primary | PENDING_REVEAL | — | — | — |

## Pairwise comparisons on matched primary endpoints

Negative delta favors the left model.

| Left | Right | Blocks | Delta NLPD | Preferred |
|---|---|---|---:|---|
| H4_pure_N_minus_5_over_4 | q2_even_scalar | p49_clean_full_curve | 17.545916 | q2_even_scalar |
| H4_pure_N_minus_5_over_4 | rank2_Jordan_log | p49_clean_full_curve | 18.032055 | rank2_Jordan_log |
| H4_pure_N_minus_5_over_4 | zero_effect | p49_clean_full_curve | -503.016098 | H4_pure_N_minus_5_over_4 |
| H4_x17_over_4 | H4_x21_over_4 | issue43_n185_n265_deltaM | 13.673035 | H4_x21_over_4 |
| H4_x17_over_4 | zero_effect | issue43_n185_n265_deltaM | 0.540595 | zero_effect |
| H4_x21_over_4 | zero_effect | issue43_n185_n265_deltaM, p50_third_fixed_coordinate | -20.705711 | H4_x21_over_4 |
| q2_even_scalar | rank2_Jordan_log | p49_clean_full_curve | 0.486139 | rank2_Jordan_log |
| q2_even_scalar | zero_effect | p49_clean_full_curve | -520.562014 | q2_even_scalar |
| rank2_Jordan_log | zero_effect | p49_clean_full_curve | -521.048153 | rank2_Jordan_log |
