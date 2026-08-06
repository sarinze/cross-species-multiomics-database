# ================================================================
# STEP 28B: DIAGNOSE ORIGINAL CONSERVED/DIVERGENT RESPONSE RULES
# ================================================================

import os
import itertools
import numpy as np
import pandas as pd


print("\n" + "=" * 78)
print("STEP 28B: DIAGNOSE ORIGINAL CONSERVED/DIVERGENT RESPONSE RULES")
print("=" * 78 + "\n")


# ----------------------------------------------------------------
# 1. SETTINGS
# ----------------------------------------------------------------

INPUT_FILE = "Final_Clustered_Matrix.csv"

OUTPUT_DIRECTORY = "Conserved_Response_Rule_Diagnostics"

SUMMARY_OUTPUT = os.path.join(
    OUTPUT_DIRECTORY,
    "Candidate_Rule_Results.csv"
)

MATCHING_OUTPUT = os.path.join(
    OUTPUT_DIRECTORY,
    "Rules_Matching_7_Conserved_15_Divergent.csv"
)

EXPECTED_CONSERVED = 7
EXPECTED_DIVERGENT = 15

THRESHOLDS = [
    0.5,
    1.0,
    1.5,
    2.0
]


# ----------------------------------------------------------------
# 2. HELPER FUNCTIONS
# ----------------------------------------------------------------

def find_column(df, possible_names, required=False):

    for column in possible_names:
        if column in df.columns:
            return column

    if required:
        raise KeyError(
            f"None of these columns were found: {possible_names}"
        )

    return None


def to_numeric(df, column):

    if column is None:
        return None

    return pd.to_numeric(
        df[column],
        errors="coerce"
    )


def response_direction(series, threshold):

    direction = pd.Series(
        0,
        index=series.index,
        dtype=int
    )

    direction.loc[series >= threshold] = 1
    direction.loc[series <= -threshold] = -1

    return direction


def classify_direct_pair(
    response_a,
    response_b,
    threshold,
    require_both_responsive=True
):
    """
    Compare two individual response columns.

    Conserved:
        same non-zero response direction.

    Divergent:
        opposite non-zero response direction.
    """

    direction_a = response_direction(
        response_a,
        threshold
    )

    direction_b = response_direction(
        response_b,
        threshold
    )

    comparable = (
        response_a.notna()
        &
        response_b.notna()
    )

    if require_both_responsive:
        responsive = (
            direction_a.ne(0)
            &
            direction_b.ne(0)
        )
    else:
        responsive = (
            direction_a.ne(0)
            |
            direction_b.ne(0)
        )

    conserved = (
        comparable
        &
        responsive
        &
        direction_a.ne(0)
        &
        direction_b.ne(0)
        &
        direction_a.eq(direction_b)
    )

    divergent = (
        comparable
        &
        responsive
        &
        direction_a.ne(0)
        &
        direction_b.ne(0)
        &
        direction_a.eq(-direction_b)
    )

    return conserved, divergent


def classify_consensus(
    pseudomonas_columns,
    klebsiella_column,
    threshold,
    consensus_method="mean",
    require_internal_agreement=False,
    min_pseudomonas_measurements=1
):
    """
    Combine PA14 responses before comparing with Klebsiella.
    """

    p_matrix = pd.concat(
        pseudomonas_columns,
        axis=1
    )

    measurement_count = p_matrix.notna().sum(axis=1)

    if consensus_method == "mean":
        p_consensus = p_matrix.mean(axis=1)

    elif consensus_method == "median":
        p_consensus = p_matrix.median(axis=1)

    elif consensus_method == "maximum_absolute":
        maximum_index = (
            p_matrix.abs()
            .idxmax(axis=1)
        )

        p_consensus = pd.Series(
            np.nan,
            index=p_matrix.index
        )

        for index in p_matrix.index:
            selected_column = maximum_index.loc[index]

            if pd.notna(selected_column):
                p_consensus.loc[index] = (
                    p_matrix.loc[index, selected_column]
                )

    else:
        raise ValueError(
            f"Unknown consensus method: {consensus_method}"
        )

    p_consensus.loc[
        measurement_count < min_pseudomonas_measurements
    ] = np.nan

    if require_internal_agreement:

        signs = np.sign(p_matrix)

        positive_present = signs.eq(1).any(axis=1)
        negative_present = signs.eq(-1).any(axis=1)

        disagreement = (
            positive_present
            &
            negative_present
        )

        p_consensus.loc[disagreement] = np.nan

    return classify_direct_pair(
        p_consensus,
        klebsiella_column,
        threshold,
        require_both_responsive=True
    )


def classify_any_pair(
    comparison_pairs,
    threshold,
    final_rule="any_conserved_any_divergent"
):
    """
    Evaluate multiple comparison pairs for each gene.

    Possible final rules:

    any_conserved_any_divergent:
        divergent if any pair is opposite;
        otherwise conserved if any pair agrees.

    majority:
        classify according to whether agreements or disagreements dominate.

    unanimous:
        conserved only if all informative pairs agree;
        divergent only if all informative pairs disagree.

    score:
        sum +1 for agreement and -1 for disagreement.
    """

    conserved_pair_results = []
    divergent_pair_results = []

    for series_a, series_b in comparison_pairs:

        conserved, divergent = classify_direct_pair(
            series_a,
            series_b,
            threshold,
            require_both_responsive=True
        )

        conserved_pair_results.append(
            conserved.astype(int)
        )

        divergent_pair_results.append(
            divergent.astype(int)
        )

    conserved_count = pd.concat(
        conserved_pair_results,
        axis=1
    ).sum(axis=1)

    divergent_count = pd.concat(
        divergent_pair_results,
        axis=1
    ).sum(axis=1)

    informative_count = (
        conserved_count
        +
        divergent_count
    )

    if final_rule == "any_conserved_any_divergent":

        divergent = divergent_count.ge(1)

        conserved = (
            divergent_count.eq(0)
            &
            conserved_count.ge(1)
        )

    elif final_rule == "majority":

        conserved = (
            informative_count.ge(1)
            &
            conserved_count.gt(divergent_count)
        )

        divergent = (
            informative_count.ge(1)
            &
            divergent_count.gt(conserved_count)
        )

    elif final_rule == "unanimous":

        conserved = (
            informative_count.ge(1)
            &
            conserved_count.eq(informative_count)
        )

        divergent = (
            informative_count.ge(1)
            &
            divergent_count.eq(informative_count)
        )

    elif final_rule == "score":

        score = (
            conserved_count
            -
            divergent_count
        )

        conserved = score.gt(0)
        divergent = score.lt(0)

    else:
        raise ValueError(
            f"Unknown final rule: {final_rule}"
        )

    return conserved, divergent


# ----------------------------------------------------------------
# 3. LOAD DATA
# ----------------------------------------------------------------

if not os.path.exists(INPUT_FILE):
    raise FileNotFoundError(
        f"Input file not found: {INPUT_FILE}"
    )

os.makedirs(
    OUTPUT_DIRECTORY,
    exist_ok=True
)

df = pd.read_csv(INPUT_FILE)

print(f"Loaded: {INPUT_FILE}")
print(f"Rows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")


# ----------------------------------------------------------------
# 4. IDENTIFY RESPONSE COLUMNS
# ----------------------------------------------------------------

column_map = {
    "PA14_Colistin_RNAseq": find_column(
        df,
        [
            "Colistin_RNAseq_Fold_Change",
            "Colistin_RNAseq_FC",
            "Colistin_RNAseq_Log2FC"
        ]
    ),

    "PA14_Colistin_Riboseq": find_column(
        df,
        [
            "Colistin_Riboseq_Fold_Change",
            "Colistin_Riboseq_FC",
            "Colistin_Riboseq_Log2FC"
        ]
    ),

    "PA14_Tobramycin_RNAseq": find_column(
        df,
        [
            "Tobramycin_RNAseq_Fold_Change",
            "Tobramycin_RNAseq_FC",
            "Tobramycin_RNAseq_Log2FC"
        ]
    ),

    "PA14_Tobramycin_Riboseq": find_column(
        df,
        [
            "Tobramycin_Riboseq_Fold_Change",
            "Tobramycin_Riboseq_FC",
            "Tobramycin_Riboseq_Log2FC"
        ]
    ),

    "K56_Colistin": find_column(
        df,
        [
            "K56_vs_Colistin_Log2FC"
        ]
    ),

    "K56_Combination": find_column(
        df,
        [
            "K56_vs_Combination_Log2FC"
        ]
    ),

    "Colistin_vs_Combination": find_column(
        df,
        [
            "Colistin_vs_Combination_Log2FC"
        ]
    ),

    "KPPR1_TnSeq_FC": find_column(
        df,
        [
            "KPPR1_TnSeq_FC"
        ]
    ),

    "KPPR1_TnSeq_Log2FC": find_column(
        df,
        [
            "KPPR1_TnSeq_Log2FC"
        ]
    ),

    "PAO1_Tobramycin_TnSeq": find_column(
        df,
        [
            "PAO1_Tobramycin_Selection_Ratio"
        ]
    ),

    "PAO1_Persister_Mean": find_column(
        df,
        [
            "PAO1_Persister_All_Mean_SI",
            "Persister_All_Mean_SI"
        ]
    )
}

column_map = {
    name: column
    for name, column in column_map.items()
    if column is not None
}

print("\nDetected response columns:")

for name, column in column_map.items():
    print(f"  {name}: {column}")


numeric = {
    name: to_numeric(df, column)
    for name, column in column_map.items()
}


# ----------------------------------------------------------------
# 5. DEFINE BIOLOGICALLY PLAUSIBLE PAIRS
# ----------------------------------------------------------------

candidate_pairs = {}

if (
    "PA14_Colistin_RNAseq" in numeric
    and
    "K56_Colistin" in numeric
):
    candidate_pairs[
        "PA14 colistin RNAseq vs K56 colistin"
    ] = [
        (
            numeric["PA14_Colistin_RNAseq"],
            numeric["K56_Colistin"]
        )
    ]

if (
    "PA14_Colistin_Riboseq" in numeric
    and
    "K56_Colistin" in numeric
):
    candidate_pairs[
        "PA14 colistin Riboseq vs K56 colistin"
    ] = [
        (
            numeric["PA14_Colistin_Riboseq"],
            numeric["K56_Colistin"]
        )
    ]

if all(
    name in numeric
    for name in [
        "PA14_Colistin_RNAseq",
        "PA14_Colistin_Riboseq",
        "K56_Colistin"
    ]
):
    candidate_pairs[
        "Both PA14 colistin layers vs K56 colistin"
    ] = [
        (
            numeric["PA14_Colistin_RNAseq"],
            numeric["K56_Colistin"]
        ),
        (
            numeric["PA14_Colistin_Riboseq"],
            numeric["K56_Colistin"]
        )
    ]

if all(
    name in numeric
    for name in [
        "PA14_Colistin_RNAseq",
        "PA14_Colistin_Riboseq",
        "K56_Colistin",
        "K56_Combination"
    ]
):
    candidate_pairs[
        "PA14 colistin layers vs K56 colistin and combination"
    ] = list(
        itertools.product(
            [
                numeric["PA14_Colistin_RNAseq"],
                numeric["PA14_Colistin_Riboseq"]
            ],
            [
                numeric["K56_Colistin"],
                numeric["K56_Combination"]
            ]
        )
    )

if all(
    name in numeric
    for name in [
        "PA14_Colistin_RNAseq",
        "PA14_Colistin_Riboseq",
        "K56_Colistin",
        "Colistin_vs_Combination"
    ]
):
    candidate_pairs[
        "PA14 colistin layers vs all Klebsiella colistin comparisons"
    ] = list(
        itertools.product(
            [
                numeric["PA14_Colistin_RNAseq"],
                numeric["PA14_Colistin_Riboseq"]
            ],
            [
                numeric["K56_Colistin"],
                numeric["Colistin_vs_Combination"]
            ]
        )
    )

if all(
    name in numeric
    for name in [
        "PA14_Colistin_RNAseq",
        "PA14_Colistin_Riboseq",
        "K56_Colistin",
        "K56_Combination",
        "Colistin_vs_Combination"
    ]
):
    candidate_pairs[
        "PA14 colistin layers vs three Klebsiella DEG comparisons"
    ] = list(
        itertools.product(
            [
                numeric["PA14_Colistin_RNAseq"],
                numeric["PA14_Colistin_Riboseq"]
            ],
            [
                numeric["K56_Colistin"],
                numeric["K56_Combination"],
                numeric["Colistin_vs_Combination"]
            ]
        )
    )


# ----------------------------------------------------------------
# 6. TEST DIRECT PAIR AND MULTI-PAIR RULES
# ----------------------------------------------------------------

results = []

multi_pair_rules = [
    "any_conserved_any_divergent",
    "majority",
    "unanimous",
    "score"
]

for comparison_name, pairs in candidate_pairs.items():

    for threshold in THRESHOLDS:

        for final_rule in multi_pair_rules:

            conserved, divergent = classify_any_pair(
                comparison_pairs=pairs,
                threshold=threshold,
                final_rule=final_rule
            )

            conserved_n = int(conserved.sum())
            divergent_n = int(divergent.sum())

            results.append(
                {
                    "Analysis_Type": "Individual pair comparison",
                    "Comparison_Set": comparison_name,
                    "Threshold": threshold,
                    "Consensus_Method": "",
                    "Internal_Agreement_Required": "",
                    "Minimum_PA14_Measurements": "",
                    "Final_Rule": final_rule,
                    "Conserved_Count": conserved_n,
                    "Divergent_Count": divergent_n,
                    "Total_Classified": conserved_n + divergent_n,
                    "Distance_From_Expected": (
                        abs(conserved_n - EXPECTED_CONSERVED)
                        +
                        abs(divergent_n - EXPECTED_DIVERGENT)
                    )
                }
            )


# ----------------------------------------------------------------
# 7. TEST PA14 CONSENSUS RULES
# ----------------------------------------------------------------

if all(
    name in numeric
    for name in [
        "PA14_Colistin_RNAseq",
        "PA14_Colistin_Riboseq",
        "K56_Colistin"
    ]
):

    p_columns = [
        numeric["PA14_Colistin_RNAseq"],
        numeric["PA14_Colistin_Riboseq"]
    ]

    for threshold in THRESHOLDS:

        for consensus_method in [
            "mean",
            "median",
            "maximum_absolute"
        ]:

            for agreement_required in [
                False,
                True
            ]:

                for minimum_measurements in [
                    1,
                    2
                ]:

                    conserved, divergent = classify_consensus(
                        pseudomonas_columns=p_columns,
                        klebsiella_column=numeric["K56_Colistin"],
                        threshold=threshold,
                        consensus_method=consensus_method,
                        require_internal_agreement=agreement_required,
                        min_pseudomonas_measurements=minimum_measurements
                    )

                    conserved_n = int(conserved.sum())
                    divergent_n = int(divergent.sum())

                    results.append(
                        {
                            "Analysis_Type": "PA14 consensus comparison",
                            "Comparison_Set": (
                                "PA14 colistin consensus vs K56 colistin"
                            ),
                            "Threshold": threshold,
                            "Consensus_Method": consensus_method,
                            "Internal_Agreement_Required": (
                                agreement_required
                            ),
                            "Minimum_PA14_Measurements": (
                                minimum_measurements
                            ),
                            "Final_Rule": "direct consensus comparison",
                            "Conserved_Count": conserved_n,
                            "Divergent_Count": divergent_n,
                            "Total_Classified": (
                                conserved_n + divergent_n
                            ),
                            "Distance_From_Expected": (
                                abs(
                                    conserved_n
                                    -
                                    EXPECTED_CONSERVED
                                )
                                +
                                abs(
                                    divergent_n
                                    -
                                    EXPECTED_DIVERGENT
                                )
                            )
                        }
                    )


# ----------------------------------------------------------------
# 8. SAVE AND DISPLAY RESULTS
# ----------------------------------------------------------------

results_df = pd.DataFrame(results)

results_df = results_df.sort_values(
    [
        "Distance_From_Expected",
        "Conserved_Count",
        "Divergent_Count"
    ]
)

results_df.to_csv(
    SUMMARY_OUTPUT,
    index=False
)

matching_df = results_df[
    (
        results_df["Conserved_Count"]
        ==
        EXPECTED_CONSERVED
    )
    &
    (
        results_df["Divergent_Count"]
        ==
        EXPECTED_DIVERGENT
    )
].copy()

matching_df.to_csv(
    MATCHING_OUTPUT,
    index=False
)


print("\n" + "-" * 78)
print("TEN CLOSEST CANDIDATE RULES")
print("-" * 78 + "\n")

display_columns = [
    "Analysis_Type",
    "Comparison_Set",
    "Threshold",
    "Consensus_Method",
    "Internal_Agreement_Required",
    "Minimum_PA14_Measurements",
    "Final_Rule",
    "Conserved_Count",
    "Divergent_Count",
    "Distance_From_Expected"
]

print(
    results_df[
        display_columns
    ]
    .head(10)
    .to_string(index=False)
)


print("\n" + "-" * 78)
print("EXACT MATCHES")
print("-" * 78 + "\n")

if len(matching_df) > 0:

    print(
        matching_df[
            display_columns
        ].to_string(index=False)
    )

    print(
        "\nAt least one plausible rule reproduced exactly "
        "7 conserved and 15 divergent genes."
    )

else:

    print(
        "No tested rule reproduced exactly "
        "7 conserved and 15 divergent genes."
    )

    print(
        "\nThe original analysis probably used another response column, "
        "a significance filter, or a conserved-response score that is not "
        "represented by these candidate rules."
    )


print("\nSaved:")
print(f"  {SUMMARY_OUTPUT}")
print(f"  {MATCHING_OUTPUT}")

print("\n" + "=" * 78)
print("STEP 28C COMPLETE")
print("=" * 78 + "\n")