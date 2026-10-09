#!/usr/bin/env python3
"""Fetch pinned public UCI artifacts and write declarations for the C++ importer.

No mathematical analysis runs here. Original processed files, codes and dictionary
bytes are retained. Licensing/attribution: https://doi.org/10.24432/C52P4X,
CC BY 4.0; Janosi, Steinbrunn, Pfisterer and Detrano, UCI repository.
"""
import argparse
import hashlib
import json
import re
import urllib.request
from pathlib import Path

BASE = "https://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/"
PINNED = {
    "processed.cleveland.data": "a74b7efa387bc9d108d7d0115d831fe9b414b29ae7124f331b622b4efa0427c8",
    "processed.hungarian.data": "d1ad108f785768cd3d7e82dc522e6f5a61eea93cccfb3a46ee8076f73fc3d796",
    "heart-disease.names": "ec7ed5aa4ed8321097808f47b836301c2a0c36f3bb09c40b8b43b622363e6b61",
}
# Name, channel, quantitative interpretation, units, type. Unspecified source
# units remain unspecified; codes never become continuous measurements.
COLUMNS = [
    ("age", "covariate", "reported_age", "year", "quantitative"),
    ("sex", "covariate", "reported_sex_code", "source_code", "categorical"),
    ("cp", "covariate", "reported_chest_pain_code", "source_code", "categorical"),
    ("trestbps", "observation", "reported_resting_blood_pressure", "mm Hg", "quantitative"),
    ("chol", "observation", "reported_serum_cholesterol", "mg/dl", "quantitative"),
    ("fbs", "covariate", "reported_fasting_glucose_threshold_code", "source_code", "categorical"),
    ("restecg", "observation", "reported_resting_ECG_code", "source_code", "categorical"),
    ("thalach", "observation", "reported_maximum_heart_rate", "source_unit_unspecified", "quantitative"),
    ("exang", "covariate", "reported_exercise_angina_code", "source_code", "categorical"),
    ("oldpeak", "observation", "reported_ST_depression", "source_unit_unspecified", "quantitative"),
    ("slope", "observation", "reported_ST_slope_code", "source_code", "categorical"),
    ("ca", "observation", "reported_major_vessel_count", "vessel_count", "quantitative"),
    ("thal", "observation", "reported_thal_code", "source_code", "categorical"),
    ("num", "observation", "reported_diagnostic_code", "source_code", "categorical"),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("data/public/uci-heart"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for name, expected in PINNED.items():
        path = args.out / name
        if not path.exists():
            with urllib.request.urlopen(BASE + name, timeout=60) as response:
                data = response.read(2 * 1024 * 1024 + 1)
            if len(data) > 2 * 1024 * 1024 or hashlib.sha256(data).hexdigest() != expected:
                raise ValueError(f"unexpected artifact bytes: {name}")
            path.write_bytes(data)
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"local artifact changed: {name}")
    dictionary = (args.out / "heart-disease.names").read_bytes()
    columns = []
    for name, channel, quantity, units, value_type in COLUMNS:
        match = re.search(rb"(?m)^\s+\d+ " + name.encode() + rb"(?::| =)", dictionary)
        if match is None:
            raise ValueError(f"dictionary entry absent: {name}")
        next_entry = re.search(rb"(?m)^\s+\d+ [a-z]", dictionary[match.end():])
        end = match.end() + next_entry.start() if next_entry else len(dictionary)
        columns.append(dict(id=name, channel=channel, quantity=quantity, units=units,
                            value_type=value_type, dictionary_start=match.start(), dictionary_end=end))
    for site, namespace in [("cleveland", "45001"), ("hungarian", "45002")]:
        filename = f"processed.{site}.data"
        spec = dict(schema="gadget.table.import.v1", id=f"uci-heart-{site}", source_namespace=namespace,
                    source=dict(path=filename, sha256="sha256:" + PINNED[filename], delimiter=",", header=False),
                    dictionary=dict(path="heart-disease.names", sha256="sha256:" + PINNED["heart-disease.names"]),
                    columns=columns, missing_tokens=["?"], relation_columns=["sex", "cp"],
                    question="Understand reported clinical variables, support and declared categorical membership relations.",
                    population=f"retained UCI processed {site} records; row identity is not an inferred unique-person identifier",
                    selection=f"all original processed {site} rows and 14 variables; no complete-case filtering or target recoding",
                    limitations=["Observational source cohort; no population representativeness, diagnosis or causal adjustment established.",
                                 "Categorical group membership is an explicitly declared observation relation, not a biological interaction.",
                                 "No visit times, genomic/expression observations or cross-site patient correspondence are provided.",
                                 "Original codes and missing tokens are retained; inter-site code equivalence is not inferred.",
                                 "Source-stalk coordinates 900 follow relation_columns order; indicator basis and constant identity maps describe membership only."],
                    source_url="https://archive.ics.uci.edu/dataset/45/heart+disease",
                    attribution="Janosi, Steinbrunn, Pfisterer and Detrano; UCI Machine Learning Repository; DOI 10.24432/C52P4X",
                    license="CC BY 4.0")
        declaration = args.out / f"{site}-import.json"
        content = json.dumps(spec, indent=2) + "\n"
        if declaration.exists() and declaration.read_text() != content:
            raise ValueError(f"local declaration changed; preserve it and select another output: {declaration}")
        if not declaration.exists():
            declaration.write_text(content)
    print("PASS: pinned source artifacts and two declarations; execute imports with the C++ gadget CLI")


if __name__ == "__main__":
    main()
