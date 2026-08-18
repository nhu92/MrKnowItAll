from __future__ import annotations

from pathlib import Path

from .backbone import Backbone
from .fasta import FastaRecord, ungap
from .qc import best_stop_fraction


def train_anomaly_model(backbone_dir: Path, output_path: Path, seed: int = 353) -> dict:
    """Train an optional unsupervised anomaly detector on curated backbone sequences."""
    try:
        import joblib
        import numpy as np
        from sklearn.ensemble import IsolationForest
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
    except ImportError as exc:
        raise RuntimeError("Install the optional ML dependencies: pip install .[ml]") from exc

    backbone = Backbone(backbone_dir)
    rows: list[list[float]] = []
    for locus in backbone.loci.values():
        median_length = locus.median_ungapped_length
        for record in locus.records:
            clean = ungap(record.sequence)
            # Kew target-capture recoveries are frequently legitimate partial CDS. Augment the
            # curated full-length distribution with in-frame fragments so the model learns that
            # truncation within the deterministic hard gates is not itself anomalous.
            for fraction in (0.4, 0.6, 0.8, 1.0):
                width = max(3, int(len(clean) * fraction) // 3 * 3)
                start = max(0, (len(clean) - width) // 2)
                fragment = clean[start : start + width]
                ambiguity = sum(base not in "ACGT" for base in fragment) / max(1, len(fragment))
                rows.append(
                    [
                        len(fragment) / max(1, median_length),
                        ambiguity,
                        best_stop_fraction(fragment),
                    ]
                )
    model = make_pipeline(
        StandardScaler(),
        IsolationForest(n_estimators=300, contamination=0.01, random_state=seed, n_jobs=-1),
    )
    model.fit(np.asarray(rows, dtype=float))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, output_path)
    return {"training_sequences": len(rows), "loci": len(backbone.loci), "output": str(output_path)}


def anomaly_score(model_path: Path, record: FastaRecord, median_length: float) -> float:
    import joblib
    import numpy as np

    model = joblib.load(model_path)
    clean = ungap(record.sequence)
    ambiguity = sum(base not in "ACGT" for base in clean) / max(1, len(clean))
    row = np.asarray(
        [[len(clean) / max(1, median_length), ambiguity, best_stop_fraction(clean)]], dtype=float
    )
    return float(model.decision_function(row)[0])
