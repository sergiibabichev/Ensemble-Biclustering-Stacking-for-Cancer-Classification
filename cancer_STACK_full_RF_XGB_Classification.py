"""
Classification module for STACK_full biclustering-derived gene subset.

Purpose
-------
This script evaluates the same feature subset using two classifier families:
1) Random Forest (bagging-based baseline used in the original manuscript)
2) XGBoost (gradient boosting alternative requested by the reviewer)

Outputs
-------
All outputs are saved in: results_STACK_full

- best_rf_params_STACK_full.json
- best_xgb_params_STACK_full.json
- classification_report_RF_STACK_full.csv
- classification_report_XGB_STACK_full.csv
- classification_report_RF_STACK_full.png
- classification_report_XGB_STACK_full.png
- cv_scores_RF_STACK_full.csv
- cv_scores_XGB_STACK_full.csv
- metrics_summary_STACK_full.csv
- metrics_summary_STACK_full.json
- confusion_matrix_RF_STACK_full.csv
- confusion_matrix_XGB_STACK_full.csv

Notes
-----
- Confidence intervals are computed as 95% CI from repeated stratified CV scores:
  mean ± 1.96 * SD / sqrt(n_folds_total).
- XGBoost labels are encoded using LabelEncoder and then decoded for reports.
"""

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import (
    train_test_split,
    StratifiedKFold,
    RepeatedStratifiedKFold,
    cross_val_score
)
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.preprocessing import LabelEncoder

from bayes_opt import BayesianOptimization

try:
    from xgboost import XGBClassifier
except ImportError as exc:
    raise ImportError(
        "The xgboost package is required for this script. "
        "Install it with: pip install xgboost"
    ) from exc


# ---------------------------------------------------------------------
# 0) Configuration
# ---------------------------------------------------------------------
DATA_CSV = "STACK_subset_unique_genes.csv"
TAG = "STACK_full"
OUT_DIR = Path("results_STACK_full")
OUT_DIR.mkdir(parents=True, exist_ok=True)

RANDOM_STATE = 42
TEST_SIZE = 0.30

# Bayesian optimization settings
RF_INIT_POINTS = 4
RF_N_ITER = 12

XGB_INIT_POINTS = 4
XGB_N_ITER = 12

# Final confidence interval settings
CV_SPLITS = 5
CV_REPEATS = 3
N_JOBS = -1


# ---------------------------------------------------------------------
# 1) Utility functions
# ---------------------------------------------------------------------
def ci95_from_scores(scores: np.ndarray) -> tuple[float, float, float]:
    """Return mean, standard deviation, and 95% CI half-width."""
    scores = np.asarray(scores, dtype=float)
    mean_score = float(np.mean(scores))
    std_score = float(np.std(scores, ddof=1)) if len(scores) > 1 else 0.0
    ci95 = float(1.96 * std_score / np.sqrt(len(scores))) if len(scores) > 1 else 0.0
    return mean_score, std_score, ci95


def save_report_table(report_df: pd.DataFrame, output_png: Path, title: str) -> None:
    """Save classification report as a PNG table."""
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.axis("tight")
    ax.axis("off")
    ax.set_title(title, fontsize=12, pad=12)

    table = ax.table(
        cellText=report_df.values,
        colLabels=report_df.columns,
        rowLabels=report_df.index,
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.2)

    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_report_dataframe(y_true, y_pred, labels=None) -> pd.DataFrame:
    """Create classification report dataframe with accuracy and weighted F1 rows."""
    report = classification_report(
        y_true,
        y_pred,
        labels=labels,
        output_dict=True,
        zero_division=0
    )
    report_df = pd.DataFrame(report).transpose()

    accuracy = accuracy_score(y_true, y_pred)
    weighted_f1 = report_df.loc["weighted avg", "f1-score"]

    report_df.loc["Accuracy"] = [accuracy, None, None, None]
    report_df.loc["Weighted F1-Score"] = [weighted_f1, None, None, None]
    return report_df.round(4)


def evaluate_final_model(
    model,
    model_name: str,
    X_train,
    y_train,
    X_test,
    y_test,
    labels_for_report=None,
    decode_predictions=None
) -> dict:
    """Run repeated CV, train final model, evaluate test set, and save outputs."""
    cv = RepeatedStratifiedKFold(
        n_splits=CV_SPLITS,
        n_repeats=CV_REPEATS,
        random_state=RANDOM_STATE
    )

    t0 = time.time()
    cv_scores = cross_val_score(
        model,
        X_train,
        y_train,
        cv=cv,
        scoring="accuracy",
        n_jobs=N_JOBS
    )
    cv_time = time.time() - t0

    cv_mean, cv_std, cv_ci95 = ci95_from_scores(cv_scores)

    pd.DataFrame({"accuracy": cv_scores}).to_csv(
        OUT_DIR / f"cv_scores_{model_name}_{TAG}.csv",
        index=False
    )

    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    if decode_predictions is not None:
        y_test_report = decode_predictions(y_test)
        y_pred_report = decode_predictions(y_pred)
        labels_report = labels_for_report
    else:
        y_test_report = y_test
        y_pred_report = y_pred
        labels_report = labels_for_report

    test_accuracy = float(accuracy_score(y_test_report, y_pred_report))
    report_df = build_report_dataframe(
        y_test_report,
        y_pred_report,
        labels=labels_report
    )
    weighted_f1 = float(report_df.loc["Weighted F1-Score", "precision"])

    report_csv = OUT_DIR / f"classification_report_{model_name}_{TAG}.csv"
    report_png = OUT_DIR / f"classification_report_{model_name}_{TAG}.png"
    report_df.to_csv(report_csv)
    save_report_table(
        report_df,
        report_png,
        title=f"{model_name} classification report: {TAG}"
    )

    cm = confusion_matrix(y_test_report, y_pred_report, labels=labels_report)
    cm_df = pd.DataFrame(cm, index=labels_report, columns=labels_report)
    cm_df.to_csv(OUT_DIR / f"confusion_matrix_{model_name}_{TAG}.csv")

    return {
        "method": TAG,
        "classifier": model_name,
        "cv_accuracy_mean": round(cv_mean, 4),
        "cv_accuracy_std": round(cv_std, 4),
        "cv_accuracy_ci95": round(cv_ci95, 4),
        "cv_accuracy_formatted": f"{cv_mean:.4f} ± {cv_ci95:.4f}",
        "test_accuracy": round(test_accuracy, 4),
        "test_weighted_f1": round(weighted_f1, 4),
        "cv_total_folds": int(CV_SPLITS * CV_REPEATS),
        "cv_time_sec": round(cv_time, 2),
        "n_train": int(len(y_train)),
        "n_test": int(len(y_test)),
        "n_features": int(X_train.shape[1])
    }


# ---------------------------------------------------------------------
# 2) Data loading
# ---------------------------------------------------------------------
data_full = pd.read_csv(DATA_CSV)

if "Class" not in data_full.columns:
    raise ValueError("Input data must contain a 'Class' column.")

X = data_full.drop(columns=["Class"])
y = data_full["Class"].astype(str)

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=TEST_SIZE,
    random_state=RANDOM_STATE,
    stratify=y
)

p = X_train.shape[1]
class_labels = sorted(y.unique().tolist())

print(f"[INFO] Method: {TAG}")
print(f"[INFO] Data: {DATA_CSV}")
print(f"[INFO] Samples: {X.shape[0]}, features: {X.shape[1]}, classes: {len(class_labels)}")
print(f"[INFO] Train: {X_train.shape}, test: {X_test.shape}")


# ---------------------------------------------------------------------
# 3) Random Forest with Bayesian optimization
# ---------------------------------------------------------------------
max_feat_upper = min(0.3, max(0.05, 256.0 / float(p)))
max_feat_lower = 0.05 if 0.05 < max_feat_upper else max_feat_upper * 0.5

rf_bounds = {
    "n_estimators": (50, 250),
    "max_depth": (3, 24),
    "min_samples_split": (2, 10),
    "min_samples_leaf": (1, 4),
    "max_features": (max_feat_lower, max_feat_upper),
    "bootstrap": (0, 1),
    "criterion": (0, 1)
}


def rf_cv(
    n_estimators,
    max_depth,
    min_samples_split,
    min_samples_leaf,
    max_features,
    bootstrap,
    criterion
):
    try:
        bootstrap_value = bool(round(bootstrap))
        criterion_value = "gini" if round(criterion) == 0 else "entropy"

        model = RandomForestClassifier(
            n_estimators=int(n_estimators),
            max_depth=None if max_depth < 1 else int(max_depth),
            min_samples_split=int(min_samples_split),
            min_samples_leaf=int(min_samples_leaf),
            max_features=float(max_features),
            bootstrap=bootstrap_value,
            criterion=criterion_value,
            n_jobs=N_JOBS,
            random_state=RANDOM_STATE
        )

        cv = StratifiedKFold(
            n_splits=3,
            shuffle=True,
            random_state=RANDOM_STATE
        )
        score = np.mean(cross_val_score(
            model,
            X_train,
            y_train,
            cv=cv,
            scoring="accuracy",
            n_jobs=N_JOBS
        ))
        return score
    except Exception as exc:
        print(f"[RF BO ERROR] {exc}")
        return 0.0


print("\n[INFO] Starting Random Forest Bayesian optimization...")
rf_optimizer = BayesianOptimization(
    f=rf_cv,
    pbounds=rf_bounds,
    random_state=RANDOM_STATE,
    verbose=2
)
rf_optimizer.maximize(init_points=RF_INIT_POINTS, n_iter=RF_N_ITER)

rf_best_raw = rf_optimizer.max["params"]
with open(OUT_DIR / f"best_rf_params_raw_{TAG}.json", "w") as f:
    json.dump(rf_best_raw, f, indent=4)

rf_best = {
    "n_estimators": int(rf_best_raw["n_estimators"]),
    "max_depth": None if rf_best_raw["max_depth"] < 1 else int(rf_best_raw["max_depth"]),
    "min_samples_split": int(rf_best_raw["min_samples_split"]),
    "min_samples_leaf": int(rf_best_raw["min_samples_leaf"]),
    "max_features": float(rf_best_raw["max_features"]),
    "bootstrap": bool(round(rf_best_raw["bootstrap"])),
    "criterion": "gini" if round(rf_best_raw["criterion"]) == 0 else "entropy",
}

with open(OUT_DIR / f"best_rf_params_{TAG}.json", "w") as f:
    json.dump(rf_best, f, indent=4)

rf_model = RandomForestClassifier(
    **rf_best,
    n_jobs=N_JOBS,
    random_state=RANDOM_STATE
)


# ---------------------------------------------------------------------
# 4) XGBoost with Bayesian optimization
# ---------------------------------------------------------------------
label_encoder = LabelEncoder()
y_train_enc = label_encoder.fit_transform(y_train)
y_test_enc = label_encoder.transform(y_test)

xgb_bounds = {
    "n_estimators": (50, 300),
    "max_depth": (2, 8),
    "learning_rate": (0.01, 0.25),
    "subsample": (0.6, 1.0),
    "colsample_bytree": (0.3, 1.0),
    "min_child_weight": (1, 10),
    "gamma": (0.0, 5.0),
    "reg_alpha": (0.0, 2.0),
    "reg_lambda": (0.5, 5.0)
}


def make_xgb_model(params: dict) -> XGBClassifier:
    """Create XGBClassifier with version-tolerant parameters."""
    return XGBClassifier(
        objective="multi:softprob",
        num_class=len(label_encoder.classes_),
        eval_metric="mlogloss",
        tree_method="hist",
        n_estimators=int(params["n_estimators"]),
        max_depth=int(params["max_depth"]),
        learning_rate=float(params["learning_rate"]),
        subsample=float(params["subsample"]),
        colsample_bytree=float(params["colsample_bytree"]),
        min_child_weight=float(params["min_child_weight"]),
        gamma=float(params["gamma"]),
        reg_alpha=float(params["reg_alpha"]),
        reg_lambda=float(params["reg_lambda"]),
        n_jobs=N_JOBS,
        random_state=RANDOM_STATE
    )


def xgb_cv(
    n_estimators,
    max_depth,
    learning_rate,
    subsample,
    colsample_bytree,
    min_child_weight,
    gamma,
    reg_alpha,
    reg_lambda
):
    try:
        params = {
            "n_estimators": n_estimators,
            "max_depth": max_depth,
            "learning_rate": learning_rate,
            "subsample": subsample,
            "colsample_bytree": colsample_bytree,
            "min_child_weight": min_child_weight,
            "gamma": gamma,
            "reg_alpha": reg_alpha,
            "reg_lambda": reg_lambda,
        }

        model = make_xgb_model(params)

        cv = StratifiedKFold(
            n_splits=3,
            shuffle=True,
            random_state=RANDOM_STATE
        )

        # n_jobs=1 avoids nested parallelism conflicts with XGBoost.
        score = np.mean(cross_val_score(
            model,
            X_train,
            y_train_enc,
            cv=cv,
            scoring="accuracy",
            n_jobs=1
        ))
        return score
    except Exception as exc:
        print(f"[XGB BO ERROR] {exc}")
        return 0.0


print("\n[INFO] Starting XGBoost Bayesian optimization...")
xgb_optimizer = BayesianOptimization(
    f=xgb_cv,
    pbounds=xgb_bounds,
    random_state=RANDOM_STATE,
    verbose=2
)
xgb_optimizer.maximize(init_points=XGB_INIT_POINTS, n_iter=XGB_N_ITER)

xgb_best_raw = xgb_optimizer.max["params"]
with open(OUT_DIR / f"best_xgb_params_raw_{TAG}.json", "w") as f:
    json.dump(xgb_best_raw, f, indent=4)

xgb_best = {
    "n_estimators": int(xgb_best_raw["n_estimators"]),
    "max_depth": int(xgb_best_raw["max_depth"]),
    "learning_rate": float(xgb_best_raw["learning_rate"]),
    "subsample": float(xgb_best_raw["subsample"]),
    "colsample_bytree": float(xgb_best_raw["colsample_bytree"]),
    "min_child_weight": float(xgb_best_raw["min_child_weight"]),
    "gamma": float(xgb_best_raw["gamma"]),
    "reg_alpha": float(xgb_best_raw["reg_alpha"]),
    "reg_lambda": float(xgb_best_raw["reg_lambda"]),
}

with open(OUT_DIR / f"best_xgb_params_{TAG}.json", "w") as f:
    json.dump(xgb_best, f, indent=4)

xgb_model = make_xgb_model(xgb_best)


# ---------------------------------------------------------------------
# 5) Final evaluation and saving summary
# ---------------------------------------------------------------------
print("\n[INFO] Final RF evaluation with repeated stratified CV...")
rf_summary = evaluate_final_model(
    rf_model,
    "RF",
    X_train,
    y_train,
    X_test,
    y_test,
    labels_for_report=class_labels
)

print("\n[INFO] Final XGBoost evaluation with repeated stratified CV...")
xgb_summary = evaluate_final_model(
    xgb_model,
    "XGB",
    X_train,
    y_train_enc,
    X_test,
    y_test_enc,
    labels_for_report=class_labels,
    decode_predictions=lambda arr: label_encoder.inverse_transform(np.asarray(arr, dtype=int))
)

summary_df = pd.DataFrame([rf_summary, xgb_summary])
summary_csv = OUT_DIR / f"metrics_summary_{TAG}.csv"
summary_json = OUT_DIR / f"metrics_summary_{TAG}.json"

summary_df.to_csv(summary_csv, index=False)
with open(summary_json, "w") as f:
    json.dump([rf_summary, xgb_summary], f, indent=4)

print("\n[FINAL SUMMARY]")
print(summary_df.to_string(index=False))
print(f"\n[INFO] Results saved in: {OUT_DIR.resolve()}")
