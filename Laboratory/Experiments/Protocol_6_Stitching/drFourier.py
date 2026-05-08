import numpy as np
import pandas as pd
from os import path as p


## dummy classes
class FilePath:
    pass


class DirPath:
    pass


def _build_harmonic_design_matrix(n_points: int, max_functions: int) -> np.ndarray:
    angles = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
    cols = []
    for periodicity in range(1, max_functions + 1):
        cols.append(np.cos(periodicity * angles))
        cols.append(np.sin(periodicity * angles))
    return np.column_stack(cols)


def _huber_weights(residuals: np.ndarray, delta: float) -> np.ndarray:
    abs_residuals = np.abs(residuals)
    weights = np.ones_like(abs_residuals)
    mask = abs_residuals > delta
    weights[mask] = delta / abs_residuals[mask]
    return weights


def _group_soft_threshold(beta: np.ndarray, threshold: float) -> np.ndarray:
    out = np.zeros_like(beta)
    n_groups = len(beta) // 2
    for group_idx in range(n_groups):
        i = 2 * group_idx
        group = beta[i:i + 2]
        group_norm = np.linalg.norm(group)
        if group_norm <= threshold:
            continue
        out[i:i + 2] = (1.0 - threshold / group_norm) * group
    return out


def _solve_weighted_group_lasso(
    Xw: np.ndarray,
    yw: np.ndarray,
    lambda_value: float,
    max_iterations: int,
    convergence_tol: float,
    warm_start: np.ndarray | None = None,
) -> np.ndarray:
    if warm_start is None:
        beta = np.zeros(Xw.shape[1], dtype=float)
    else:
        beta = warm_start.copy()

    spectral_norm = np.linalg.norm(Xw, ord=2)
    lipschitz = max(spectral_norm * spectral_norm, 1e-8)
    step = 1.0 / lipschitz

    for _ in range(max_iterations):
        grad = Xw.T @ (Xw @ beta - yw)
        candidate = beta - step * grad
        updated = _group_soft_threshold(candidate, step * lambda_value)
        delta = np.linalg.norm(updated - beta)
        scale = max(np.linalg.norm(beta), 1.0)
        beta = updated
        if delta / scale < convergence_tol:
            break
    return beta


def _robust_group_lasso_fit(
    X: np.ndarray,
    y: np.ndarray,
    lambda_value: float,
    huber_delta: float,
    max_outer_iterations: int,
    max_inner_iterations: int,
    convergence_tol: float,
    warm_start: np.ndarray | None = None,
) -> np.ndarray:
    beta = np.zeros(X.shape[1], dtype=float) if warm_start is None else warm_start.copy()
    for _ in range(max_outer_iterations):
        residuals = y - X @ beta
        weights = _huber_weights(residuals, huber_delta)
        sqrt_w = np.sqrt(weights)
        Xw = X * sqrt_w[:, None]
        yw = y * sqrt_w
        updated = _solve_weighted_group_lasso(
            Xw, yw, lambda_value, max_inner_iterations, convergence_tol, warm_start=beta
        )
        delta = np.linalg.norm(updated - beta)
        scale = max(np.linalg.norm(beta), 1.0)
        beta = updated
        if delta / scale < convergence_tol:
            break
    return beta


def _robust_refit_active_groups(
    X: np.ndarray,
    y: np.ndarray,
    active_groups: list[int],
    huber_delta: float,
    max_iterations: int,
    convergence_tol: float,
) -> np.ndarray:
    beta = np.zeros(X.shape[1], dtype=float)
    if len(active_groups) == 0:
        return beta

    active_cols = []
    for group_idx in active_groups:
        i = 2 * group_idx
        active_cols.extend([i, i + 1])

    Xa = X[:, active_cols]
    beta_active = np.zeros(len(active_cols), dtype=float)

    for _ in range(max_iterations):
        residuals = y - Xa @ beta_active
        weights = _huber_weights(residuals, huber_delta)
        W = weights[:, None]
        lhs = Xa.T @ (W * Xa) + 1e-8 * np.eye(Xa.shape[1])
        rhs = Xa.T @ (weights * y)
        updated = np.linalg.solve(lhs, rhs)
        delta = np.linalg.norm(updated - beta_active)
        scale = max(np.linalg.norm(beta_active), 1.0)
        beta_active = updated
        if delta / scale < convergence_tol:
            break

    beta[active_cols] = beta_active
    return beta


def _select_active_groups(beta: np.ndarray, threshold: float = 1e-6) -> list[int]:
    active = []
    n_groups = len(beta) // 2
    for group_idx in range(n_groups):
        i = 2 * group_idx
        if np.linalg.norm(beta[i:i + 2]) > threshold:
            active.append(group_idx)
    return active


def _reconstruct_signal_from_beta(X: np.ndarray, beta: np.ndarray) -> np.ndarray:
    return X @ beta


def _bic_score(observed: np.ndarray, predicted: np.ndarray, n_active_groups: int) -> float:
    n_points = len(observed)
    rss = np.sum((observed - predicted) ** 2)
    rss = max(rss, 1e-12)
    n_params = max(2 * n_active_groups, 1)
    return n_points * np.log(rss / n_points) + n_params * np.log(n_points)


def _params_from_beta(beta: np.ndarray) -> pd.DataFrame:
    rows = []
    n_groups = len(beta) // 2
    for group_idx in range(n_groups):
        i = 2 * group_idx
        cos_coef = beta[i]
        sin_coef = beta[i + 1]
        amplitude = float(np.sqrt(cos_coef * cos_coef + sin_coef * sin_coef))
        if amplitude <= 1e-8:
            continue
        phase = float(np.degrees(np.arctan2(sin_coef, cos_coef)))
        rows.append({
            "Amplitude": amplitude,
            "Period": group_idx + 1,
            "Phase": phase,
        })

    if len(rows) == 0:
        rows = [{"Amplitude": 0.0, "Period": 1, "Phase": 0.0}]

    param_df = pd.DataFrame(rows)
    param_df.sort_values(by="Amplitude", ascending=False, inplace=True, ignore_index=True)
    return param_df


def _cosine_components_from_params(param_df: pd.DataFrame, n_points: int) -> dict:
    angle = np.radians(np.linspace(0, 360, n_points, endpoint=False))
    components = {}
    for _, row in param_df.iterrows():
        amplitude = float(row["Amplitude"])
        periodicity = abs(float(row["Period"]))
        phase = np.radians(float(row["Phase"]))
        component = amplitude * (1 + np.cos(periodicity * angle - phase))
        components[periodicity] = component
    return components


def _signal_from_params(param_df: pd.DataFrame, n_points: int) -> np.ndarray:
    angle = np.radians(np.linspace(0, 360, n_points, endpoint=False))
    signal = np.zeros_like(angle)
    for _, row in param_df.iterrows():
        amplitude = float(row["Amplitude"])
        periodicity = abs(float(row["Period"]))
        phase = np.radians(float(row["Phase"]))
        signal += amplitude * (1 + np.cos(periodicity * angle - phase))
    signal = signal - signal.min()
    return signal


def _fit_active_groups_least_squares(X: np.ndarray, y: np.ndarray, active_groups: list[int], ridge: float = 1e-8) -> np.ndarray:
    beta = np.zeros(X.shape[1], dtype=float)
    if len(active_groups) == 0:
        return beta

    active_cols = []
    for group_idx in active_groups:
        i = 2 * group_idx
        active_cols.extend([i, i + 1])

    Xa = X[:, active_cols]
    lhs = Xa.T @ Xa + ridge * np.eye(Xa.shape[1])
    rhs = Xa.T @ y
    beta_active = np.linalg.solve(lhs, rhs)
    beta[active_cols] = beta_active
    return beta


def _stabilize_overfit_beta(
    beta: np.ndarray,
    qm_amplitude: float,
    max_term_amplitude_factor: float,
    max_total_amplitude_factor: float,
) -> np.ndarray:
    stabilized = beta.copy()
    term_cap = max(1e-8, max_term_amplitude_factor * qm_amplitude)
    total_cap = max(1e-8, max_total_amplitude_factor * qm_amplitude)

    group_norms = []
    for group_idx in range(len(stabilized) // 2):
        i = 2 * group_idx
        group = stabilized[i:i + 2]
        group_norm = float(np.linalg.norm(group))
        if group_norm > term_cap:
            stabilized[i:i + 2] = group * (term_cap / group_norm)
            group_norm = term_cap
        group_norms.append(group_norm)

    current_total = sum(group_norms)
    if current_total > total_cap:
        stabilized *= total_cap / current_total
    return stabilized


def _overfit_prune_harmonic_fit(qm_torsion_energy: np.ndarray, max_cosine_functions: int, config: dict) -> tuple[pd.DataFrame, dict, float]:
    fit_cfg = config["parameterFittingInfo"]
    n_points = len(qm_torsion_energy)
    max_possible_functions = max(1, n_points // 2)
    initial_functions = fit_cfg.get("overfitInitialFunctions", 12)
    initial_functions = max(max_cosine_functions, initial_functions)
    initial_functions = min(initial_functions, max_possible_functions)

    min_functions = fit_cfg.get("overfitMinFunctions", 1)
    min_functions = max(1, min(min_functions, max_cosine_functions))
    prune_tolerance = fit_cfg.get("overfitPruneMaeIncreaseTolerance", 0.02)
    ridge = fit_cfg.get("overfitRidge", 1e-4)
    max_term_amplitude_factor = fit_cfg.get("overfitMaxTermAmplitudeFactor", 1.5)
    max_total_amplitude_factor = fit_cfg.get("overfitMaxTotalAmplitudeFactor", 3.0)

    y = qm_torsion_energy - np.min(qm_torsion_energy)
    y_centered = y - np.mean(y)
    qm_amplitude = float(np.max(y) - np.min(y))
    X = _build_harmonic_design_matrix(n_points, initial_functions)

    active_groups = list(range(initial_functions))
    beta = _fit_active_groups_least_squares(X, y_centered, active_groups, ridge=ridge)
    beta = _stabilize_overfit_beta(
        beta,
        qm_amplitude=qm_amplitude,
        max_term_amplitude_factor=max_term_amplitude_factor,
        max_total_amplitude_factor=max_total_amplitude_factor,
    )
    current_mae = float(np.mean(np.abs(y_centered - X @ beta)))

    while len(active_groups) > min_functions:
        best_group_to_remove = None
        best_beta = None
        best_mae = np.inf
        best_delta_mae = np.inf

        for group_idx in active_groups:
            candidate_groups = [g for g in active_groups if g != group_idx]
            candidate_beta = _fit_active_groups_least_squares(X, y_centered, candidate_groups, ridge=ridge)
            candidate_beta = _stabilize_overfit_beta(
                candidate_beta,
                qm_amplitude=qm_amplitude,
                max_term_amplitude_factor=max_term_amplitude_factor,
                max_total_amplitude_factor=max_total_amplitude_factor,
            )
            candidate_mae = float(np.mean(np.abs(y_centered - X @ candidate_beta)))
            delta_mae = candidate_mae - current_mae
            if delta_mae < best_delta_mae:
                best_delta_mae = delta_mae
                best_mae = candidate_mae
                best_group_to_remove = group_idx
                best_beta = candidate_beta

        if best_group_to_remove is None:
            break

        if best_delta_mae <= prune_tolerance or len(active_groups) > max_cosine_functions:
            active_groups.remove(best_group_to_remove)
            beta = best_beta
            current_mae = best_mae
        else:
            break

    param_df = _params_from_beta(beta)
    if len(param_df) > max_cosine_functions:
        param_df = param_df.sort_values(by="Amplitude", ascending=False).head(max_cosine_functions).reset_index(drop=True)

    cosine_components = _cosine_components_from_params(param_df, n_points)
    reconstructed_signal = _signal_from_params(param_df, n_points)
    mean_average_error = float(np.mean(np.abs(reconstructed_signal - y)))
    return param_df, cosine_components, mean_average_error


def _sparse_harmonic_fit(qm_torsion_energy: np.ndarray, max_cosine_functions: int, config: dict) -> tuple[pd.DataFrame, dict, float]:
    fit_cfg = config["parameterFittingInfo"]
    l2_damping = fit_cfg["l2DampingFactor"] if fit_cfg["l2DampingFactor"] is not None else 0.0
    huber_delta = fit_cfg.get("robustHuberDelta", 1.0)
    lambda_count = fit_cfg.get("robustLambdaCount", 25)
    lambda_ratio = fit_cfg.get("robustLambdaRatio", 0.02)
    max_outer_iterations = fit_cfg.get("robustOuterIterations", 15)
    max_inner_iterations = fit_cfg.get("robustInnerIterations", 300)
    convergence_tol = fit_cfg.get("robustConvergenceTolerance", 1e-6)

    y = qm_torsion_energy - np.min(qm_torsion_energy)
    y_centered = y - np.mean(y)
    n_points = len(y_centered)
    X = _build_harmonic_design_matrix(n_points, max_cosine_functions)

    grad_at_zero = X.T @ y_centered
    lambda_max = 0.0
    for group_idx in range(max_cosine_functions):
        i = 2 * group_idx
        group_norm = np.linalg.norm(grad_at_zero[i:i + 2])
        lambda_max = max(lambda_max, group_norm)
    lambda_max = max(lambda_max, 1e-6)
    lambda_min = lambda_max * lambda_ratio
    lambda_path = np.geomspace(lambda_max, lambda_min, num=lambda_count)

    best_beta = np.zeros(X.shape[1], dtype=float)
    best_bic = np.inf
    best_n_groups = np.inf
    warm_beta = np.zeros(X.shape[1], dtype=float)

    for lambda_value in lambda_path:
        beta_lasso = _robust_group_lasso_fit(
            X,
            y_centered,
            lambda_value,
            huber_delta,
            max_outer_iterations=max_outer_iterations,
            max_inner_iterations=max_inner_iterations,
            convergence_tol=convergence_tol,
            warm_start=warm_beta,
        )
        warm_beta = beta_lasso
        active_groups = _select_active_groups(beta_lasso)
        beta_refit = _robust_refit_active_groups(
            X, y_centered, active_groups, huber_delta, max_iterations=max_outer_iterations, convergence_tol=convergence_tol
        )
        if l2_damping > 0:
            beta_refit = beta_refit / (1 + l2_damping * np.abs(beta_refit))

        reconstructed = _reconstruct_signal_from_beta(X, beta_refit)
        bic = _bic_score(y_centered, reconstructed, len(active_groups))
        if bic < best_bic or (np.isclose(bic, best_bic) and len(active_groups) < best_n_groups):
            best_bic = bic
            best_n_groups = len(active_groups)
            best_beta = beta_refit

    param_df = _params_from_beta(best_beta)
    cosine_components = _cosine_components_from_params(param_df, n_points)
    reconstructed_signal = _signal_from_params(param_df, n_points)
    mean_average_error = float(np.mean(np.abs(reconstructed_signal - y)))
    return param_df, cosine_components, mean_average_error


def fourier_transform_protocol(qmTorsionEnergy, torsionTag, torsionFittingDir, config, maxCosineFunctions: int | None = None):
    configured_max = config["parameterFittingInfo"]["maxCosineFunctions"]
    if maxCosineFunctions is None:
        maxCosineFunctions = configured_max
    else:
        maxCosineFunctions = min(maxCosineFunctions, configured_max)

    fitting_protocol = config["parameterFittingInfo"].get("fittingProtocol", "ROBUST_SPARSE_HARMONIC")
    if fitting_protocol == "OVERFIT_PRUNE_HARMONIC":
        param_df, cosine_components, mean_average_error = _overfit_prune_harmonic_fit(
            qmTorsionEnergy, maxCosineFunctions, config
        )
    else:
        param_df, cosine_components, mean_average_error = _sparse_harmonic_fit(
            qmTorsionEnergy, maxCosineFunctions, config
        )

    out_csv: FilePath = p.join(torsionFittingDir, f"{torsionTag}.csv")
    param_df.to_csv(out_csv)

    return param_df, cosine_components, mean_average_error


if __name__ == "__main__":
    raise NotImplementedError
