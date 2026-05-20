# Stitching fit score

This fitting stage now uses a composite score instead of plain MAE alone.
The goal is to reward not just small energy differences, but also the right
overall shape of the torsion profile.

Lower is better.

## 1. Energy normalization

For a QM profile `E_qm` and MM profile `E_mm`, both arrays are shifted so
their minimum value becomes zero. This removes any constant energy offset and
lets the comparison focus on the profile shape.

$$
E'_{qm} = E_{qm} - \min(E_{qm})
$$

$$
E'_{mm} = E_{mm} - \min(E_{mm})
$$

The peak-to-peak amplitude is the height of the curve from its lowest point to
its highest point:

$$
A_{qm} = \max(E'_{qm}) - \min(E'_{qm})
$$

$$
A_{mm} = \max(E'_{mm}) - \min(E'_{mm})
$$

The QM amplitude is used as the normalization scale for the rest of the score:

$$
S = \max(A_{qm}, \varepsilon)
$$

Here, `\varepsilon` is a tiny positive number that prevents division by zero
when the QM curve is essentially flat.

## 2. Stationary points

Stationary points are detected as local maxima and local minima using peak
finding on the energy curve and on its negation.

$$
\text{stationary points} = \text{peaks}(E) \cup \text{peaks}(-E)
$$

The code works on a 10 degree scan grid, so each index is converted to an
angle using:

$$
\theta_i = 10 \cdot i
$$

## 3. Composite score terms

The final score is the simple average of four equally weighted terms. Each one
looks at a different part of the fit, so a profile has to do well overall
rather than only in one narrow sense.

### a) Stationary-point location score

The QM and MM stationary-point angles are paired in order and compared with a
circular angular distance, so points near `0°` and `360°` are treated as being
close to each other.

$$
d(\theta_{qm}, \theta_{mm}) =
\min\left(
|\theta_{qm} - \theta_{mm}| \bmod 360,\;
360 - \left(|\theta_{qm} - \theta_{mm}| \bmod 360\right)
\right)
$$

The location score is the mean of those distances after normalizing by half a
turn:

$$
L = \operatorname{mean}\left(\frac{d(\theta_{qm}, \theta_{mm})}{180}\right)
$$

This keeps the value roughly in the range `[0, 1]`.

### b) Amplitude score

The amplitude term measures whether the MM curve has the same overall height as
the QM curve:

$$
P = \frac{|A_{qm} - A_{mm}|}{S}
$$

### c) Stationary-point count score

The stationary-point count term checks whether the MM curve has the same number
of peaks and troughs as the QM curve:

$$
C = \frac{|N_{qm} - N_{mm}|}{\max(N_{qm}, N_{mm}, 1)}
$$

where `N_qm` and `N_mm` are the numbers of stationary points.

### d) Normalized MAE score

The MAE term is still the familiar average absolute difference, but it is
normalized by the QM amplitude so that large-energy profiles do not dominate
just because they are larger in magnitude.

$$
M = \frac{\operatorname{mean}(|E'_{qm} - E'_{mm}|)}{S}
$$

## 4. Final composite score

The final score is the average of the four pieces:

$$
\text{Score} = \frac{L + P + C + M}{4}
$$

## 5. How the plot annotations map to the score

Each `fitting_shuffle_*` plot shows two red annotation boxes, one on the total
energy panel and one on the torsion energy panel.

The **total** box reports the current values for the MM total-vs-QM total
comparison:

- `Score` is the composite score used for convergence on the total curve.
- `Loc` is the stationary-point location score.
- `Amp` is the amplitude score.
- `Count` is the stationary-point count score.
- `nMAE` is the normalized MAE score.
- `MAE` is the plain mean absolute error, shown for reference only.

The **torsion** box reports the same fields, but for the fitted torsion curve
after Fourier reconstruction. That means the torsion box describes how well the
current parameter set reproduces the QM torsion profile, while the total box
describes the raw MM total-energy profile used alongside it.

In both boxes, the component values are the ingredients of the composite score,
so a low `Score` means the curve is matching in shape, height, stationary-point
count, and normalized average deviation.

## 6. Convergence tolerance

The tolerance comes from `parameterFittingInfo.converganceTolerance`.

This same tolerance is used for both convergence checks:

1. Per-torsion freezing: a torsion is marked converged when both its torsion score and total score are below the tolerance.
2. Shuffle-level convergence: the run is treated as converged when the RMS score across active torsions is below the tolerance.

If `converganceTolerance` is `None`, convergence checking is disabled and the
run continues until the maximum shuffle limit is reached.
