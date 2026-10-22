//! Non-negative least squares (Lawson–Hanson), shared by the trace decoder
//! (`dna-trace`) and the run-length decision (`dna-consensus`), so both fit
//! by one method.

/// Gradient at or below which no column enters the active set, and value at
/// or below which a passive coefficient leaves it.
const TOLERANCE: f64 = 1e-10;

/// Non-negative least squares by the Lawson–Hanson active-set method over
/// `design`'s columns (each a vector over the observed rows), at most
/// `max_outer` columns entering: the entering column is the first of the
/// largest gradient, and a passive coefficient at or below the tolerance
/// leaves the set. The passive least-squares problems take the minimum-norm
/// solution.
#[must_use]
pub fn nnls(design: &[Vec<f64>], observed: &[f64], max_outer: usize) -> Vec<f64> {
    let columns = design.len();
    let mut solution = vec![0.0; columns];
    let mut passive = vec![false; columns];
    // Each column's non-zero rows: a one-hot design is three-quarters zeros,
    // and a zero adds nothing to a fit or a gradient.
    let sparse = design
        .iter()
        .map(|column| {
            column
                .iter()
                .enumerate()
                .filter(|(_, value)| **value != 0.0)
                .map(|(row, &value)| (row, value))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let gradient = |solution: &[f64]| -> Vec<f64> {
        // The fit, column by column in column order, over the positive
        // fractions only.
        let mut fitted = vec![-0.0; observed.len()];
        for (column, &x) in sparse.iter().zip(solution) {
            if x != 0.0 {
                for &(row, value) in column {
                    fitted[row] += value * x;
                }
            }
        }
        let residual = observed
            .iter()
            .zip(&fitted)
            .map(|(y, fit)| y - fit)
            .collect::<Vec<_>>();
        sparse
            .iter()
            .map(|column| {
                column
                    .iter()
                    .map(|&(row, value)| value * residual[row])
                    .sum()
            })
            .collect()
    };
    let mut factor = Factor::default();
    let mut descent = gradient(&solution);
    let mut outer = 0;
    while outer < max_outer {
        let Some(entering) =
            (0..columns)
                .filter(|&column| !passive[column])
                .fold(None::<usize>, |best, column| match best {
                    Some(kept) if descent[column] <= descent[kept] => Some(kept),
                    _ => Some(column),
                })
        else {
            break;
        };
        if descent[entering] <= TOLERANCE {
            break;
        }
        outer += 1;
        passive[entering] = true;
        loop {
            let chosen = (0..columns)
                .filter(|&column| passive[column])
                .collect::<Vec<_>>();
            let least = least_squares(design, observed, &chosen, &mut factor);
            let mut trial = vec![0.0; columns];
            for (&column, &value) in chosen.iter().zip(&least) {
                trial[column] = value;
            }
            if chosen.iter().all(|&column| trial[column] > TOLERANCE) {
                solution = trial;
                break;
            }
            let step = chosen
                .iter()
                .filter(|&&column| trial[column] <= TOLERANCE)
                .map(|&column| solution[column] / (solution[column] - trial[column]))
                .fold(f64::INFINITY, f64::min);
            for column in 0..columns {
                solution[column] += step * (trial[column] - solution[column]);
            }
            for column in 0..columns {
                passive[column] = passive[column] && solution[column] > TOLERANCE;
                if !passive[column] {
                    solution[column] = 0.0;
                }
            }
            if !passive.iter().any(|&flag| flag) {
                break;
            }
        }
        descent = gradient(&solution);
    }
    solution
}

/// Minimum-norm least squares over `chosen` columns of `design`, by a
/// one-sided Jacobi singular value decomposition: singular values below
/// machine epsilon times the larger dimension times the largest one count as
/// zero, as LAPACK's minimum-norm solver treats them. The decode's columns
/// can be linearly dependent, so the minimum-norm solution is the one that
/// fixes the active set.
fn least_squares(
    design: &[Vec<f64>],
    observed: &[f64],
    chosen: &[usize],
    factor: &mut Factor,
) -> Vec<f64> {
    let size = chosen.len();
    if size == 0 {
        return Vec::new();
    }
    // A full-rank passive set has one least-squares solution, which QR finds
    // fast; only a rank-deficient one needs the minimum-norm solve.
    if let Some(solution) = full_rank_solution(design, observed, chosen, factor) {
        return solution;
    }
    let mut left = chosen
        .iter()
        .map(|&column| design[column].clone())
        .collect::<Vec<_>>();
    let mut right = (0..size)
        .map(|i| {
            (0..size)
                .map(|j| if i == j { 1.0 } else { 0.0 })
                .collect::<Vec<f64>>()
        })
        .collect::<Vec<_>>();
    for _ in 0..80 {
        let mut rotated = false;
        for p in 0..size {
            for q in p + 1..size {
                let alpha = left[p].iter().map(|value| value * value).sum::<f64>();
                let beta = left[q].iter().map(|value| value * value).sum::<f64>();
                let gamma = left[p]
                    .iter()
                    .zip(&left[q])
                    .map(|(a, b)| a * b)
                    .sum::<f64>();
                if gamma.abs() <= f64::EPSILON * (alpha * beta).sqrt() || gamma == 0.0 {
                    continue;
                }
                rotated = true;
                let zeta = (beta - alpha) / (2.0 * gamma);
                let tangent = zeta.signum() / (zeta.abs() + (1.0 + zeta * zeta).sqrt());
                let cos_theta = 1.0 / (1.0 + tangent * tangent).sqrt();
                let sin_theta = cos_theta * tangent;
                for matrix in [&mut left, &mut right] {
                    let (first, second) = matrix.split_at_mut(q);
                    for (x, y) in first[p].iter_mut().zip(second[0].iter_mut()) {
                        let (old_x, old_y) = (*x, *y);
                        *x = cos_theta * old_x - sin_theta * old_y;
                        *y = sin_theta * old_x + cos_theta * old_y;
                    }
                }
            }
        }
        if !rotated {
            break;
        }
    }
    let values = left
        .iter()
        .map(|column| column.iter().map(|value| value * value).sum::<f64>().sqrt())
        .collect::<Vec<_>>();
    let largest = values.iter().copied().fold(0.0, f64::max);
    let cutoff = f64::EPSILON * observed.len().max(size) as f64 * largest;
    let mut solution = vec![0.0; size];
    for (k, &value) in values.iter().enumerate() {
        if value <= cutoff || value == 0.0 {
            continue;
        }
        let projection = left[k]
            .iter()
            .zip(observed)
            .map(|(a, y)| a * y)
            .sum::<f64>()
            / (value * value);
        for (cell, v) in solution.iter_mut().zip(&right[k]) {
            *cell += projection * v;
        }
    }
    solution
}

/// A modified Gram–Schmidt QR of some columns, in order, kept between
/// solves: the factor of a column depends only on the columns before it, so
/// a solve whose columns share a prefix with the last one's reuses it.
#[derive(Default)]
struct Factor {
    /// The columns factored, in order.
    columns: Vec<usize>,
    /// The orthonormal basis, one vector per column.
    basis: Vec<Vec<f64>>,
    /// Column `k` of the upper triangle: its `k + 1` entries.
    upper: Vec<Vec<f64>>,
}

/// The least-squares solution over `chosen` columns by modified Gram–Schmidt
/// QR, `None` when a column is, within the rank cutoff of [`least_squares`],
/// a combination of the ones before it.
fn full_rank_solution(
    design: &[Vec<f64>],
    observed: &[f64],
    chosen: &[usize],
    factor: &mut Factor,
) -> Option<Vec<f64>> {
    let size = chosen.len();
    let largest = chosen
        .iter()
        .map(|&column| {
            design[column]
                .iter()
                .map(|value| value * value)
                .sum::<f64>()
                .sqrt()
        })
        .fold(0.0, f64::max);
    let cutoff = 1e-8 * largest.max(f64::MIN_POSITIVE);
    let shared = factor
        .columns
        .iter()
        .zip(chosen)
        .take_while(|(kept, wanted)| kept == wanted)
        .count();
    factor.columns.truncate(shared);
    factor.basis.truncate(shared);
    factor.upper.truncate(shared);
    if factor
        .upper
        .iter()
        .enumerate()
        .any(|(k, column)| column[k] <= cutoff)
    {
        return None;
    }
    for &column in &chosen[shared..] {
        let mut vector = design[column].clone();
        let mut entries = Vec::with_capacity(factor.basis.len() + 1);
        for q in &factor.basis {
            let dot = q.iter().zip(&vector).map(|(a, b)| a * b).sum::<f64>();
            entries.push(dot);
            for (cell, value) in vector.iter_mut().zip(q) {
                *cell -= dot * value;
            }
        }
        let norm = vector.iter().map(|value| value * value).sum::<f64>().sqrt();
        entries.push(norm);
        if norm <= cutoff {
            return None;
        }
        for cell in &mut vector {
            *cell /= norm;
        }
        factor.columns.push(column);
        factor.basis.push(vector);
        factor.upper.push(entries);
    }
    let projected = factor
        .basis
        .iter()
        .map(|q| q.iter().zip(observed).map(|(a, y)| a * y).sum::<f64>())
        .collect::<Vec<_>>();
    let mut solution = vec![0.0; size];
    for k in (0..size).rev() {
        let mut sum = projected[k];
        for (column, value) in factor.upper.iter().zip(&solution).skip(k + 1) {
            sum -= column[k] * value;
        }
        solution[k] = sum / factor.upper[k][k];
    }
    Some(solution)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(left: &[f64], right: &[f64]) -> bool {
        left.len() == right.len() && left.iter().zip(right).all(|(a, b)| (a - b).abs() <= 1e-9)
    }

    #[test]
    fn recovers_a_non_negative_combination() {
        let design = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![1.0, 1.0, 1.0],
        ];
        let solution = nnls(&design, &[5.0, 3.0, 3.0], 10);
        assert!(close(&solution, &[2.0, 0.0, 3.0]), "{solution:?}");
    }

    #[test]
    fn holds_a_negative_least_squares_coefficient_at_zero() {
        // Unconstrained, the fit is 3 of the first column and -1 of the
        // second; the best non-negative fit drops the second.
        let design = vec![vec![1.0, 0.0], vec![1.0, 1.0]];
        let solution = nnls(&design, &[2.0, -1.0], 10);
        assert!(close(&solution, &[2.0, 0.0]), "{solution:?}");
        // An observation no column correlates with positively fits nothing.
        assert!(close(&nnls(&[vec![1.0, 1.0]], &[-1.0, -1.0], 10), &[0.0]));
    }

    #[test]
    fn enters_at_most_max_outer_columns() {
        let design = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        assert!(close(&nnls(&design, &[1.0, 2.0], 0), &[0.0, 0.0]));
        // The column of the largest gradient enters first.
        assert!(close(&nnls(&design, &[1.0, 2.0], 1), &[0.0, 2.0]));
        assert!(close(&nnls(&design, &[1.0, 2.0], 2), &[1.0, 2.0]));
        assert!(nnls(&[], &[1.0], 5).is_empty());
    }

    #[test]
    fn a_tie_enters_the_first_column() {
        let design = vec![vec![1.0, 0.0], vec![1.0, 0.0]];
        assert!(close(&nnls(&design, &[2.0, 0.0], 10), &[2.0, 0.0]));
    }

    #[test]
    fn a_rank_deficient_passive_set_takes_the_minimum_norm_solution() {
        let design = vec![vec![1.0, 0.0], vec![1.0, 0.0]];
        let mut factor = Factor::default();
        assert!(full_rank_solution(&design, &[2.0, 0.0], &[0, 1], &mut factor).is_none());
        let solution = least_squares(&design, &[2.0, 0.0], &[0, 1], &mut Factor::default());
        assert!(close(&solution, &[1.0, 1.0]), "{solution:?}");
        assert!(least_squares(&design, &[2.0, 0.0], &[], &mut factor).is_empty());
    }

    #[test]
    fn a_reused_factor_solves_as_a_fresh_one() {
        let design = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![1.0, 1.0, 1.0],
        ];
        let observed = [5.0, 3.0, 3.0];
        let mut factor = Factor::default();
        let first = full_rank_solution(&design, &observed, &[0, 2], &mut factor);
        let reused = full_rank_solution(&design, &observed, &[0, 1, 2], &mut factor);
        let fresh = full_rank_solution(&design, &observed, &[0, 1, 2], &mut Factor::default());
        assert!(first.is_some());
        assert_eq!(reused, fresh);
        assert!(close(&fresh.unwrap_or_default(), &[2.0, 0.0, 3.0]));
    }
}
