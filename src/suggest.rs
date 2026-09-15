//! Edit-distance typo suggestions ("did you mean sqrt?").

pub fn levenshtein(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    for (i, ca) in a.chars().enumerate() {
        let mut cur = vec![i + 1; b.len() + 1];
        for (j, &cb) in b.iter().enumerate() {
            let cost = usize::from(ca != cb);
            cur[j + 1] = (prev[j] + cost).min(prev[j + 1] + 1).min(cur[j] + 1);
        }
        prev = cur;
    }
    prev[b.len()]
}

/// Closest candidate to `name`, if it's close enough to plausibly be a typo.
/// Never suggests for a distance equal to the name's length (so `x` doesn't
/// become "did you mean e?").
pub fn suggest(name: &str, candidates: &[&str]) -> Option<String> {
    let lower = name.to_lowercase();
    let (dist, best) = candidates
        .iter()
        .map(|&c| (levenshtein(&lower, c), c))
        .min_by_key(|&(d, _)| d)?;
    let len = name.chars().count();
    let limit = if len <= 3 { 1 } else { 2 };
    (dist <= limit && dist < len).then(|| best.to_string())
}
