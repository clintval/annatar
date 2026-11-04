#[must_use]
pub fn genotype_code(alleles: &[Option<usize>], alt_index: usize) -> Option<i8> {
    if alleles.is_empty() || alleles.iter().any(Option::is_none) {
        return None;
    }
    let copies = alleles.iter().filter(|a| **a == Some(alt_index)).count();
    Some(match copies {
        0 => 1,
        1 if alleles.len() == 1 => 3,
        1 => 2,
        _ => 3,
    })
}

#[must_use]
pub fn covered(gt_called: bool, ad: &[Option<i32>]) -> bool {
    gt_called || ad.iter().flatten().sum::<i32>() > 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genotype_codes() {
        assert_eq!(genotype_code(&[Some(0), Some(0)], 1), Some(1));
        assert_eq!(genotype_code(&[Some(0), Some(1)], 1), Some(2));
        assert_eq!(genotype_code(&[Some(1), Some(1)], 1), Some(3));
        assert_eq!(genotype_code(&[Some(0), Some(2)], 1), Some(1));
        assert_eq!(genotype_code(&[Some(1)], 1), Some(3));
        assert_eq!(
            genotype_code(&[Some(0), Some(0), Some(0), Some(1)], 1),
            Some(2)
        );
        assert_eq!(genotype_code(&[None, None], 1), None);
        assert_eq!(genotype_code(&[Some(1), None], 1), None);
        assert_eq!(genotype_code(&[], 1), None);
    }

    #[test]
    fn covered_by_gt_or_depth() {
        assert!(covered(true, &[]));
        assert!(covered(false, &[Some(5), None]));
        assert!(!covered(false, &[Some(0), None]));
        assert!(!covered(false, &[]));
    }
}
