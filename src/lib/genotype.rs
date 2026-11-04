#[must_use]
pub fn copies_of_alt(alleles: &[Option<usize>], alt_index: usize) -> i8 {
    let copies = alleles.iter().filter(|a| **a == Some(alt_index)).count();
    i8::try_from(copies).unwrap_or(i8::MAX)
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(&[Some(0), Some(1)], 1)]
    #[case(&[None, Some(1)], 1)]
    #[case(&[Some(1), None], 1)]
    #[case(&[Some(1), Some(1)], 2)]
    #[case(&[Some(1), Some(1), Some(1)], 3)]
    #[case(&[Some(1)], 1)]
    #[case(&[Some(0), Some(0), Some(0), Some(1)], 1)]
    #[case(&[Some(0), Some(0)], 0)]
    #[case(&[Some(0), Some(2)], 0)]
    #[case(&[Some(0)], 0)]
    #[case(&[None, None], 0)]
    #[case(&[], 0)]
    fn copies_count_this_alt_in_the_call(#[case] alleles: &[Option<usize>], #[case] copies: i8) {
        assert_eq!(copies_of_alt(alleles, 1), copies);
    }

    #[test]
    fn copies_saturate_within_the_stored_type() {
        assert_eq!(copies_of_alt(&[Some(1); 200], 1), i8::MAX);
    }
}
