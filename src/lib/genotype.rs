#[must_use]
pub fn copies_of_alt(alleles: &[Option<usize>], alt_index: usize) -> i8 {
    let copies = alleles.iter().filter(|a| **a == Some(alt_index)).count();
    i8::try_from(copies).unwrap_or(i8::MAX)
}

/// `X` never claims more than the reads show, and a missing count cannot cap anything.
#[must_use]
pub fn capped_copies(copies: i8, reads: Option<i32>) -> i8 {
    reads
        .and_then(|reads| i8::try_from(reads).ok())
        .map_or(copies, |reads| copies.min(reads))
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(2, Some(1), 1)]
    #[case(2, Some(2), 2)]
    #[case(2, Some(30), 2)]
    #[case(3, Some(2), 2)]
    #[case(1, Some(0), 0)]
    #[case(1, None, 1)]
    #[case(2, None, 2)]
    #[case(0, Some(4), 0)]
    #[case(1, Some(1000), 1)]
    fn copies_are_capped_at_the_alt_count(
        #[case] copies: i8,
        #[case] reads: Option<i32>,
        #[case] expected: i8,
    ) {
        assert_eq!(capped_copies(copies, reads), expected);
    }

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
