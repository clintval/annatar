use spillsort::{Config, Sorter, key};

/// The positions of variants in order of contig rank, start, end, REF and ALT.
pub fn variant_order(
    ranks: &[i64],
    starts: &[i64],
    ends: &[i64],
    references: &[String],
    alternates: &[String],
) -> Result<Vec<usize>, spillsort::Error> {
    let mut sorter = Sorter::new(Config::default())?;
    let mut sort_key = Vec::new();
    let spans = ranks.iter().zip(starts).zip(ends);
    let variants = spans.zip(references).zip(alternates);
    for (row, ((((rank, start), end), reference), alternate)) in variants.enumerate() {
        sort_key.clear();
        key::push_int(&mut sort_key, *rank);
        key::push_int(&mut sort_key, *start);
        key::push_int(&mut sort_key, *end);
        key::push_str(&mut sort_key, reference);
        key::push_str(&mut sort_key, alternate);
        sorter.push(&sort_key, &row.to_le_bytes())?;
    }
    let mut merged = sorter.finish()?;
    let mut order = Vec::with_capacity(ranks.len());
    while let Some((_, value)) = merged.next_record()? {
        let row = value.try_into().expect("a value is one row number");
        order.push(usize::from_le_bytes(row));
    }
    Ok(order)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn text(values: &[&str]) -> Vec<String> {
        values.iter().map(ToString::to_string).collect()
    }

    #[test]
    fn rank_comes_before_start_then_end_then_ref_then_alt() {
        let order = variant_order(
            &[1, 0, 1, 1, 1],
            &[5, 100, 5, 5, 20],
            &[6, 101, 7, 6, 21],
            &text(&["C", "A", "AT", "C", "A"]),
            &text(&["T", "G", "A", "A", "G"]),
        )
        .unwrap();
        assert_eq!(order, [1, 3, 0, 2, 4]);
    }

    #[test]
    fn coordinates_compare_as_numbers() {
        let order = variant_order(
            &[0, 0, 0],
            &[1000, 20, 100],
            &[1001, 21, 101],
            &text(&["A", "A", "A"]),
            &text(&["G", "G", "G"]),
        )
        .unwrap();
        assert_eq!(order, [1, 2, 0]);
    }

    #[test]
    fn nothing_sorts_to_nothing() {
        assert_eq!(
            variant_order(&[], &[], &[], &[], &[]).unwrap(),
            Vec::<usize>::new()
        );
    }
}
