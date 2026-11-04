#[derive(Clone, Debug, PartialEq)]
pub enum TextColumn {
    Float(Vec<f64>),
    Text(Vec<Option<String>>),
}

#[allow(clippy::needless_pass_by_value)]
pub fn type_column(values: Vec<Option<String>>) -> TextColumn {
    let parts: Vec<Vec<&str>> = values.iter().map(|v| present_parts(v.as_deref())).collect();

    if parts.iter().all(Vec::is_empty) {
        return TextColumn::Text(vec![None; parts.len()]);
    }

    let numbers: Option<Vec<f64>> = parts
        .iter()
        .map(|parts| {
            parts
                .iter()
                .map(|part| part.parse::<f64>().ok().filter(|n| n.is_finite()))
                .try_fold(f64::NAN, |largest, n| n.map(|n| largest.max(n)))
        })
        .collect();

    match numbers {
        Some(numbers) => TextColumn::Float(numbers),
        None => TextColumn::Text(
            parts
                .into_iter()
                .map(|parts| {
                    let mut distinct: Vec<&str> = Vec::new();
                    for part in parts {
                        if !distinct.contains(&part) {
                            distinct.push(part);
                        }
                    }
                    (!distinct.is_empty()).then(|| distinct.join("&"))
                })
                .collect(),
        ),
    }
}

fn present_parts(value: Option<&str>) -> Vec<&str> {
    value
        .map(|v| {
            v.split('&')
                .filter(|p| !p.is_empty() && *p != ".")
                .collect()
        })
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn some(values: &[Option<&str>]) -> Vec<Option<String>> {
        values.iter().map(|v| v.map(str::to_string)).collect()
    }

    fn floats(column: TextColumn) -> Vec<f64> {
        match column {
            TextColumn::Float(values) => values,
            TextColumn::Text(_) => panic!("expected a float column"),
        }
    }

    #[test]
    fn numbers_keep_the_largest_part_and_missing_is_nan() {
        let values = floats(type_column(some(&[Some("0.1&0.7"), None, Some(".")])));
        assert_eq!(values[0], 0.7);
        assert!(values[1].is_nan());
        assert!(values[2].is_nan());
    }

    #[test]
    fn missing_parts_are_dropped_before_typing_numbers() {
        let values = floats(type_column(some(&[Some("0.1&.")])));
        assert_eq!(values, [0.1]);
    }

    #[test]
    fn text_keeps_distinct_parts_in_order() {
        let column = type_column(some(&[Some("A&B&A"), Some("C")]));
        assert_eq!(column, TextColumn::Text(some(&[Some("A&B"), Some("C")])));
    }

    #[test]
    fn text_drops_missing_parts_and_empty_results_are_missing() {
        let column = type_column(some(&[Some("A&."), Some("&"), Some("x")]));
        assert_eq!(
            column,
            TextColumn::Text(some(&[Some("A"), None, Some("x")]))
        );
    }

    #[test]
    fn one_non_number_makes_the_column_text() {
        let column = type_column(some(&[Some("1"), Some("x")]));
        assert_eq!(column, TextColumn::Text(some(&[Some("1"), Some("x")])));
    }

    #[test]
    fn a_column_without_values_is_an_empty_categorical() {
        let column = type_column(some(&[None, Some("."), Some("")]));
        assert_eq!(column, TextColumn::Text(vec![None, None, None]));
    }

    #[test]
    fn non_finite_words_make_the_column_text() {
        for word in ["nan", "inf", "infinity"] {
            let column = type_column(some(&[Some("1"), Some(word)]));
            assert!(matches!(column, TextColumn::Text(_)), "{word}");
        }
    }
}
