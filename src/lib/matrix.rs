use nalgebra_sparse::{CscMatrix, CsrMatrix};

use crate::fields::Kind;
use crate::h5ad::LayerMatrix;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Cell {
    Int(Option<i32>),
    Float(Option<f32>),
}

pub struct Matrices {
    pub x: CsrMatrix<i8>,
    pub layers: Vec<(String, LayerMatrix)>,
    pub missing: Vec<(String, u64)>,
}

struct Columns<T> {
    indptr: Vec<usize>,
    indices: Vec<usize>,
    data: Vec<T>,
}

impl<T: nalgebra_sparse::na::Scalar> Columns<T> {
    fn new() -> Self {
        Self {
            indptr: vec![0],
            indices: Vec::new(),
            data: Vec::new(),
        }
    }

    fn push(&mut self, row: u32, value: T) {
        self.indices.push(row as usize);
        self.data.push(value);
    }

    fn end_column(&mut self) {
        self.indptr.push(self.data.len());
    }

    fn into_csr(self, n_obs: usize) -> CsrMatrix<T> {
        let n_vars = self.indptr.len() - 1;
        let csc = CscMatrix::try_from_csc_data(n_obs, n_vars, self.indptr, self.indices, self.data)
            .expect("rows must be sorted, unique and below n_obs within each column");
        CsrMatrix::from(&csc)
    }
}

enum LayerColumns {
    Int(Columns<i32>),
    Float(Columns<f32>),
}

struct Layer {
    name: String,
    columns: LayerColumns,
    missing: u64,
}

pub struct SparseBuilder {
    n_obs: usize,
    x: Columns<i8>,
    layers: Vec<Layer>,
}

impl SparseBuilder {
    /// # Panics
    /// Panics if a layer is neither `Kind::Integer` nor `Kind::Float`.
    #[must_use]
    pub fn new(n_obs: usize, layers: &[(String, Kind)]) -> Self {
        let layers = layers
            .iter()
            .map(|(name, kind)| Layer {
                name: name.clone(),
                columns: match kind {
                    Kind::Integer => LayerColumns::Int(Columns::new()),
                    Kind::Float => LayerColumns::Float(Columns::new()),
                    Kind::Flag | Kind::Text => panic!("layer {name} must be integer or float"),
                },
                missing: 0,
            })
            .collect();
        Self {
            n_obs,
            x: Columns::new(),
            layers,
        }
    }

    /// # Panics
    /// Panics if `cells` does not match the layers, or a cell's kind differs from its layer's.
    pub fn push_column(&mut self, rows: &[u32], x: &[Option<i8>], cells: &[Vec<Cell>]) {
        assert_eq!(cells.len(), self.layers.len(), "one cell vector per layer");
        for (&row, code) in rows.iter().zip(x) {
            if let Some(code) = code {
                self.x.push(row, *code);
            }
        }
        self.x.end_column();
        for (layer, column) in self.layers.iter_mut().zip(cells) {
            assert_eq!(column.len(), rows.len(), "one cell per row");
            for (&row, cell) in rows.iter().zip(column) {
                match (&mut layer.columns, cell) {
                    (LayerColumns::Int(c), Cell::Int(Some(v))) => c.push(row, *v),
                    (LayerColumns::Float(c), Cell::Float(Some(v))) => c.push(row, *v),
                    (LayerColumns::Int(_), Cell::Int(None)) => layer.missing += 1,
                    (LayerColumns::Float(c), Cell::Float(None)) => {
                        c.push(row, f32::NAN);
                        layer.missing += 1;
                    }
                    _ => panic!("cell kind does not match layer {}", layer.name),
                }
            }
            match &mut layer.columns {
                LayerColumns::Int(c) => c.end_column(),
                LayerColumns::Float(c) => c.end_column(),
            }
        }
    }

    #[must_use]
    pub fn finish(self) -> Matrices {
        let n_obs = self.n_obs;
        let mut layers = Vec::with_capacity(self.layers.len());
        let mut missing = Vec::new();
        for layer in self.layers {
            if layer.missing > 0 {
                missing.push((layer.name.clone(), layer.missing));
            }
            let matrix = match layer.columns {
                LayerColumns::Int(c) => LayerMatrix::I32(c.into_csr(n_obs)),
                LayerColumns::Float(c) => LayerMatrix::F32(c.into_csr(n_obs)),
            };
            layers.push((layer.name, matrix));
        }
        Matrices {
            x: self.x.into_csr(n_obs),
            layers,
            missing,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn int_layers() -> Vec<(String, Kind)> {
        vec![("AD_alt".to_string(), Kind::Integer)]
    }

    fn csr_entry<T: Copy>(m: &CsrMatrix<T>, row: usize, col: usize) -> Option<T> {
        let r = m.row(row);
        let at = r.col_indices().iter().position(|&c| c == col)?;
        Some(r.values()[at])
    }

    #[test]
    fn sparse_builder_keeps_explicit_zeros() {
        let mut b = SparseBuilder::new(3, &int_layers());
        b.push_column(
            &[0, 2],
            &[Some(1), Some(2)],
            &[vec![Cell::Int(Some(0)), Cell::Int(Some(4))]],
        );
        let m = b.finish();
        let LayerMatrix::I32(ad) = &m.layers[0].1 else {
            panic!("wrong kind")
        };
        assert_eq!(ad.nnz(), 2);
        assert_eq!(csr_entry(ad, 0, 0), Some(0));
        assert_eq!(csr_entry(ad, 2, 0), Some(4));
        assert_eq!(csr_entry(ad, 1, 0), None);
        assert_eq!(m.x.nnz(), 2);
        assert_eq!(m.missing.len(), 0);
    }

    #[test]
    fn missing_int_is_unstored_and_counted() {
        let mut b = SparseBuilder::new(2, &int_layers());
        b.push_column(
            &[0, 1],
            &[Some(1), Some(2)],
            &[vec![Cell::Int(None), Cell::Int(Some(3))]],
        );
        let m = b.finish();
        let LayerMatrix::I32(ad) = &m.layers[0].1 else {
            panic!("wrong kind")
        };
        assert_eq!(ad.nnz(), 1);
        assert_eq!(m.missing, vec![("AD_alt".to_string(), 1)]);
    }

    #[test]
    fn missing_float_is_nan() {
        let mut b = SparseBuilder::new(1, &[("AF".to_string(), Kind::Float)]);
        b.push_column(&[0], &[Some(1)], &[vec![Cell::Float(None)]]);
        let m = b.finish();
        let LayerMatrix::F32(af) = &m.layers[0].1 else {
            panic!("wrong kind")
        };
        assert!(csr_entry(af, 0, 0).unwrap().is_nan());
        assert_eq!(m.missing, vec![("AF".to_string(), 1)]);
    }

    #[test]
    fn x_only_where_called() {
        let mut b = SparseBuilder::new(2, &int_layers());
        b.push_column(
            &[0, 1],
            &[None, Some(2)],
            &[vec![Cell::Int(Some(1)), Cell::Int(Some(1))]],
        );
        let m = b.finish();
        assert_eq!(m.x.nnz(), 1);
        assert_eq!(csr_entry(&m.x, 1, 0), Some(2));
    }

    #[test]
    fn depth_only_cell_is_covered_with_missing_alt() {
        let ad = [Some(5), None];
        assert!(crate::genotype::covered(false, &ad));
        let layers = vec![
            ("AD_ref".to_string(), Kind::Integer),
            ("AD_alt".to_string(), Kind::Integer),
        ];
        let mut b = SparseBuilder::new(1, &layers);
        b.push_column(
            &[0],
            &[None],
            &[vec![Cell::Int(Some(5))], vec![Cell::Int(None)]],
        );
        let m = b.finish();
        let LayerMatrix::I32(r) = &m.layers[0].1 else {
            panic!("wrong kind")
        };
        assert_eq!(csr_entry(r, 0, 0), Some(5));
        assert_eq!(m.missing, vec![("AD_alt".to_string(), 1)]);
    }

    #[test]
    fn columns_accumulate_across_pushes() {
        let mut b = SparseBuilder::new(2, &int_layers());
        b.push_column(&[1], &[Some(1)], &[vec![Cell::Int(Some(7))]]);
        b.push_column(
            &[0, 1],
            &[Some(3), Some(2)],
            &[vec![Cell::Int(Some(8)), Cell::Int(Some(9))]],
        );
        let m = b.finish();
        assert_eq!((m.x.nrows(), m.x.ncols()), (2, 2));
        assert_eq!(csr_entry(&m.x, 1, 0), Some(1));
        assert_eq!(csr_entry(&m.x, 0, 1), Some(3));
        let LayerMatrix::I32(ad) = &m.layers[0].1 else {
            panic!("wrong kind")
        };
        assert_eq!(csr_entry(ad, 1, 1), Some(9));
    }
}
