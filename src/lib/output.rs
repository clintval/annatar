use nalgebra_sparse::CsrMatrix;

pub const X_ENCODING: &str = "number of copies of this ALT in the called genotype, capped at \
    this ALT's read count; hom-ref and uncalled are not stored";

pub enum LayerMatrix {
    I32(CsrMatrix<i32>),
    F32(CsrMatrix<f32>),
}

pub enum UnsValue {
    Text(String),
    List(Vec<String>),
}

pub enum Column {
    Int(Vec<i64>),
    NullableInt(Vec<Option<i64>>),
    Float(Vec<f64>),
    Flag(Vec<bool>),
    NullableFlag(Vec<Option<bool>>),
    Text(Vec<Option<String>>),
}

pub struct Output {
    pub obs_names: Vec<String>,
    pub obs: Vec<(String, Column)>,
    pub var_names: Vec<String>,
    pub var: Vec<(String, Column)>,
    pub x: CsrMatrix<i8>,
    pub layers: Vec<(String, LayerMatrix)>,
    pub uns: Vec<(String, UnsValue)>,
}
