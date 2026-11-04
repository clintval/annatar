import warnings

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import scipy.sparse as sp


def test_fixture_encodings(fixture_path):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        warnings.filterwarnings("ignore", message=".*future.infer_string.*")
        a = ad.read_h5ad(fixture_path)
        assert a.shape == (3, 2) and a.X.dtype == np.int8 and sp.issparse(a.X)
        assert a.layers["AD_alt"].nnz == 2 and a.layers["AD_alt"][0, 0] == 0
        assert np.isnan(a.layers["VAF"].toarray()[0, 0])
        assert isinstance(a.var["SYMBOL"].dtype, pd.CategoricalDtype)
        assert pd.isna(a.var["SYMBOL"].iloc[1])
        assert np.isnan(a.var["REVEL_score"].iloc[1])
        assert str(a.var["AC"].dtype) == "Int64" and pd.isna(a.var["AC"].iloc[1])
        assert a.var["INDEL"].dtype == bool
        assert list(a.uns["x_encoding"])[1].startswith("1:")
        assert a.uns["source_vcf"] == "in.vcf"
    md.MuData({"snv": a})
