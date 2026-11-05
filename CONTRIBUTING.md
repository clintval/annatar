# Contributing

## Toolchain

Use [`pixi`](https://pixi.sh/) for a reproducible environment including a Rust toolchain, Python with `anndata`, `mudata` and `pytest`, and `bcftools`.

## Building

```bash
pixi run build
./target/release/annatar --help
```

The `dist` profile (`cargo build --profile dist`) creates a release optimized, stripped binary for distribution.

## Python Package

The Python package lives in `python/annatar/` and its extension, `annatar._native`, in `src/lib/python.rs`.
The extension is built by [`maturin`](https://www.maturin.rs) with the `python` feature and without the default `cli` feature, which holds the h5ad writer, so the wheel carries no HDF5.
Build it into the pixi environment and run the Python tests with:

```bash
pixi run develop
pixi run test-readback
```

## Mandatory Pre-submission Checks

All must pass before a change is submitted:

```bash
pixi run fmt-check
pixi run lint
pixi run lint-python
pixi run test
pixi run build
pixi run test-readback
```

- The `fmt-check` task verifies formatting with `cargo fmt -- --check`; `pixi run fmt` applies it.
- The `lint` task runs clippy on all targets with warnings denied.
- The `lint-python` task checks the Python formatting, lint and types with `ruff` and `mypy`.
- The `test` task runs the Rust unit and integration tests.
- The `build` task compiles the release binary.
- The `test-readback` task builds the release binary and the extension first, writes the fixture h5ad files, reads every output back with Python `anndata`, and checks `annatar.read_vcf` against the h5ad of the command on every fixture.

## Releases

Push a tag that matches the version in `Cargo.toml` (such as `0.1.0`) from the `main` branch.
The publish workflow tests, builds binaries for linux x86_64 and macOS arm64, builds and checks the sdist and wheels, publishes them to PyPI, and creates the GitHub release with a changelog from `git-cliff`.
The release job needs a repository environment named `github`.
PyPI needs a trusted publisher set up before the first release, for the project `annatar`, the workflow `publish_annatar.yml` and the environment `pypi`, and the repository needs an environment with that name.
The Linux archive is built on `ubuntu-latest` and needs a matching glibc where it runs.
