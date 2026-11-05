# Contributing

## Toolchain

Use [`pixi`](https://pixi.sh/) for a reproducible environment including a Rust toolchain, Python with `anndata`, `mudata` and `pytest`, and `bcftools`.

## Building

The Rust engine lives in `src/lib/`, and the Python package and the `annatar` command live in `python/annatar/`.
The extension, `annatar._native`, is built by [`maturin`](https://www.maturin.rs) with the `python` feature from `src/lib/python.rs`.
Build it into the pixi environment and try the command with:

```bash
pixi run develop
pixi run annatar --help
```

Install from GitHub with Rust:

```bash
pip install 'git+https://github.com/clintval/annatar'
```

## Mandatory Pre-submission Checks

All must pass before a change is submitted:

```bash
pixi run fmt-check
pixi run lint
pixi run lint-python
pixi run test
pixi run test-readback
```

- The `fmt-check` task verifies formatting with `cargo fmt -- --check`; `pixi run fmt` applies it.
- The `lint` task runs clippy on all targets with warnings denied.
- The `lint-python` task checks the Python formatting, lint and types with `ruff` and `mypy`.
- The `test` task runs the Rust unit and integration tests.
- The `test-readback` task builds the extension first, drives the real `annatar` command, reads every output back with Python `anndata`, and checks `annatar.read_vcf` against the h5ad of the command on every fixture and `annatar.concat` against `annatar merge`.

## Releases

Push a tag that matches the version in `Cargo.toml` (such as `0.1.0`) from the `main` branch.
The publish workflow tests, builds and checks the sdist and wheels, publishes them to PyPI, and creates the GitHub release with a changelog from `git-cliff`.
The release job needs a repository environment named `github`.
PyPI needs a trusted publisher set up before the first release, for the project `annatar`, the workflow `publish_annatar.yml` and the environment `pypi`, and the repository needs an environment with that name.
