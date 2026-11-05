# Contributing

## Toolchain

Use [`pixi`](https://pixi.sh/) for a reproducible environment including a Rust toolchain, Python with `anndata`, `mudata` and `pytest`, and `bcftools`.

## Building

```bash
pixi run build
./target/release/annatar --help
```

The `dist` profile (`cargo build --profile dist`) creates a release optimized, stripped binary for distribution.

## Mandatory Pre-submission Checks

All must pass before a change is submitted:

```bash
pixi run fmt-check
pixi run lint
pixi run test
pixi run build
pixi run test-readback
```

- The `fmt-check` task verifies formatting with `cargo fmt -- --check`; `pixi run fmt` applies it.
- The `lint` task runs clippy on all targets with warnings denied.
- The `test` task runs the Rust unit and integration tests.
- The `build` task compiles the release binary.
- The `test-readback` task builds the release binary first, writes the fixture h5ad files, and reads every output back with Python `anndata`.

## Releases

Push a tag that matches the version in `Cargo.toml` (such as `0.1.0`) from the `main` branch.
The publish workflow tests, builds binaries for linux x86_64 and macOS arm64, and creates the GitHub release with a changelog from `git-cliff`.
The release job needs a repository environment named `github`.
The Linux archive is built on `ubuntu-latest` and needs a matching glibc where it runs.
