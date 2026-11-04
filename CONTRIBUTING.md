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

- `fmt-check` verifies formatting with `cargo fmt -- --check`; `pixi run fmt` applies it.
- `lint` runs clippy on all targets with warnings denied.
- `test` runs the Rust unit and integration tests.
- `build` compiles the release binary.
- `test-readback` builds the release binary first, writes the fixture h5ad files, and reads every output back with Python `anndata`.

## Releases

Push a tag that matches the version in `Cargo.toml` (such as `0.1.0`) from the `main` branch.
The publish workflow tests, builds binaries for linux x86_64 and macOS arm64, and creates the GitHub release with a changelog from `git-cliff`.
