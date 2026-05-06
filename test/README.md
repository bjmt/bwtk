# bwtk tests

## Running

```bash
# From the repo root:
make test            # full suite (66 tests)
make test-smoke      # single round-trip sanity check

# From this directory:
bash test_suite.sh
bash test.sh

# With ASan (requires make debug first):
make debug && BWTK_DEBUG=1 bash test/test_suite.sh
```

## Structure

- **`test.sh`** — smoke test: `bg2bw` → `adjust -B` → `diff golden/`. Exits 0/1.
- **`test_suite.sh`** — 66 named tests organized in five categories:
  - **Tier 1** — silent data corruption regressions
  - **Tier 2** — crash/UB regressions
  - **Tier 3** — resource handling (partial-file cleanup, signal handling)
  - **Feature** — happy-path coverage for every subcommand option
  - **Stats** — numerical correctness (float-tolerant via `assert_bg_close`)
- **`golden/`** — fixtures for the smoke test only; per-test fixtures are built inline at runtime.

## Adding a test

1. Define `test_<name>() { local d="$1"; ... }` anywhere in `test_suite.sh`.
2. Append `<name> \` to the registry loop near the bottom of the file.
3. Use the framework helpers: `write_bg`, `write_bed`, `write_chromsizes`, `make_bw`, `dump_bw`, `assert_eq`, `assert_contains`, `assert_bg_close`, `assert_file_exists`.
4. Return non-zero (or call `return 1`) to fail; return 0 (or fall off the end) to pass.
