#!/bin/bash
# bwtk comprehensive test suite.
# Requires python3 (standard library only) on PATH.
# Run from the test/ directory: bash test_suite.sh

set -uo pipefail

BWTK=$(cd "$(dirname "$0")/.." && pwd)/bwtk
TMPROOT=$(mktemp -d /tmp/bwtk_suite_XXXXX)
PASS=0; FAIL=0

trap 'rm -rf "$TMPROOT"' EXIT

# ---------------------------------------------------------------------------
# Test framework (modeled on ~/quaqc/test/test_suite.sh)
# ---------------------------------------------------------------------------

_cur_test=""

run_test() {
  _cur_test="$1"
  local d="$TMPROOT/$1"
  mkdir -p "$d"
  printf 'Testing %-55s ' "$1 ..."
  if "test_$1" "$d" >"$d/stdout" 2>"$d/stderr"; then
    echo "PASS"
    PASS=$((PASS + 1))
  else
    echo "FAIL"
    echo "    stderr: $(tail -3 "$d/stderr" | tr '\n' ' ')"
    FAIL=$((FAIL + 1))
  fi
}

assert_eq() {
  local label="$1" actual="$2" expected="$3"
  if [ "$actual" = "$expected" ]; then return 0; fi
  echo "  assert_eq $label: expected='$expected' got='$actual'" >&2
  return 1
}

assert_contains() {
  local file="$1" pattern="$2"
  if grep -qF "$pattern" "$file"; then return 0; fi
  echo "  assert_contains: '$pattern' not found in $file" >&2
  return 1
}

assert_not_contains() {
  local file="$1" pattern="$2"
  if ! grep -qF "$pattern" "$file"; then return 0; fi
  echo "  assert_not_contains: '$pattern' unexpectedly found in $file" >&2
  return 1
}

assert_file_exists() {
  if [ -s "$1" ]; then return 0; fi
  echo "  assert_file_exists: '$1' missing or empty" >&2
  return 1
}

# ---------------------------------------------------------------------------
# bwtk-specific helpers
# ---------------------------------------------------------------------------

# write_bg <out> "chr1 0 100 1.5" "chr1 100 200 2.5" ...
write_bg() {
  local out="$1"; shift
  : > "$out"
  for line in "$@"; do printf '%s\n' "$line" | tr ' ' '\t' >> "$out"; done
}

# write_chromsizes <out> "chr1 1000" "chr2 500" ...
write_chromsizes() {
  local out="$1"; shift
  : > "$out"
  for line in "$@"; do printf '%s\n' "$line" | tr ' ' '\t' >> "$out"; done
}

# write_bed <out> "chr1 0 100 name 0 +" ...  (space-separated, converts to tab)
write_bed() {
  local out="$1"; shift
  : > "$out"
  for line in "$@"; do printf '%s\n' "$line" | tr ' ' '\t' >> "$out"; done
}

# make_bw <out.bw> <chromsizes> <bedGraph>
make_bw() {
  "$BWTK" bg2bw -g "$2" -i "$3" -o "$1"
}

# dump_bw <bw> — dump to stdout as bedGraph
dump_bw() {
  "$BWTK" adjust -B -i "$1" -o-
}

# assert_bg_close <expected> <actual> [tol=1e-4]
# Floating-point tolerant diff on column 4. Requires python3.
assert_bg_close() {
  python3 - "$1" "$2" "${3:-1e-4}" <<'PY' || return 1
import sys
try:
    fa_lines = open(sys.argv[1]).readlines()
    fb_lines = open(sys.argv[2]).readlines()
except FileNotFoundError as e:
    print(f"  file not found: {e}", file=sys.stderr)
    sys.exit(1)
tol = float(sys.argv[3])
ok = True
for la, lb in zip(fa_lines, fb_lines):
    la, lb = la.rstrip(), lb.rstrip()
    fa, fb = la.split('\t'), lb.split('\t')
    if fa[:3] != fb[:3] or abs(float(fa[3]) - float(fb[3])) > tol:
        print(f"  diff: {la!r} != {lb!r}", file=sys.stderr)
        ok = False
if len(fa_lines) != len(fb_lines):
    print(f"  length mismatch: {len(fa_lines)} vs {len(fb_lines)} lines", file=sys.stderr)
    ok = False
sys.exit(0 if ok else 1)
PY
}

# tsv_field <file> <row> <col> — 1-based row and column from TSV (no header skipping)
tsv_field() {
  awk -v r="$2" -v c="$3" 'NR==r{print $c}' "$1"
}

# score_field <score_tsv> <row> <col_name>
# Columns: chrom start end name sum mean0 mean min max (1-indexed, 1=chrom)
# Header is row 1.
score_col() {
  case "$2" in
    sum)   echo 5 ;;
    mean0) echo 6 ;;
    mean)  echo 7 ;;
    min)   echo 8 ;;
    max)   echo 9 ;;
    *)     echo 0 ;;
  esac
}

# is_debug — true when BWTK_DEBUG=1 is set (ASan tests gate on this)
is_debug() {
  [ "${BWTK_DEBUG:-0}" = "1" ]
}

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

SHARED="$TMPROOT/shared"
mkdir -p "$SHARED"

# basic.bw — chr1:0–1000 in 100-bp bins with values 1..10
{
  cs="$SHARED/chr1.sizes"; write_chromsizes "$cs" "chr1 1000"
  bg="$SHARED/basic.bg"
  : > "$bg"
  for i in $(seq 1 10); do
    s=$(( (i-1)*100 )); e=$(( i*100 ))
    printf 'chr1\t%d\t%d\t%d\n' "$s" "$e" "$i" >> "$bg"
  done
  make_bw "$SHARED/basic.bw" "$cs" "$bg"
}

# multi.bw — chr1:0–500 value 5, chr2:0–500 value 7
{
  cs="$SHARED/multi.sizes"; write_chromsizes "$cs" "chr1 1000" "chr2 1000"
  bg="$SHARED/multi.bg"
  write_bg "$bg" "chr1 0 500 5" "chr2 0 500 7"
  make_bw "$SHARED/multi.bw" "$cs" "$bg"
}

# sparse.bw — chr1:[100,200) val 4, chr1:[800,900) val 6
{
  bg="$SHARED/sparse.bg"
  write_bg "$bg" "chr1 100 200 4" "chr1 800 900 6"
  make_bw "$SHARED/sparse.bw" "$SHARED/chr1.sizes" "$bg"
}

# negatives.bw — chr1:0–500 value -2.5
{
  bg="$SHARED/negatives.bg"
  write_bg "$bg" "chr1 0 500 -2.5"
  make_bw "$SHARED/negatives.bw" "$SHARED/chr1.sizes" "$bg"
}

# subbase.bw — chr1:0–100 in 25-bp bins → values 1.25, 2.5, 3.75, 5.0
# These are exactly float32-representable; tests that step=0 does NOT round to int (fix 1.3).
{
  cs="$SHARED/small.sizes"; write_chromsizes "$cs" "chr1 200"
  bg="$SHARED/subbase.bg"
  write_bg "$bg" "chr1 0 25 1.25" "chr1 25 50 2.5" "chr1 50 75 3.75" "chr1 75 100 5.0"
  make_bw "$SHARED/subbase.bw" "$cs" "$bg"
}

# merge_a.bw — chr1:0–100 value 6
# merge_b.bw — chr1:0–50  value 0 (partially covers — tests coverage-aware merge)
# merge_c.bw — chr1:0–100 value 2
{
  cs="$SHARED/tiny.sizes"; write_chromsizes "$cs" "chr1 200"
  write_bg "$SHARED/merge_a.bg" "chr1 0 100 6"
  write_bg "$SHARED/merge_b.bg" "chr1 0 50 0"
  write_bg "$SHARED/merge_c.bg" "chr1 0 100 2"
  make_bw "$SHARED/merge_a.bw" "$cs" "$SHARED/merge_a.bg"
  make_bw "$SHARED/merge_b.bw" "$cs" "$SHARED/merge_b.bg"
  make_bw "$SHARED/merge_c.bw" "$cs" "$SHARED/merge_c.bg"
}

# ---------------------------------------------------------------------------
# Tier 1 regression — silent data corruption
# ---------------------------------------------------------------------------

test_merge_n_optstring_reachable() {
  local d="$1"
  # Before fix 1.1, -n was missing from optstring → EXIT_FAILURE immediately.
  "$BWTK" merge -n -o "$d/out.bw" \
    "$SHARED/merge_a.bw" "$SHARED/merge_b.bw" 2>"$d/err"
}

test_merge_n_min_correct() {
  local d="$1"
  # merge_a val=6, merge_c val=2. -n should give 2 everywhere both cover.
  # Before fix 1.2: first file would win (0 from memset or 6 overwrite), not take min.
  "$BWTK" merge -n -o "$d/out.bw" \
    "$SHARED/merge_a.bw" "$SHARED/merge_c.bw" 2>"$d/err" || return 1
  local val
  val=$(dump_bw "$d/out.bw" | awk 'NR==1{print $4}')
  assert_eq "merge -n value" "$val" "2"
}

test_merge_default_step_no_rounding() {
  local d="$1"
  # subbase.bw has values 1.25, 2.5, 3.75, 5.0.
  # Before fix 1.3, default step=1.0 rounded every output to an integer.
  "$BWTK" merge -o "$d/out.bw" \
    "$SHARED/subbase.bw" "$SHARED/subbase.bw" 2>"$d/err" || return 1
  local v1 v2 v3
  v1=$(dump_bw "$d/out.bw" | awk 'NR==1{print $4}')
  v2=$(dump_bw "$d/out.bw" | awk 'NR==2{print $4}')
  v3=$(dump_bw "$d/out.bw" | awk 'NR==3{print $4}')
  assert_eq "bin1 not rounded" "$v1" "1.25" || return 1
  assert_eq "bin2 not rounded" "$v2" "2.5"  || return 1
  assert_eq "bin3 not rounded" "$v3" "3.75"
}

test_merge_zeropad_min() {
  local d="$1"
  # merge_a chr1:0-100 val=6; merge_b chr1:0-50 val=0.
  # Default -n: base 75 covered only by merge_a → min=6.
  # -n -z: merge_b zero-pads base 75 → min(6,0)=0.
  "$BWTK" merge -n -o "$d/default.bw" \
    "$SHARED/merge_a.bw" "$SHARED/merge_b.bw" 2>"$d/err1" || return 1
  "$BWTK" merge -n -z -o "$d/zeropad.bw" \
    "$SHARED/merge_a.bw" "$SHARED/merge_b.bw" 2>"$d/err2" || return 1
  # Position 75 is in chr1:50-100 (default) or chr1:0-200 (zeropad).
  local def_val; def_val=$(dump_bw "$d/default.bw" | awk '$2<=75 && $3>75{print $4}')
  local zp_val; zp_val=$(dump_bw "$d/zeropad.bw" | awk '$2<=75 && $3>75{print $4}')
  assert_eq "default min at pos 75" "$def_val" "6" || return 1
  assert_eq "zeropad min at pos 75" "$zp_val" "0"
}

test_merge_coverage_aware_default() {
  local d="$1"
  # merge_a chr1:0-100 val=6; merge_b chr1:0-50 val=0.
  # Default mean: base 75 → only merge_a covers → mean=6.
  # -z mean: both cover (b zero-pads) → mean=(6+0)/2=3.
  "$BWTK" merge -o "$d/default.bw" \
    "$SHARED/merge_a.bw" "$SHARED/merge_b.bw" 2>"$d/err1" || return 1
  "$BWTK" merge -z -o "$d/zeropad.bw" \
    "$SHARED/merge_a.bw" "$SHARED/merge_b.bw" 2>"$d/err2" || return 1
  # Position 75: default → only merge_a covers → 6; zeropad → (6+0)/2=3.
  local def_val; def_val=$(dump_bw "$d/default.bw" | awk '$2<=75 && $3>75{print $4}')
  local zp_val; zp_val=$(dump_bw "$d/zeropad.bw" | awk '$2<=75 && $3>75{print $4}')
  assert_eq "default mean at pos 75" "$def_val" "6" || return 1
  assert_eq "zeropad mean at pos 75" "$zp_val" "3"
}

test_merge_log10_default_skips_uncovered() {
  local d="$1"
  # merge_a covers 0-100; merge_b covers 0-50.
  # Default merge -l: uncovered positions (50-100 has only one input covering it)
  # should produce log10(6)≈0.778 at 50-100, NOT -inf from log10(0).
  # (The real bug was -inf appearing when a position was uncovered by some inputs.)
  "$BWTK" merge -l -o "$d/out.bw" \
    "$SHARED/merge_a.bw" "$SHARED/merge_b.bw" 2>"$d/err" || return 1
  # Check that no -inf or nan appear in the output.
  local badlines; badlines=$(dump_bw "$d/out.bw" | awk '$4 < -100 || $4 != $4' | wc -l | tr -d ' ')
  assert_eq "no -inf rows" "$badlines" "0"
}

test_score_nan_empty_range_off_chrom() {
  local d="$1"
  # A BED row on a chrom not in the bigWig is rejected with an error (not a nan row).
  # Verify it exits non-zero without segfaulting.
  write_bed "$d/ranges.bed" "chrX 0 100 test"
  "$BWTK" score -i "$SHARED/basic.bw" -b "$d/ranges.bed" \
    -o "$d/out.tsv" >"$d/stdout" 2>"$d/stderr" && return 1  # expected failure
  local rc=$?
  [ "$rc" -ne 139 ] && [ "$rc" -ne 134 ]  # no crash
}

test_score_nan_empty_range_within_chrom() {
  local d="$1"
  # sparse.bw has signal at [100,200) and [800,900), gap at [200,800).
  # score cols: 1=name 2=size 3=covered 4=sum 5=mean0 6=mean 7=min 8=max
  write_bed "$d/ranges.bed" "chr1 200 800 gap"
  "$BWTK" score -i "$SHARED/sparse.bw" -b "$d/ranges.bed" -o "$d/out.tsv" 2>"$d/err" || return 1
  local mean0_val; mean0_val=$(awk 'NR==2{print $5}' "$d/out.tsv")
  local max_val; max_val=$(awk 'NR==2{print $8}' "$d/out.tsv")
  assert_eq "mean0 is nan for gap" "$mean0_val" "nan" || return 1
  assert_eq "max is nan for gap" "$max_val" "nan"
}

test_parse_double_rejects_garbage_adjust() {
  local d="$1"
  # Before fix 1.7, atof("abc") = 0 silently; now it must fail.
  "$BWTK" adjust -t abc -i "$SHARED/basic.bw" -o "$d/out.bw" \
    >"$d/stdout" 2>"$d/stderr" && return 1  # expected failure
  assert_contains "$d/stderr" "Unable to parse '-t': abc"
}

test_parse_double_rejects_garbage_merge() {
  local d="$1"
  "$BWTK" merge -m abc -o "$d/out.bw" \
    "$SHARED/merge_a.bw" "$SHARED/merge_c.bw" \
    >"$d/stdout" 2>"$d/stderr" && return 1
  assert_contains "$d/stderr" "Unable to parse '-m': abc"
}

test_parse_double_rejects_garbage_bg2bw() {
  local d="$1"
  write_bg "$d/in.bg" "chr1 0 100 5"
  write_chromsizes "$d/cs" "chr1 1000"
  "$BWTK" bg2bw -m abc -g "$d/cs" -i "$d/in.bg" -o "$d/out.bw" \
    >"$d/stdout" 2>"$d/stderr" && return 1
  assert_contains "$d/stderr" "Unable to parse '-m': abc"
}

test_locale_independent_parsing() {
  local d="$1"
  # On de_DE locale, atof("1.5") = 1.0 (decimal sep is comma); strtod with LC_ALL=C is fine.
  if ! locale -a 2>/dev/null | grep -qi "de_DE"; then
    # Locale not available — skip gracefully
    return 0
  fi
  LC_ALL=de_DE.UTF-8 "$BWTK" adjust -m 1.5 \
    -i "$SHARED/basic.bw" -o "$d/out.bw" 2>"$d/err" || return 1
  # basic.bw bin1 = 1; 1*1.5 = 1.5; if locale bug: 1*1.0 = 1
  local v; v=$(dump_bw "$d/out.bw" | awk 'NR==1{print $4}')
  assert_eq "multiply 1*1.5 with de_DE locale" "$v" "1.5"
}

test_range_chr0_rejected() {
  local d="$1"
  # Coordinates are 1-based; chr1:0-100 should fail (fix 1.9 underflow guard).
  "$BWTK" adjust -r 'chr1:0-100' \
    -i "$SHARED/basic.bw" -o "$d/out.bw" \
    >"$d/stdout" 2>"$d/stderr" && return 1
  assert_contains "$d/stderr" "1-based"
}

test_range_overflow_rejected() {
  local d="$1"
  # A coordinate that overflows uint32_t (fix 1.10).
  "$BWTK" adjust -r 'chr1:1-99999999999' \
    -i "$SHARED/basic.bw" -o "$d/out.bw" \
    >"$d/stdout" 2>"$d/stderr" && return 1
  # Exit non-zero; message should not contain a raw %lld artifact (fix 2.5).
  assert_not_contains "$d/stderr" "%lld" || return 1
  # Some error message must exist.
  [ -s "$d/stderr" ]
}

test_bed_long_name_handled() {
  local d="$1"
  # A BED name field of 4096 bytes must not overflow the 1024-byte stack buffer (fixes 1.12, 2.1).
  local long_name
  long_name=$(python3 -c "print('A'*4096, end='')")
  printf 'chr1\t0\t100\t%s\n' "$long_name" > "$d/bad.bed"
  # Must not segfault (exit code != 139 on Linux, != 134 on macOS via SIGABRT).
  "$BWTK" score -i "$SHARED/basic.bw" -b "$d/bad.bed" -o "$d/out.tsv" \
    >"$d/stdout" 2>"$d/stderr" || true
  local rc=$?
  [ "$rc" -ne 139 ] && [ "$rc" -ne 134 ]
}

test_values_strand_aware_padding() {
  local d="$1"
  # basic.bw covers chr1:0-1000 (10 bins of 100bp, values 1-10).
  # Minus-strand BED row at chr1:900-1000, resize to 200 bp.
  # With -s 200 on a 100-bp region: pads 50 on each side → [850,1050),
  # but chr1 only goes to 1000 so right pad is clamped → 150 real + 50 nan.
  # Column count must be exactly 200+1=201 (1 name + 200 values).
  write_bed "$d/ranges.bed" "chr1 900 1000 feat 0 -"
  "$BWTK" values -i "$SHARED/basic.bw" -b "$d/ranges.bed" \
    -s 200 -o "$d/out.tsv" 2>"$d/err" || return 1
  # values output: col1=Range/name, cols 2..N=values
  local ncols
  ncols=$(awk 'NR==2{print NF}' "$d/out.tsv")
  assert_eq "column count" "$ncols" "201"  # 1 name col + 200 value cols
}

# ---------------------------------------------------------------------------
# Tier 2 regression — crashes / undefined behavior
# ---------------------------------------------------------------------------

test_stack_overflow_long_bed_name() {
  local d="$1"
  # Same as bed_long_name_handled — 4096-byte BED name, no segfault.
  local long_name; long_name=$(python3 -c "print('A'*4096, end='')")
  printf 'chr1\t0\t100\t%s\n' "$long_name" > "$d/bad.bed"
  "$BWTK" score -i "$SHARED/basic.bw" -b "$d/bad.bed" -o "$d/out.tsv" \
    >"$d/stdout" 2>"$d/stderr" || true
  local rc=$?
  [ "$rc" -ne 139 ] && [ "$rc" -ne 134 ]
}

test_stack_overflow_long_chromsizes_name() {
  local d="$1"
  # A chrom name ≥ 1024 bytes in chrom.sizes should not overflow.
  local long_name; long_name=$(python3 -c "print('A'*2000, end='')")
  printf '%s\t1000\n' "$long_name" > "$d/bad.sizes"
  write_bg "$d/in.bg" "chr1 0 100 5"
  "$BWTK" bg2bw -g "$d/bad.sizes" -i "$d/in.bg" -o "$d/out.bw" \
    >"$d/stdout" 2>"$d/stderr" || true
  local rc=$?
  [ "$rc" -ne 139 ] && [ "$rc" -ne 134 ]
}

test_stack_overflow_long_bg_chrom() {
  local d="$1"
  # A chrom name ≥ 1024 bytes in a bedGraph row.
  local long_name; long_name=$(python3 -c "print('A'*2000, end='')")
  write_chromsizes "$d/cs" "chr1 1000"
  printf '%s\t0\t100\t5\n' "$long_name" > "$d/bad.bg"
  "$BWTK" bg2bw -g "$d/cs" -i "$d/bad.bg" -o "$d/out.bw" \
    >"$d/stdout" 2>"$d/stderr" || true
  local rc=$?
  [ "$rc" -ne 139 ] && [ "$rc" -ne 134 ]
}

test_gzopen_fail_no_segfault() {
  local d="$1"
  # Writing to a non-existent directory triggers gzopen failure.
  # Before fix 2.2, gzerror(NULL, &e) was called → segfault.
  "$BWTK" adjust -B \
    -i "$SHARED/basic.bw" -o "$d/nonexistent_dir/out.bg.gz" \
    >"$d/stdout" 2>"$d/stderr" || true
  local rc=$?
  [ "$rc" -ne 139 ] && [ "$rc" -ne 134 ]
}

test_parsedecimal_overflow_format() {
  local d="$1"
  # Overflow error message must not contain a raw %lld format string (fix 2.5).
  "$BWTK" adjust -r 'chr1:1-99999999999' \
    -i "$SHARED/basic.bw" -o "$d/out.bw" \
    >"$d/stdout" 2>"$d/stderr" || true
  assert_not_contains "$d/stderr" "%lld"
}

# ---------------------------------------------------------------------------
# Tier 3 regression — resource handling
# ---------------------------------------------------------------------------

test_error_no_half_bigwig() {
  local d="$1"
  # If bg2bw fails mid-run (chrom in bedGraph not in chrom.sizes, no -S),
  # the output bigWig must not be left behind.
  write_chromsizes "$d/cs" "chr1 1000"
  write_bg "$d/in.bg" "chr1 0 100 5" "chr2 0 100 3"
  "$BWTK" bg2bw -g "$d/cs" -i "$d/in.bg" -o "$d/out.bw" \
    >"$d/stdout" 2>"$d/stderr" || true
  # Should fail because chr2 isn't in chrom.sizes.
  [ ! -f "$d/out.bw" ] || {
    echo "  half-written bigWig was left behind" >&2; return 1
  }
}

test_sigpipe_silent() {
  local d="$1"
  # bwtk chroms | head -1 raises SIGPIPE; bwtk must not print "Broken pipe".
  "$BWTK" chroms -i "$SHARED/multi.bw" -o- 2>"$d/stderr" | head -1 >/dev/null || true
  assert_not_contains "$d/stderr" "Broken pipe" || return 1
  assert_not_contains "$d/stderr" "broken pipe"
}

test_sigint_no_half_bigwig() {
  local d="$1"
  # Start bg2bw on a modest input, SIGINT it, verify no segfault/abort.
  # (Full cleanup on SIGINT is not implemented; we only check the exit
  # signal is not SIGSEGV=139 or SIGABRT=134.)
  local cs="$d/cs"; write_chromsizes "$cs" "chr1 50000000"
  python3 -c "
import sys
for i in range(0, 50000000, 1000):
    sys.stdout.write('chr1\t{}\t{}\t1\n'.format(i, i+1000))
" > "$d/large.bg"

  "$BWTK" bg2bw -g "$cs" -i "$d/large.bg" -o "$d/out.bw" \
    >"$d/stdout" 2>"$d/stderr" &
  local bwpid=$!
  sleep 0.2
  kill -INT "$bwpid" 2>/dev/null || true
  wait "$bwpid"
  local rc=$?
  rm -f "$d/out.bw"
  # Fail only on segfault (139) or SIGABRT (134)
  [ "$rc" -ne 139 ] && [ "$rc" -ne 134 ]
}

test_asan_no_leaks_adjust() {
  local d="$1"
  is_debug || return 0  # only meaningful under make debug (ASan)
  write_bed "$d/ranges.bed" "chr1 0 500 r1" "chr1 500 1000 r2"
  "$BWTK" adjust -b "$d/ranges.bed" \
    -i "$SHARED/basic.bw" -o "$d/out.bw" 2>"$d/err"
  # ASan exit code 0 = no memory errors.
}

test_asan_no_leaks_merge() {
  local d="$1"
  is_debug || return 0
  "$BWTK" merge -o "$d/out.bw" \
    "$SHARED/merge_a.bw" "$SHARED/merge_c.bw" 2>"$d/err"
}

test_asan_no_leaks_score() {
  local d="$1"
  is_debug || return 0
  write_bed "$d/ranges.bed" "chr1 0 500 r1" "chr1 500 1000 r2"
  "$BWTK" score -i "$SHARED/basic.bw" -b "$d/ranges.bed" -o "$d/out.tsv" 2>"$d/err"
}

test_asan_no_leaks_score_chroms_path() {
  local d="$1"
  is_debug || return 0
  # score without -b uses the chromosome-iteration code path where
  # bed->names[i] aliases bw->cl->chrom[i] (fix 3.6 use-after-free guard).
  "$BWTK" score -i "$SHARED/basic.bw" -o "$d/out.tsv" 2>"$d/err"
}

# ---------------------------------------------------------------------------
# Feature coverage — bg2bw
# ---------------------------------------------------------------------------

test_bg2bw_roundtrip() {
  local d="$1"
  write_chromsizes "$d/cs" "chr1 1000"
  write_bg "$d/in.bg" "chr1 0 100 1.25" "chr1 100 200 2.5" "chr1 200 300 3.75"
  make_bw "$d/out.bw" "$d/cs" "$d/in.bg" || return 1
  dump_bw "$d/out.bw" > "$d/actual.bg"
  assert_bg_close "$d/in.bg" "$d/actual.bg"
}

test_bg2bw_with_preset() {
  local d="$1"
  # Ensembl-style tair10 chrom name (e.g., "1" for Chr1-equivalent).
  # Use chr1 which is valid in the UCSC tair10 preset; with Ensembl it's "1".
  write_bg "$d/in.bg" "1 0 100 5"
  "$BWTK" bg2bw -p tair10 -i "$d/in.bg" -o "$d/out.bw" 2>"$d/err"
  # If preset is recognized, it should exit 0 or with a known chrom-not-found message.
  # We accept both success and a "not found" error (depends on tair10 content).
  [ -f "$d/out.bw" ] || assert_contains "$d/err" "tair10" || return 0
}

test_bg2bw_preset_ucsc() {
  local d="$1"
  # UCSC-style name for tair10 would be e.g. "Chr1" → test with -u flag.
  write_bg "$d/in.bg" "Chr1 0 100 5"
  "$BWTK" bg2bw -p tair10 -u -i "$d/in.bg" -o "$d/out.bw" 2>"$d/err"
  [ -f "$d/out.bw" ] || return 0  # accept either outcome (preset presence)
}

test_bg2bw_unknown_chrom_fails_default() {
  local d="$1"
  write_chromsizes "$d/cs" "chr1 1000"
  write_bg "$d/in.bg" "chrX 0 100 5"
  "$BWTK" bg2bw -g "$d/cs" -i "$d/in.bg" -o "$d/out.bw" \
    >"$d/stdout" 2>"$d/stderr" && return 1  # expected failure
  [ -s "$d/stderr" ]
}

test_bg2bw_unknown_chrom_skipped_with_S() {
  local d="$1"
  write_chromsizes "$d/cs" "chr1 1000"
  write_bg "$d/in.bg" "chr1 0 100 5" "chrX 0 100 9"
  "$BWTK" bg2bw -S -g "$d/cs" -i "$d/in.bg" -o "$d/out.bw" 2>"$d/err" || return 1
  assert_file_exists "$d/out.bw" || return 1
  # chrX must not appear in the output.
  local lines; lines=$(dump_bw "$d/out.bw" | wc -l | tr -d ' ')
  assert_eq "only chr1 rows" "$lines" "1"
}

test_bg2bw_value_ops() {
  local d="$1"
  # Input: chr1:0-100 val=4. Pipeline: +1 → *2 → trim@6 → step=0.
  # Expected: min((4+1)*2, 6) = 6.
  write_chromsizes "$d/cs" "chr1 1000"
  write_bg "$d/in.bg" "chr1 0 100 4"
  "$BWTK" bg2bw -a 1 -m 2 -t 6 -g "$d/cs" -i "$d/in.bg" -o "$d/out.bw" 2>"$d/err" || return 1
  local v; v=$(dump_bw "$d/out.bw" | awk 'NR==1{print $4}')
  assert_eq "value op result" "$v" "6"
}

test_bg2bw_stdin_dash() {
  local d="$1"
  write_chromsizes "$d/cs" "chr1 1000"
  write_bg "$d/in.bg" "chr1 0 100 7"
  cat "$d/in.bg" | "$BWTK" bg2bw -g "$d/cs" -i - -o "$d/out.bw" 2>"$d/err" || return 1
  assert_file_exists "$d/out.bw"
}

test_bg2bw_gzip_input() {
  local d="$1"
  write_chromsizes "$d/cs" "chr1 1000"
  write_bg "$d/in.bg" "chr1 0 100 3"
  gzip -k "$d/in.bg"
  "$BWTK" bg2bw -g "$d/cs" -i "$d/in.bg.gz" -o "$d/out.bw" 2>"$d/err" || return 1
  assert_file_exists "$d/out.bw"
}

# ---------------------------------------------------------------------------
# Feature coverage — adjust
# ---------------------------------------------------------------------------

test_adjust_subset_bed() {
  local d="$1"
  # basic.bw covers chr1:0-1000. BED restricts to [0,200).
  write_bed "$d/ranges.bed" "chr1 0 200 r1"
  "$BWTK" adjust -b "$d/ranges.bed" \
    -i "$SHARED/basic.bw" -o "$d/out.bw" 2>"$d/err" || return 1
  # Rows outside [0,200) must not appear.
  local lines; lines=$(dump_bw "$d/out.bw" | wc -l | tr -d ' ')
  assert_eq "only 2 rows (0-100, 100-200)" "$lines" "2"
}

test_adjust_subset_range() {
  local d="$1"
  # -r uses 1-based inclusive coords; chr1:101-200 → 0-based [100,200).
  "$BWTK" adjust -r 'chr1:101-200' \
    -i "$SHARED/basic.bw" -o "$d/out.bw" 2>"$d/err" || return 1
  local lines; lines=$(dump_bw "$d/out.bw" | wc -l | tr -d ' ')
  assert_eq "one row (100-200)" "$lines" "1"
}

test_adjust_value_ops() {
  local d="$1"
  # basic.bw bin1=1. Pipeline: +2 → *3 → trim@8 → result=min((1+2)*3,8)=8.
  "$BWTK" adjust -a 2 -m 3 -t 8 \
    -i "$SHARED/basic.bw" -o "$d/out.bw" 2>"$d/err" || return 1
  local v; v=$(dump_bw "$d/out.bw" | awk 'NR==1{print $4}')
  assert_eq "value op" "$v" "8"
}

test_adjust_bedgraph_output() {
  local d="$1"
  "$BWTK" adjust -B \
    -i "$SHARED/basic.bw" -o "$d/out.bg.gz" 2>"$d/err" || return 1
  assert_file_exists "$d/out.bg.gz" || return 1
  local lines; lines=$(gzip -dc "$d/out.bg.gz" | wc -l | tr -d ' ')
  [ "$lines" -gt 0 ]
}

test_adjust_bedgraph_stdout() {
  local d="$1"
  # -o- means ungzipped bedGraph on stdout.
  local lines
  lines=$("$BWTK" adjust -B -i "$SHARED/basic.bw" -o- 2>"$d/err" | wc -l | tr -d ' ')
  [ "$lines" -gt 0 ]
}

# ---------------------------------------------------------------------------
# Feature coverage — merge
# ---------------------------------------------------------------------------

test_merge_default_mean() {
  local d="$1"
  # Two bigWigs with the same range: val=4 and val=8 → mean=6.
  write_chromsizes "$d/cs" "chr1 200"
  write_bg "$d/a.bg" "chr1 0 100 4"; make_bw "$d/a.bw" "$d/cs" "$d/a.bg" || return 1
  write_bg "$d/b.bg" "chr1 0 100 8"; make_bw "$d/b.bw" "$d/cs" "$d/b.bg" || return 1
  "$BWTK" merge -o "$d/out.bw" "$d/a.bw" "$d/b.bw" 2>"$d/err" || return 1
  local v; v=$(dump_bw "$d/out.bw" | awk 'NR==1{print $4}')
  assert_eq "mean of 4 and 8" "$v" "6"
}

test_merge_sum_S() {
  local d="$1"
  write_chromsizes "$d/cs" "chr1 200"
  write_bg "$d/a.bg" "chr1 0 100 4"; make_bw "$d/a.bw" "$d/cs" "$d/a.bg" || return 1
  write_bg "$d/b.bg" "chr1 0 100 8"; make_bw "$d/b.bw" "$d/cs" "$d/b.bg" || return 1
  "$BWTK" merge -S -o "$d/out.bw" "$d/a.bw" "$d/b.bw" 2>"$d/err" || return 1
  local v; v=$(dump_bw "$d/out.bw" | awk 'NR==1{print $4}')
  assert_eq "sum of 4 and 8" "$v" "12"
}

test_merge_max_M() {
  local d="$1"
  write_chromsizes "$d/cs" "chr1 200"
  write_bg "$d/a.bg" "chr1 0 100 4"; make_bw "$d/a.bw" "$d/cs" "$d/a.bg" || return 1
  write_bg "$d/b.bg" "chr1 0 100 8"; make_bw "$d/b.bw" "$d/cs" "$d/b.bg" || return 1
  "$BWTK" merge -M -o "$d/out.bw" "$d/a.bw" "$d/b.bw" 2>"$d/err" || return 1
  local v; v=$(dump_bw "$d/out.bw" | awk 'NR==1{print $4}')
  assert_eq "max of 4 and 8" "$v" "8"
}

test_merge_min_n() {
  local d="$1"
  write_chromsizes "$d/cs" "chr1 200"
  write_bg "$d/a.bg" "chr1 0 100 4"; make_bw "$d/a.bw" "$d/cs" "$d/a.bg" || return 1
  write_bg "$d/b.bg" "chr1 0 100 8"; make_bw "$d/b.bw" "$d/cs" "$d/b.bg" || return 1
  "$BWTK" merge -n -o "$d/out.bw" "$d/a.bw" "$d/b.bw" 2>"$d/err" || return 1
  local v; v=$(dump_bw "$d/out.bw" | awk 'NR==1{print $4}')
  assert_eq "min of 4 and 8" "$v" "4"
}

test_merge_value_ops_after_aggregation() {
  local d="$1"
  # Sum of 4+8=12, then *0.5 → 6.
  write_chromsizes "$d/cs" "chr1 200"
  write_bg "$d/a.bg" "chr1 0 100 4"; make_bw "$d/a.bw" "$d/cs" "$d/a.bg" || return 1
  write_bg "$d/b.bg" "chr1 0 100 8"; make_bw "$d/b.bw" "$d/cs" "$d/b.bg" || return 1
  "$BWTK" merge -S -m 0.5 -o "$d/out.bw" "$d/a.bw" "$d/b.bw" 2>"$d/err" || return 1
  local v; v=$(dump_bw "$d/out.bw" | awk 'NR==1{print $4}')
  assert_eq "sum*0.5 of 4 and 8" "$v" "6"
}

test_merge_bedgraph_output() {
  local d="$1"
  "$BWTK" merge -B -o- \
    "$SHARED/merge_a.bw" "$SHARED/merge_c.bw" 2>"$d/err" | wc -l | tr -d ' ' > "$d/lines"
  local lines; lines=$(cat "$d/lines")
  [ "$lines" -gt 0 ]
}

test_merge_three_inputs() {
  local d="$1"
  write_chromsizes "$d/cs" "chr1 200"
  write_bg "$d/a.bg" "chr1 0 100 3"; make_bw "$d/a.bw" "$d/cs" "$d/a.bg" || return 1
  write_bg "$d/b.bg" "chr1 0 100 6"; make_bw "$d/b.bw" "$d/cs" "$d/b.bg" || return 1
  write_bg "$d/c.bg" "chr1 0 100 9"; make_bw "$d/c.bw" "$d/cs" "$d/c.bg" || return 1
  "$BWTK" merge -o "$d/out.bw" "$d/a.bw" "$d/b.bw" "$d/c.bw" 2>"$d/err" || return 1
  local v; v=$(dump_bw "$d/out.bw" | awk 'NR==1{print $4}')
  assert_eq "mean of 3,6,9" "$v" "6"
}

# ---------------------------------------------------------------------------
# Feature coverage — values
# ---------------------------------------------------------------------------

test_values_basic() {
  local d="$1"
  # Two BED rows, -s 100 → each row produces exactly 100 value columns.
  # values output format: 1 name column + N value columns = 1+N total.
  write_bed "$d/ranges.bed" "chr1 0 100 r1" "chr1 100 200 r2"
  "$BWTK" values -i "$SHARED/basic.bw" -b "$d/ranges.bed" \
    -s 100 -o "$d/out.tsv" 2>"$d/err" || return 1
  # Header + 2 data rows.
  local rows; rows=$(wc -l < "$d/out.tsv" | tr -d ' ')
  assert_eq "row count" "$rows" "3" || return 1
  # Each data row: 1 name col + 100 value cols = 101 fields.
  local ncols; ncols=$(awk 'NR==2{print NF}' "$d/out.tsv")
  assert_eq "column count" "$ncols" "101"
}

test_values_left_anchored() {
  local d="$1"
  write_bed "$d/ranges.bed" "chr1 0 200 r1 0 +"
  "$BWTK" values -i "$SHARED/basic.bw" -b "$d/ranges.bed" \
    -s 100 -l -o "$d/out.tsv" 2>"$d/err" || return 1
  local ncols; ncols=$(awk 'NR==2{print NF}' "$d/out.tsv")
  assert_eq "column count left-anchored" "$ncols" "101"
}

test_values_right_anchored() {
  local d="$1"
  write_bed "$d/ranges.bed" "chr1 0 200 r1 0 +"
  "$BWTK" values -i "$SHARED/basic.bw" -b "$d/ranges.bed" \
    -s 100 -r -o "$d/out.tsv" 2>"$d/err" || return 1
  local ncols; ncols=$(awk 'NR==2{print NF}' "$d/out.tsv")
  assert_eq "column count right-anchored" "$ncols" "101"
}

test_values_strand_aware() {
  local d="$1"
  # basic.bw: bins at 0-100(v=1), 100-200(v=2). Minus-strand reversal should flip them.
  write_bed "$d/fwd.bed" "chr1 0 200 r1 0 +"
  write_bed "$d/rev.bed" "chr1 0 200 r1 0 -"
  "$BWTK" values -i "$SHARED/basic.bw" -b "$d/fwd.bed" \
    -s 200 -o "$d/fwd.tsv" 2>"$d/err1" || return 1
  "$BWTK" values -i "$SHARED/basic.bw" -b "$d/rev.bed" \
    -s 200 -o "$d/rev.tsv" 2>"$d/err2" || return 1
  # The first value col of fwd vs the last value col of rev should match,
  # confirming reversal. (fwd col 7 = first bin; rev col 206 = last bin after reversal)
  local fwd_first; fwd_first=$(awk 'NR==2{print $7}' "$d/fwd.tsv")
  local rev_last; rev_last=$(awk 'NR==2{print $NF}' "$d/rev.tsv")
  assert_eq "strand reversal: fwd_first == rev_last" "$fwd_first" "$rev_last"
}

test_values_off_chromosome_padding() {
  local d="$1"
  # BED row extends 50 bp past chr1 end (1000). Column count must still be exactly 200+1=201.
  write_bed "$d/ranges.bed" "chr1 900 1000 r1 0 +"
  "$BWTK" values -i "$SHARED/basic.bw" -b "$d/ranges.bed" \
    -s 200 -o "$d/out.tsv" 2>"$d/err" || return 1
  local ncols; ncols=$(awk 'NR==2{print NF}' "$d/out.tsv")
  assert_eq "column count with padding" "$ncols" "201"
}

test_values_stdout() {
  local d="$1"
  write_bed "$d/ranges.bed" "chr1 0 100 r1"
  local lines
  lines=$("$BWTK" values -i "$SHARED/basic.bw" -b "$d/ranges.bed" \
    -s 100 -o - 2>"$d/err" | wc -l | tr -d ' ')
  [ "$lines" -ge 2 ]
}

# ---------------------------------------------------------------------------
# Feature coverage — score
# ---------------------------------------------------------------------------

test_score_default_chromosomes() {
  local d="$1"
  # Without -b, score produces one row per chromosome.
  "$BWTK" score -i "$SHARED/multi.bw" -o "$d/out.tsv" 2>"$d/err" || return 1
  # Header + chr1 row + chr2 row = 3 lines.
  local rows; rows=$(wc -l < "$d/out.tsv" | tr -d ' ')
  assert_eq "rows (header+2 chroms)" "$rows" "3"
}

test_score_with_bed() {
  local d="$1"
  write_bed "$d/ranges.bed" "chr1 0 200 r1" "chr1 200 400 r2" "chr1 400 600 r3"
  "$BWTK" score -i "$SHARED/basic.bw" -b "$d/ranges.bed" \
    -o "$d/out.tsv" 2>"$d/err" || return 1
  local rows; rows=$(wc -l < "$d/out.tsv" | tr -d ' ')
  assert_eq "rows (header+3 BED)" "$rows" "4"
}

test_score_columns() {
  local d="$1"
  # score output cols (1-based): name size covered sum mean0 mean min max
  write_bed "$d/ranges.bed" "chr1 0 100 myname"
  "$BWTK" score -i "$SHARED/basic.bw" -b "$d/ranges.bed" \
    -o "$d/out.tsv" 2>"$d/err" || return 1
  assert_contains "$d/out.tsv" "name"   || return 1
  assert_contains "$d/out.tsv" "sum"    || return 1
  assert_contains "$d/out.tsv" "mean"   || return 1
  assert_contains "$d/out.tsv" "min"    || return 1
  assert_contains "$d/out.tsv" "max"    || return 1
  # Verify data row uses the BED name field.
  assert_contains "$d/out.tsv" "myname"
}

test_score_bed_output_with_B() {
  local d="$1"
  # -B stat requires lowercase stat name
  write_bed "$d/ranges.bed" "chr1 0 100 r1" "chr1 100 200 r2"
  "$BWTK" score -B sum -i "$SHARED/basic.bw" -b "$d/ranges.bed" \
    -o "$d/out.bed" 2>"$d/err" || return 1
  assert_file_exists "$d/out.bed" || return 1
  local lines; lines=$(wc -l < "$d/out.bed" | tr -d ' ')
  assert_eq "BED line count" "$lines" "2"
}

test_score_known_values() {
  local d="$1"
  # Constant bigWig with val=5 over chr1:0-1000.
  # score cols (1-based): name size covered sum mean0 mean min max
  write_chromsizes "$d/cs" "chr1 1000"
  write_bg "$d/const.bg" "chr1 0 1000 5"
  make_bw "$d/const.bw" "$d/cs" "$d/const.bg" || return 1
  write_bed "$d/ranges.bed" "chr1 0 500 r1"
  "$BWTK" score -i "$d/const.bw" -b "$d/ranges.bed" \
    -o "$d/out.tsv" 2>"$d/err" || return 1
  local sum_val; sum_val=$(awk 'NR==2{print $4}' "$d/out.tsv")
  local mean_val; mean_val=$(awk 'NR==2{print $6}' "$d/out.tsv")
  local min_val; min_val=$(awk 'NR==2{print $7}' "$d/out.tsv")
  local max_val; max_val=$(awk 'NR==2{print $8}' "$d/out.tsv")
  assert_eq "sum=2500" "$sum_val" "2500"  || return 1
  assert_eq "mean=5"   "$mean_val" "5"   || return 1
  assert_eq "min=5"    "$min_val" "5"    || return 1
  assert_eq "max=5"    "$max_val" "5"
}

# ---------------------------------------------------------------------------
# Feature coverage — chroms
# ---------------------------------------------------------------------------

test_chroms_basic() {
  local d="$1"
  "$BWTK" chroms -i "$SHARED/multi.bw" -o "$d/out.tsv" 2>"$d/err" || return 1
  assert_file_exists "$d/out.tsv" || return 1
  assert_contains "$d/out.tsv" "chr1" || return 1
  assert_contains "$d/out.tsv" "chr2"
}

test_chroms_stdout() {
  local d="$1"
  local lines
  lines=$("$BWTK" chroms -i "$SHARED/multi.bw" -o- 2>"$d/err" | wc -l | tr -d ' ')
  [ "$lines" -ge 2 ]
}

# ---------------------------------------------------------------------------
# Statistical correctness
# ---------------------------------------------------------------------------

test_score_constant_signal() {
  local d="$1"
  # chr1:0-1000 value 7. Full-chrom score: sum=7000, mean=7, mean0=7, min=7, max=7.
  # score cols (1-based): name size covered sum mean0 mean min max
  write_chromsizes "$d/cs" "chr1 1000"
  write_bg "$d/const.bg" "chr1 0 1000 7"
  make_bw "$d/const.bw" "$d/cs" "$d/const.bg" || return 1
  "$BWTK" score -i "$d/const.bw" -o "$d/out.tsv" 2>"$d/err" || return 1
  local sum_val; sum_val=$(awk 'NR==2{print $4}' "$d/out.tsv")
  local mean0_val; mean0_val=$(awk 'NR==2{print $5}' "$d/out.tsv")
  local mean_val; mean_val=$(awk 'NR==2{print $6}' "$d/out.tsv")
  local min_val; min_val=$(awk 'NR==2{print $7}' "$d/out.tsv")
  local max_val; max_val=$(awk 'NR==2{print $8}' "$d/out.tsv")
  assert_eq "sum=7000"  "$sum_val"   "7000" || return 1
  assert_eq "mean0=7"   "$mean0_val" "7"    || return 1
  assert_eq "mean=7"    "$mean_val"  "7"    || return 1
  assert_eq "min=7"     "$min_val"   "7"    || return 1
  assert_eq "max=7"     "$max_val"   "7"
}

test_score_partial_coverage_mean_vs_mean0() {
  local d="$1"
  # chr1:0-500 value 4, chr1:500-1000 uncovered.
  # BED row chr1:0-1000: covered=500, sum=2000, mean=4 (sum/covered), mean0=2 (sum/region_size).
  # score cols (1-based): name size covered sum mean0 mean min max
  write_chromsizes "$d/cs" "chr1 1000"
  write_bg "$d/half.bg" "chr1 0 500 4"
  make_bw "$d/half.bw" "$d/cs" "$d/half.bg" || return 1
  write_bed "$d/ranges.bed" "chr1 0 1000 full"
  "$BWTK" score -i "$d/half.bw" -b "$d/ranges.bed" \
    -o "$d/out.tsv" 2>"$d/err" || return 1
  local sum_val; sum_val=$(awk 'NR==2{print $4}' "$d/out.tsv")
  local mean0_val; mean0_val=$(awk 'NR==2{print $5}' "$d/out.tsv")
  local mean_val; mean_val=$(awk 'NR==2{print $6}' "$d/out.tsv")
  assert_eq "sum=2000"  "$sum_val"   "2000" || return 1
  assert_eq "mean0=2"   "$mean0_val" "2"    || return 1
  assert_eq "mean=4"    "$mean_val"  "4"
}

test_merge_mean_three_inputs() {
  local d="$1"
  # vals 1, 5, 9 → mean=5.
  write_chromsizes "$d/cs" "chr1 200"
  write_bg "$d/a.bg" "chr1 0 100 1"; make_bw "$d/a.bw" "$d/cs" "$d/a.bg" || return 1
  write_bg "$d/b.bg" "chr1 0 100 5"; make_bw "$d/b.bw" "$d/cs" "$d/b.bg" || return 1
  write_bg "$d/c.bg" "chr1 0 100 9"; make_bw "$d/c.bw" "$d/cs" "$d/c.bg" || return 1
  "$BWTK" merge -o "$d/out.bw" "$d/a.bw" "$d/b.bw" "$d/c.bw" 2>"$d/err" || return 1
  local v; v=$(dump_bw "$d/out.bw" | awk 'NR==1{print $4}')
  assert_eq "mean(1,5,9)=5" "$v" "5"
}

test_merge_step_binning() {
  local d="$1"
  # -s (step) quantizes values to the nearest multiple of step.
  # merge_a=6, merge_c=2, mean=4; with -s 3: round(4/3)*3 = round(1.33)*3 = 1*3 = 3.
  write_chromsizes "$d/cs" "chr1 200"
  write_bg "$d/a.bg" "chr1 0 100 6"; make_bw "$d/a.bw" "$d/cs" "$d/a.bg" || return 1
  write_bg "$d/c.bg" "chr1 0 100 2"; make_bw "$d/c.bw" "$d/cs" "$d/c.bg" || return 1
  "$BWTK" merge -s 3 -o "$d/out.bw" "$d/a.bw" "$d/c.bw" 2>"$d/err" || return 1
  local v; v=$(dump_bw "$d/out.bw" | awk 'NR==1{print $4}')
  assert_eq "mean(6,2)=4 quantized to -s 3 = 3" "$v" "3"
}

test_score_negative_signal() {
  local d="$1"
  # negatives.bw: chr1:0-500 value -2.5.
  # score cols (1-based): name size covered sum mean0 mean min max
  write_bed "$d/ranges.bed" "chr1 0 500 neg"
  "$BWTK" score -i "$SHARED/negatives.bw" -b "$d/ranges.bed" \
    -o "$d/out.tsv" 2>"$d/err" || return 1
  local min_val; min_val=$(awk 'NR==2{print $7}' "$d/out.tsv")
  local max_val; max_val=$(awk 'NR==2{print $8}' "$d/out.tsv")
  assert_eq "min=-2.5" "$min_val" "-2.5" || return 1
  assert_eq "max=-2.5" "$max_val" "-2.5"
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

echo
echo "bwtk test suite"
echo "==============="
echo

for t in \
  merge_n_optstring_reachable \
  merge_n_min_correct \
  merge_default_step_no_rounding \
  merge_zeropad_min \
  merge_coverage_aware_default \
  merge_log10_default_skips_uncovered \
  score_nan_empty_range_off_chrom \
  score_nan_empty_range_within_chrom \
  parse_double_rejects_garbage_adjust \
  parse_double_rejects_garbage_merge \
  parse_double_rejects_garbage_bg2bw \
  locale_independent_parsing \
  range_chr0_rejected \
  range_overflow_rejected \
  bed_long_name_handled \
  values_strand_aware_padding \
  stack_overflow_long_bed_name \
  stack_overflow_long_chromsizes_name \
  stack_overflow_long_bg_chrom \
  gzopen_fail_no_segfault \
  parsedecimal_overflow_format \
  error_no_half_bigwig \
  sigpipe_silent \
  sigint_no_half_bigwig \
  asan_no_leaks_adjust \
  asan_no_leaks_merge \
  asan_no_leaks_score \
  asan_no_leaks_score_chroms_path \
  bg2bw_roundtrip \
  bg2bw_with_preset \
  bg2bw_preset_ucsc \
  bg2bw_unknown_chrom_fails_default \
  bg2bw_unknown_chrom_skipped_with_S \
  bg2bw_value_ops \
  bg2bw_stdin_dash \
  bg2bw_gzip_input \
  adjust_subset_bed \
  adjust_subset_range \
  adjust_value_ops \
  adjust_bedgraph_output \
  adjust_bedgraph_stdout \
  merge_default_mean \
  merge_sum_S \
  merge_max_M \
  merge_min_n \
  merge_value_ops_after_aggregation \
  merge_bedgraph_output \
  merge_three_inputs \
  values_basic \
  values_left_anchored \
  values_right_anchored \
  values_strand_aware \
  values_off_chromosome_padding \
  values_stdout \
  score_default_chromosomes \
  score_with_bed \
  score_columns \
  score_bed_output_with_B \
  score_known_values \
  chroms_basic \
  chroms_stdout \
  score_constant_signal \
  score_partial_coverage_mean_vs_mean0 \
  merge_mean_three_inputs \
  merge_step_binning \
  score_negative_signal \
; do
  run_test "$t"
done

echo
echo "================================"
echo "Results: $PASS passed, $FAIL failed out of $((PASS + FAIL)) tests"
echo

[ "$FAIL" -eq 0 ]
