#!/usr/bin/env bash
# Re-checks every theorem in check_axioms.lean, and everything it depends on, with kolla.
#
#   LEAN4EXPORT=/path/to/lean4export KOLLA_JAR=/path/to/kolla.jar ./kolla-check.sh
#
# lean4export has to be built for the same Lean release as lean-toolchain.
set -euo pipefail

cd "$(dirname "$0")"

: "${LEAN4EXPORT:?set LEAN4EXPORT to a lean4export binary built for $(cat lean-toolchain)}"
: "${KOLLA_JAR:?set KOLLA_JAR to the kolla CLI assembly jar}"

profile="$(sed 's/.*://' lean-toolchain)"
work="$(mktemp -d -t pctc-export)"
trap 'rm -rf "$work"' EXIT
export_file="$work/pctc.ndjson"

theorems=()
while IFS= read -r theorem; do theorems+=("$theorem"); done \
  < <(sed -n 's/^#print axioms //p' check_axioms.lean)

lake build
lake env "$LEAN4EXPORT" Pctc -- "${theorems[@]}" > "$export_file"

java -Xmx12g -jar "$KOLLA_JAR" --profile "$profile" -q "$export_file"

axiom_args=()
for theorem in "${theorems[@]}"; do axiom_args+=(--axioms "$theorem"); done
java -Xmx12g -jar "$KOLLA_JAR" --profile "$profile" "${axiom_args[@]}" "$export_file"
