#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(git -C "$script_dir" rev-parse --show-toplevel)"

dry_run=()
case "$#" in
    0) ;;
    1)
        if [[ "$1" == "--dry-run" ]]; then
            dry_run=(--dry-run)
        else
            printf 'Usage: bash %s [--dry-run]\n' "$0" >&2
            exit 2
        fi
        ;;
    *)
        printf 'Usage: bash %s [--dry-run]\n' "$0" >&2
        exit 2
        ;;
esac

sources=()
while IFS= read -r source || [[ -n "$source" ]]; do
    source="${source%%#*}"
    source="${source#"${source%%[![:space:]]*}"}"
    source="${source%"${source##*[![:space:]]}"}"
    [[ -z "$source" ]] && continue
    sources+=("$source")
done < "$script_dir/src_files.txt"

for source in "${sources[@]}"; do
    if [[ ! -f "$repo_root/$source" ]]; then
        printf 'Missing DICE mirror source: %s\n' "$repo_root/$source" >&2
        exit 1
    fi
done

for source in "${sources[@]}"; do
    rsync -avhP "${dry_run[@]}" "$repo_root/$source" "$script_dir/"
done
