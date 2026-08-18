#!/usr/bin/env bash
set -euo pipefail

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "ERROR: initialize Git first with: git init" >&2
  exit 2
fi

status=0
max_bytes=$((10 * 1024 * 1024))
mapfile -d '' files < <(
  git ls-files --cached --others --exclude-standard -z
)

for file in "${files[@]}"; do
  [[ -f "$file" ]] || continue
  if [[ "$file" == examples/* ]]; then
    echo "BLOCKED example file: $file" >&2
    status=1
  fi
  size=$(stat -c %s "$file")
  if (( size > max_bytes )); then
    echo "BLOCKED file larger than 10 MiB: $file ($size bytes)" >&2
    status=1
  fi
done

text_files=()
for file in "${files[@]}"; do
  [[ -f "$file" ]] || continue
  if LC_ALL=C grep -Iq . "$file"; then
    text_files+=("$file")
  fi
done

if (("${#text_files[@]}" > 0)); then
  if rg -n '/cluster[0-9]+/|/home/[^/]+/|AKIA[0-9A-Z]{16}|BEGIN (RSA |OPENSSH )?PRIVATE KEY' "${text_files[@]}"; then
    echo "BLOCKED: possible local path, credential, or private key shown above." >&2
    status=1
  fi
fi

if (( status != 0 )); then
  echo "Public-file audit failed." >&2
  exit "$status"
fi

echo "Public-file audit passed: ${#files[@]} files, examples/ excluded."
