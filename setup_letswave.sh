#!/usr/bin/env bash
set -euo pipefail
LW_ROOT="${MOAZ_LW_ROOT:-/home/vant7ee/scratch/FacerecognitionData/letswave6_moaz_pin}"
COMMIT="7932ccfd570dc55ba1c5ae8e062147afce615496"
if [[ ! -e "$LW_ROOT" ]]; then
  git clone https://github.com/NOCIONS/letswave6.git "$LW_ROOT"
fi
test -d "$LW_ROOT/.git"
if [[ -n "$(git -C "$LW_ROOT" status --porcelain)" ]]; then
  echo "Letswave directory has changes; refusing to alter it: $LW_ROOT" >&2
  exit 1
fi
git -C "$LW_ROOT" checkout --detach "$COMMIT"
echo "Pinned native Letswave: $LW_ROOT"
