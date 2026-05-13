#!/usr/bin/env bash
set -euo pipefail

python_bin="${1:-.envs/ldsc-neale/bin/python}"

"$python_bin" - <<'PY'
import bitarray
import numpy
import pandas
import scipy

print("python_ok")
print("numpy", numpy.__version__)
print("scipy", scipy.__version__)
print("pandas", pandas.__version__)
print("bitarray", bitarray.__version__)
PY
