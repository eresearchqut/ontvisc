#!/bin/bash --login
set +euo pipefail
conda activate ontvisc
set -euo pipefail

echo "CMD $@"
exec $@
