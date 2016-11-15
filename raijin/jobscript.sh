#!/bin/bash
# properties = {properties}

set -euo pipefail

module load adapterremoval sra-toolkit sickle zstd khmer seqhax

{exec_job}
