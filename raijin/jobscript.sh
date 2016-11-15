#!/bin/bash
# properties = {properties}

set -euo pipefail

module load adapterremoval sra-toolkit sickle zstd khmer seqhax

export TMPDIR=$PBS_JOBFS

{exec_job}
