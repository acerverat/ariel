#!/bin/bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"
docker build -f Rascall_Dockerfile -t acerverat/rascall:1.0 .
