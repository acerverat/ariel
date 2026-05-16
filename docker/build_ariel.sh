#!/bin/bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"
docker build -f Dockerfile -t ariel-env --no-cache .
