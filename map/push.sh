#!/usr/bin/env bash
# push.sh — commit + push the tags pipeline to GitHub, logging every reupload.
#
# Usage:
#   ./push.sh "commit message"      # message optional; defaults to "update pipeline"
#
# Records each push in config/reuploads.csv (timestamp, commit, branch, message),
# so the upload history stays auditable. The commit SHA logged here is the same one
# the launcher stamps into every output run dir (pipeline_version.json), tying any
# given result back to the exact pipeline version that produced it.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"

MSG="${1:-update pipeline}"
LOG="config/reuploads.csv"
BRANCH="$(git rev-parse --abbrev-ref HEAD)"

# Stage + commit code changes, if there are any.
git add -A
if ! git diff --cached --quiet; then
    git commit -m "$MSG"
fi

git push -u origin "$BRANCH"
SHA="$(git rev-parse HEAD)"
TS="$(date '+%Y-%m-%d %H:%M:%S')"

# Append to the timestamped reupload log, then commit the log entry on its own.
[ -f "$LOG" ] || echo "timestamp,commit,branch,message" > "$LOG"
printf '%s,%s,%s,%s\n' "$TS" "$SHA" "$BRANCH" "${MSG//,/;}" >> "$LOG"
git add "$LOG"
git commit -m "log: reupload ${SHA:0:12}"
git push origin "$BRANCH"

echo "pushed $SHA ($BRANCH) — logged in $LOG"
