#!/usr/bin/env bash
# Serve the repository root on localhost only and print the GWAS browser URL.
# Usage: bash code/00_common/gwas_browser/serve.sh [PORT]
set -euo pipefail
root="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
port="${1:-8765}"
if [ ! -f "$root/output_data/gwas_browser/manifest.json" ]; then
  echo "No prepared data. Run first:" >&2
  echo "  python3 code/00_common/gwas_browser/prepare_tracks.py code/00_common/gwas_browser/tracks.json" >&2
  exit 1
fi
echo "GWAS browser: http://127.0.0.1:${port}/code/00_common/gwas_browser/index.html"
echo "Ctrl-C to stop."
# no-cache: the browser revalidates (cheap 304s) so re-prepared data is never served stale from cache
exec python3 - "$port" "$root" <<'PY'
import functools, http.server, sys

class Handler(http.server.SimpleHTTPRequestHandler):
    protocol_version = "HTTP/1.1"  # keep-alive: the viewer fetches dozens of files at startup

    def end_headers(self):
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def log_message(self, fmt, *args):
        if args and str(args[1]).startswith(("4", "5")) and str(args[1]) != "404":
            super().log_message(fmt, *args)

class Server(http.server.ThreadingHTTPServer):
    request_queue_size = 128  # default 5 drops connections during the startup burst
    daemon_threads = True

port, root = int(sys.argv[1]), sys.argv[2]
server = Server(("127.0.0.1", port), functools.partial(Handler, directory=root))
try:
    server.serve_forever()
except KeyboardInterrupt:
    pass
PY
