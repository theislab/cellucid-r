"""Serve a web build folder the way www.cellucid.com does, for the web cache tests.

    python web_source.py FOLDER PORT_FILE

Each file is answered with the content type its inventory entry declares (the
inventory itself as application/json; charset=utf-8) and an exact
Content-Length. As GitHub Pages does, a client that accepts gzip gets the body
gzipped, with Content-Encoding: gzip and the compressed length. The chosen port is written to PORT_FILE; the server stops when
FOLDER/.stop exists.
"""

import gzip
import json
import os
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

folder, port_file = sys.argv[1:3]
inventory_name = "cellucid-web-assets.json"


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        path = self.path.lstrip("/")
        with open(os.path.join(folder, inventory_name), "rb") as handle:
            types = {asset["path"]: asset["content_type"] for asset in json.load(handle)["assets"]}
        types[inventory_name] = "application/json; charset=utf-8"
        target = os.path.join(folder, path)
        if path not in types or not os.path.isfile(target):
            self.send_response(404)
            self.send_header("Content-Length", "0")
            self.end_headers()
            return
        with open(target, "rb") as handle:
            body = handle.read()
        gzipped = "gzip" in self.headers.get("Accept-Encoding", "")
        if gzipped:
            body = gzip.compress(body, mtime=0)
        self.send_response(200)
        self.send_header("Content-Type", types[path])
        if gzipped:
            self.send_header("Content-Encoding", "gzip")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *args):
        pass


server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)


def watch():
    import time

    while not os.path.exists(os.path.join(folder, ".stop")):
        time.sleep(0.1)
    server.shutdown()


threading.Thread(target=watch, daemon=True).start()
with open(port_file + ".tmp", "w") as handle:
    handle.write(str(server.server_address[1]))
os.replace(port_file + ".tmp", port_file)
server.serve_forever()
