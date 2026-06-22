#!/bin/bash
# Brings up a virtual display so GUI apps (gmsh's viewer) can run inside
# the container with no real screen attached, and exposes it as a web
# page via noVNC -- open http://localhost:6080/vnc.html from any
# browser, on any host OS, to see it.

Xvfb "$DISPLAY" -screen 0 1280x800x24 &
sleep 1

x11vnc -display "$DISPLAY" -nopw -forever -shared -quiet &
sleep 1

websockify --web=/usr/share/novnc 6080 localhost:5900 &

# Shortcut so STEP files can be loaded as /Desktop/foo.step instead of
# /host_home/Desktop/foo.step -- a symlink, not mounting the whole home
# directory at container root, since that would hide system folders
# like /bin and /usr. Ask for more shortcuts (Documents, Downloads...)
# the same way if needed.
ln -sfn /host_home/Desktop /Desktop

# scene.json lives in web/, which is bind-mounted from the host -- so
# it survives a container restart even though everything else in the
# container doesn't. Without this, reopening the viewer after a fresh
# 'docker compose up' would still show whatever model was loaded/solved
# in a previous session instead of starting from a blank state.
rm -f /workspace/web/scene.json

# Serves the simple browser viewer (web/index.html + the mesh.obj /
# labels.json that 'web' regenerates each time) -- open
# http://localhost:8081/ from any browser.
( cd /workspace/web && python3 -m http.server 8081 ) &

# A real interactive shell in the browser (ttyd + xterm.js), embedded
# directly in the web viewer page -- open http://localhost:7681/ (or
# just look at the bottom panel of http://localhost:8081/). Starts by
# running ./openFEM directly, same as typing it in a normal terminal;
# quitting it ('quit' in the CLI, or it not being built yet) drops to a
# plain bash prompt instead of closing the session. --writable is
# needed or ttyd defaults to read-only. No auth: fine for a port only
# ever forwarded to localhost. -t flags match the viewer page's own
# dark palette (see web/index.html's :root colors) so the terminal
# doesn't look like a bolted-on default xterm.js theme.
( cd /workspace && ttyd --writable \
    -t fontFamily='"JetBrains Mono, monospace"' \
    -t fontSize=13 \
    -t 'theme={"background":"#141a22","foreground":"#eef1f6","cursor":"#5b8def","selectionBackground":"rgba(91,141,239,0.3)"}' \
    -p 7681 bash -c 'cd /workspace/build && ./openFEM; exec bash' ) &

exec "$@"
