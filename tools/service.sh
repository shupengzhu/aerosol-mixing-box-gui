#!/bin/bash
# Keep a shared SCRAM Box server running (no root needed; cron acts as the watchdog).
#
#   tools/service.sh start | stop | restart | status | install-cron | set-code <code>
#
# Environment (defaults in brackets):
#   SCRAMBOX_APP       unpacked release to serve        [directory containing this script/..]
#   SCRAMBOX_USER_DIR  shared data: users, sessions, access_code   [$HOME/.scrambox-shared]
#   SCRAMBOX_PORT      port on 127.0.0.1                 [8765]
#   SCRAMBOX_WORKERS   concurrent jobs                   [8]
# The server listens on 127.0.0.1 only; publish it with an HTTPS reverse proxy, e.g.
#   sudo tailscale funnel --bg 8765
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
APP=${SCRAMBOX_APP:-$(cd "$HERE/.." && pwd)}
export SCRAMBOX_USER_DIR=${SCRAMBOX_USER_DIR:-$HOME/.scrambox-shared}
export SCRAMBOX_WORKERS=${SCRAMBOX_WORKERS:-8}
PORT=${SCRAMBOX_PORT:-8765}
mkdir -p "$SCRAMBOX_USER_DIR"
PIDF=$SCRAMBOX_USER_DIR/server.pid
LOG=$SCRAMBOX_USER_DIR/server.log

running() {
  [ -f "$PIDF" ] || return 1
  local pid; pid=$(cat "$PIDF")
  kill -0 "$pid" 2>/dev/null && grep -q "scrambox" "/proc/$pid/cmdline" 2>/dev/null
}

start() {
  if running; then echo "running (pid $(cat "$PIDF"))"; return 0; fi
  [ -f "$SCRAMBOX_USER_DIR/access_code" ] || [ -n "${SCRAMBOX_ACCESS_CODE:-}" ] || \
    { echo "no access code: run '$0 set-code <code>' first (the server would be open to anyone)"; return 1; }
  if [ -f "$LOG" ] && [ "$(stat -c %s "$LOG")" -gt 50000000 ]; then mv "$LOG" "$LOG.1"; fi
  cd "$APP" || return 1
  setsid nohup python3 -m scrambox serve --host 127.0.0.1 --port "$PORT" >> "$LOG" 2>&1 < /dev/null &
  echo $! > "$PIDF"
  sleep 2
  if running; then echo "started (pid $(cat "$PIDF")) app=$APP data=$SCRAMBOX_USER_DIR port=$PORT"
  else echo "failed to start, see $LOG"; tail -5 "$LOG"; return 1; fi
}

stop() {
  if running; then kill "$(cat "$PIDF")"; sleep 1; echo stopped; else echo "not running"; fi
  rm -f "$PIDF"
}

case "${1:-status}" in
  start) start ;;
  stop) stop ;;
  restart) stop; start ;;
  status) if running; then echo "running (pid $(cat "$PIDF")) app=$APP port=$PORT";
          curl -s -o /dev/null -w "http %{http_code}\n" "http://127.0.0.1:$PORT/api/session"; else echo "not running"; fi ;;
  set-code)
    [ -n "${2:-}" ] || { echo "usage: $0 set-code <code>"; exit 1; }
    umask 077; printf '%s\n' "$2" > "$SCRAMBOX_USER_DIR/access_code"; echo "access code stored in $SCRAMBOX_USER_DIR/access_code" ;;
  install-cron)
    line="SCRAMBOX_APP=$APP SCRAMBOX_USER_DIR=$SCRAMBOX_USER_DIR SCRAMBOX_PORT=$PORT $HERE/service.sh start >/dev/null 2>&1"
    ( crontab -l 2>/dev/null | grep -v "service.sh start" ; echo "@reboot $line"; echo "*/2 * * * * $line" ) | crontab -
    echo "cron watchdog installed:"; crontab -l | grep service.sh ;;
  *) echo "usage: $0 start|stop|restart|status|install-cron|set-code <code>"; exit 1 ;;
esac
