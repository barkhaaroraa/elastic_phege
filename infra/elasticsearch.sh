#!/usr/bin/env bash
# Single-node Elasticsearch for local backbone development, via rootless podman.
#
# Security is disabled and the port is bound to loopback only -- this is a local
# dev cluster, never a deployment. Heap is capped at 1 GB because the host has
# ~7.5 GB total and ES will happily take more than it needs.
set -euo pipefail

ES_VERSION="${ES_VERSION:-9.5.1}"
CONTAINER="${CONTAINER:-phageforge-es}"
PORT="${PORT:-9200}"
HEAP="${HEAP:-1g}"
IMAGE="docker.elastic.co/elasticsearch/elasticsearch:${ES_VERSION}"

cmd="${1:-up}"

case "$cmd" in
  up)
    if podman container exists "$CONTAINER" 2>/dev/null; then
      if [ "$(podman inspect -f '{{.State.Running}}' "$CONTAINER")" = "true" ]; then
        echo "already running: $CONTAINER"
      else
        echo "starting existing container: $CONTAINER"
        podman start "$CONTAINER" >/dev/null
      fi
    else
      echo "creating container: $CONTAINER (es $ES_VERSION, heap $HEAP)"
      podman run -d \
        --name "$CONTAINER" \
        -p "127.0.0.1:${PORT}:9200" \
        -e "discovery.type=single-node" \
        -e "xpack.security.enabled=false" \
        -e "xpack.license.self_generated.type=trial" \
        -e "ES_JAVA_OPTS=-Xms${HEAP} -Xmx${HEAP}" \
        -e "bootstrap.memory_lock=false" \
        -v phageforge-esdata:/usr/share/elasticsearch/data \
        "$IMAGE" >/dev/null
    fi

    printf 'waiting for cluster'
    for _ in $(seq 1 90); do
      status=$(curl -fsS "http://localhost:${PORT}/_cluster/health" 2>/dev/null \
               | sed -n 's/.*"status":"\([a-z]*\)".*/\1/p') || status=""
      if [ "$status" = "green" ] || [ "$status" = "yellow" ]; then
        echo ""
        echo "cluster is ${status} on http://localhost:${PORT}"
        exit 0
      fi
      printf '.'
      sleep 2
    done
    echo ""
    echo "ERROR: cluster did not become healthy in 180s" >&2
    echo "--- last 40 log lines ---" >&2
    podman logs --tail 40 "$CONTAINER" >&2 || true
    exit 1
    ;;

  down)
    podman stop "$CONTAINER" >/dev/null 2>&1 && echo "stopped $CONTAINER" || echo "not running"
    ;;

  destroy)
    podman rm -f "$CONTAINER" >/dev/null 2>&1 || true
    podman volume rm phageforge-esdata >/dev/null 2>&1 || true
    echo "removed container and data volume"
    ;;

  logs)
    podman logs -f "$CONTAINER"
    ;;

  status)
    curl -fsS "http://localhost:${PORT}/_cluster/health?pretty" || echo "unreachable"
    ;;

  *)
    echo "usage: $0 {up|down|destroy|logs|status}" >&2
    exit 2
    ;;
esac
