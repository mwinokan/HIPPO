#!/usr/bin/env bash
# Start a throwaway Postgres with the RDKit cartridge and the DesignDB schema,
# for the Postgres test tier (`pytest -m postgres`).
#
# Needs `postgres`, `initdb`, `pg_ctl` and `psql` on PATH, from a Postgres that
# has the RDKit cartridge installed. The quickest way is a conda env with the
# same versions as images/xchem-designdb/Dockerfile:
#
#   micromamba create -n designdb-test -c conda-forge \
#       postgresql=18 rdkit-postgresql=2025.09.5
#   micromamba activate designdb-test
#   tests/postgres/start_db.sh
#
# (Don't use `micromamba run`: it waits for the server this script leaves
# running, so it never returns.)
#
# Then run the tests with the connection settings this script prints, e.g.
#
#   DB_HOST=localhost POSTGRES_PORT=54329 DB_NAME=designdb_test \
#       DB_USER=hippo DB_PASSWORD=hippo uv run pytest -m postgres
#
# Stop it with `pg_ctl -D <data dir> stop`; the data dir is disposable.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SCHEMA="${ROOT}/images/xchem-designdb/init-db/01_schema.sql"

PGDATA_DIR="${PGDATA_DIR:-${ROOT}/.pgtest}"
PORT="${POSTGRES_PORT:-54329}"
DB_NAME="${DB_NAME:-designdb_test}"
DB_USER="${DB_USER:-hippo}"
DB_PASSWORD="${DB_PASSWORD:-hippo}"

if [ -e "${PGDATA_DIR}" ]; then
    echo "${PGDATA_DIR} already exists; stop the old server and remove it first:" >&2
    echo "  pg_ctl -D ${PGDATA_DIR} stop; rm -rf ${PGDATA_DIR}" >&2
    exit 1
fi

pwfile="$(mktemp)"
trap 'rm -f "${pwfile}"' EXIT
printf '%s\n' "${DB_PASSWORD}" > "${pwfile}"

initdb --pgdata="${PGDATA_DIR}" --username="${DB_USER}" --pwfile="${pwfile}" \
    --auth=scram-sha-256 --encoding=UTF8 --no-instructions > /dev/null

# TCP on localhost only, no Unix socket: the default socket dir
# (/var/run/postgresql) may not be writable, and a socket in a deep data dir can
# exceed the 107-byte path limit
pg_ctl --pgdata="${PGDATA_DIR}" --log="${PGDATA_DIR}/server.log" --wait \
    -o "-p ${PORT} -c listen_addresses=localhost -c unix_socket_directories=''" \
    start > /dev/null

export PGPASSWORD="${DB_PASSWORD}"
export PGOPTIONS="-c client_min_messages=warning"
psql_() {
    psql --host=localhost --port="${PORT}" --username="${DB_USER}" \
        --no-psqlrc --quiet --set=ON_ERROR_STOP=1 "$@"
}

psql_ --dbname=postgres --command="CREATE DATABASE ${DB_NAME}" > /dev/null

# The schema's last statement, CREATE FUNCTION designdb.end_score_values_load(),
# is indented with non-breaking spaces (U+00A0) and fails to parse, so a fresh
# load -- like the deployed DB -- has begin_score_values_load() but not end_.
# The schema isn't ours to change, so mirror production: load without
# ON_ERROR_STOP, tolerate exactly that one failure, and fail on anything else.
errors="$(psql_ --set=ON_ERROR_STOP=0 --dbname="${DB_NAME}" --file="${SCHEMA}" \
    2>&1 > /dev/null | grep 'ERROR:' || true)"
end_fn_exists="$(psql_ --dbname="${DB_NAME}" --tuples-only --no-align --command \
    "SELECT to_regprocedure('designdb.end_score_values_load()') IS NOT NULL")"

if [ -n "${errors}" ]; then
    if [ "$(printf '%s\n' "${errors}" | wc -l)" -eq 1 ] && [ "${end_fn_exists}" = f ]; then
        echo "warning: designdb.end_score_values_load() not created (known schema" \
            "issue, matches production)" >&2
    else
        echo "schema load failed:" >&2
        printf '%s\n' "${errors}" >&2
        exit 1
    fi
fi

echo "DesignDB test database ready:"
echo "  DB_HOST=localhost POSTGRES_PORT=${PORT} DB_NAME=${DB_NAME}" \
    "DB_USER=${DB_USER} DB_PASSWORD=${DB_PASSWORD}"
echo "Stop with: pg_ctl -D ${PGDATA_DIR} stop"
