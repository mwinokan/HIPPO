PYTHON = uv run python

SRC = hippo
TESTS = tests

.PHONY: help install lint format typecheck check test coverage test-postgres ci build clean

help:
	@echo "Available commands:"
	@echo "  make install   Install dependencies"
	@echo "  make lint      Run ruff lint"
	@echo "  make format    Run ruff format"
	@echo "  make typecheck Run mypy"
	@echo "  make check     Run all checks"
	@echo "  make test      Run tests"
	@echo "  make coverage  Run tests with a coverage report"
	@echo "  make test-postgres  Run tests against Postgres (see tests/postgres/start_db.sh)"
	@echo "  make ci        Simulate CI run"
	@echo "  make build     Build package"
	@echo "  make clean     Remove build artifacts"

install:
	uv sync --frozen

lint:
	uv run pre-commit run ruff --all-files

format:
	uv run pre-commit run ruff-format --all-files

typecheck:
	uv run pre-commit run mypy --all-files

check:
	uv run pre-commit run --all-files

# Only the container-free SQLite tier for now; the pre-refactor tests are not
# yet ported (see `sqlite` marker in pyproject.toml).
test:
	uv run pytest -m sqlite

coverage:
	uv run pytest -m sqlite --cov --cov-report=term --cov-report=xml

# The SQLite tests plus the Postgres-only ones, against a DesignDB Postgres
# started by tests/postgres/start_db.sh (connection from the DB_* env vars).
test-postgres:
	HIPPO_TEST_DB=postgres uv run pytest -m "sqlite or postgres" $(PYTEST_ARGS)

ci: check test

build:
	uv run python -m build

clean:
	rm -rf build dist *.egg-info .pytest_cache .mypy_cache .ruff_cache .coverage coverage.xml
