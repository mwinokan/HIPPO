PYTHON = uv run python

SRC = hippo
TESTS = tests

.PHONY: help install lint format typecheck check test ci build clean

help:
	@echo "Available commands:"
	@echo "  make install   Install dependencies"
	@echo "  make lint      Run ruff lint"
	@echo "  make format    Run ruff format"
	@echo "  make typecheck Run mypy"
	@echo "  make check     Run all checks"
	@echo "  make test      Run tests"
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

test:
	uv run pytest

ci: check test

build:
	uv run python -m build

clean:
	rm -rf build dist *.egg-info .pytest_cache .mypy_cache .ruff_cache
