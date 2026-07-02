# HIPPO Onboarding Guide

> **H**it **I**nteraction **P**rofiling and **P**rocurement **O**ptimisation for fragment-based drug discovery

## Project Overview

| | |
|---|---|
| **Package** | `xchem-hippo` |
| **Languages** | Python, SQL, YAML, Dockerfile, Shell |
| **Frameworks** | Django (ORM-only), Docker/Compose, pytest, GitHub Actions |
| **Domain** | Fragment-based drug discovery (FBDD) |

HIPPO manages the lifecycle of fragment hits: ingesting crystallographic poses from Fragalysis, computing protein-ligand interactions, and optimising synthesis routes for hit progression. It uses Django purely as an ORM (not as a web framework) backed by PostgreSQL with the RDKit cartridge for chemistry-aware queries.

---

## Architecture Layers

```
┌──────────────────────────────────────────────────┐
│  Orchestration   animal.py, client.py, bootstrap │
├──────────────────────────────────────────────────┤
│  Services        ingestion, download, interaction│
│                  pose, compound, reaction, route  │
│                  recipe, generation, scoring      │
├──────────────────────────────────────────────────┤
│  Domain          Components (Compound, Pose,     │
│                  Reaction, Quote, Recipe)         │
│                  Sets (PoseSet, CompoundSet, etc.)│
├──────────────────────────────────────────────────┤
│  Data            Django ORM models, managers,    │
│                  PostgreSQL schema + RDKit ext    │
├──────────────────────────────────────────────────┤
│  Utility         SMILES/InChI processing, chem   │
│                  validation, Fragalysis parsers   │
├──────────────────────────────────────────────────┤
│  Infrastructure  Docker, CI/CD, config           │
└──────────────────────────────────────────────────┘
```

### Layer Descriptions

1. **Orchestration** — Top-level `HIPPO` class (`Animal`) and client managers that coordinate services and expose the public API for hit ingestion, quoting, and synthesis workflows.

2. **Service** — Business logic: compound lookup, data ingestion, interaction fingerprinting, pose creation, reaction management, recipe scoring, and synthesis route persistence.

3. **Domain** — Rich domain objects (`Compound`, `Pose`, `Reaction`, `Quote`, `Recipe`) wrapping ORM instances with chemistry-aware behavior, plus collection classes (`PoseSet`, `CompoundSet`) for analytics.

4. **Data** — Django ORM models, custom managers, and PostgreSQL schema definitions. The RDKit cartridge enables substructure searching and fingerprint similarity directly in SQL.

5. **Utility** — Shared helpers for SMILES canonicalization, InChI-Key generation, tautomer hashing, chemistry validation, and Fragalysis parsing.

6. **Infrastructure** — Docker containers (including RDKit+PostgreSQL cartridge builds), Compose orchestration, CI/CD workflows, and PyPI release configuration.

---

## Key Concepts

| Concept | Meaning |
|---------|---------|
| **Hit** | A fragment compound confirmed to bind a protein target |
| **Pose** | A 3D orientation of a compound bound to the target (from X-ray crystallography) |
| **Interaction** | A geometric protein-ligand contact (H-bond, hydrophobic, pi-stack) |
| **Reaction** | A chemical transformation with product and reactant compounds |
| **Route** | A sequence of reactions forming a synthesis pathway |
| **Recipe** | A scored combination of routes optimized for cost, novelty, etc. |
| **Fragalysis** | External platform providing crystallographic fragment screening data |
| **RDKit Cartridge** | PostgreSQL extension for molecular structure storage and substructure search |
| **Tautomer Hash** | Canonical compound identity that's insensitive to tautomeric forms |
| **Animal** | The central orchestrator object (a HIPPO is the "animal") |

### Design Patterns

- **Django as ORM-only** — `bootstrap.py` calls `django.setup()` with minimal config; no web server involved.
- **Service layer** — All mutation logic lives in service modules; domain components are behavior-rich but don't write to the DB directly.
- **Set classes** — `CompoundSet`, `PoseSet`, etc. wrap querysets with domain analytics (fingerprints, clustering, scoring).
- **Component wrappers** — `Compound`, `Pose`, `Reaction` wrap model instances and add chemistry methods.

---

## Guided Tour

Follow this path to understand the codebase from entry point to deployment:

### 1. Project README
Read `README.md` for domain vocabulary (hits, poses, interactions, recipes) and installation.

### 2. Package Entry Point
`hippo/__init__.py` exposes `load_hippo` — the single factory function users call.

### 3. Django Bootstrap
`hippo/bootstrap.py` configures Django settings and creates the HIPPO instance. Django is used purely as an ORM here.

### 4. The Animal Orchestrator
`hippo/designdb/animal.py` is the central "god object" wiring all services together. Its surface area reveals HIPPO's full capabilities.

### 5. Data Model
`hippo/designdb/models.py` defines the schema as Django models (16 incoming import edges — the most-imported file). Compounds → Poses → Interactions → Reactions → Routes.

### 6. Database Schema
`images/xchem-designdb/01_schema.sql` defines PostgreSQL tables with RDKit cartridge extensions for molecular storage and substructure search.

### 7. Data Ingestion
`hippo/designdb/services/ingestion.py` + `services/download.py` — the primary data entry point, parsing Fragalysis downloads and SDF files into structured records.

### 8. Interaction Profiling
`hippo/designdb/services/interaction.py` + `hippo/designdb/interactions.py` — detects H-bonds, hydrophobic contacts, and pi-stacking from 3D coordinates. The "IP" in HIPPO.

### 9. Synthesis Planning
`hippo/designdb/services/recipe.py`, `services/route.py`, `services/recipe_score.py` — the "PO" in HIPPO. Generates and scores synthesis route combinations.

### 10. Domain Components
`hippo/designdb/components/` — `Compound` adds SMILES manipulation and RDKit molecule generation; `Pose` adds spatial operations.

### 11. Collections & Analytics
`hippo/designdb/sets/` — `CompoundSet`, `PoseSet`, `InteractionSet` with fingerprint generation, community detection, and visualization.

### 12. Client Managers
`hippo/designdb/client.py` — user-facing API surface (`RecipeManager`, `ScorerManager`, etc.) composing services and sets into high-level operations.

### 13. Chemistry Utilities
`hippo/designdb/utils.py` + `utils_chem.py` — SMILES canonicalization, InChI-Key generation, tautomer hashing.

### 14. Infrastructure
`docker-compose.yaml` + `images/postgres/Dockerfile` — local dev environment with PostgreSQL compiled with the RDKit cartridge.

### 15. Release Pipeline
`.github/workflows/release.yaml` + `pyproject.toml` — builds and publishes `hippo-db` to PyPI on tag.

---

## File Map

### Orchestration

| File | Purpose |
|------|---------|
| `hippo/__init__.py` | Package entry; exposes `load_hippo` factory |
| `hippo/bootstrap.py` | Configures Django, instantiates HIPPO animal |
| `hippo/designdb/animal.py` | Central orchestrator wiring all services |
| `hippo/designdb/client.py` | User-facing manager classes (Recipe, Scorer, Generator) |
| `hippo/ta_auth_connector.py` | HTTP client for Target Access Authenticator |

### Services

| File | Purpose |
|------|---------|
| `services/ingestion.py` | Parses Fragalysis, SDF, Syndirella data into DB |
| `services/download.py` | Downloads target archives from Fragalysis API |
| `services/interaction.py` | Geometric protein-ligand interaction detection |
| `services/compound.py` | Compound lookup via tautomer-insensitive hashes |
| `services/pose.py` | Pose creation with RMSD deduplication |
| `services/reaction.py` | Chemical reaction creation and linking |
| `services/route.py` | Synthesis route persistence |
| `services/recipe.py` | Recipe factory methods (from reactions, compounds) |
| `services/recipe_score.py` | Multi-attribute recipe scoring |
| `services/generation.py` | Random recipe/selection generators |
| `services/quote.py` | Procurement quoting |
| `services/subsite.py` | Binding subsite management |

### Domain Components

| File | Purpose |
|------|---------|
| `components/compound.py` | Compound wrapper: SMILES, scaffolds, descriptors |
| `components/pose.py` | Pose wrapper: ligand molecule, protein, fingerprints |
| `components/reaction.py` | Reaction wrapper: product/reactant relationships |
| `components/price.py` | Price value object with currency handling |
| `components/quote.py` | Quote for procurement pricing |

### Domain Sets

| File | Purpose |
|------|---------|
| `sets/compound.py` | CompoundSet: set operations, filtering, export |
| `sets/pose.py` | PoseSet: scoring, community detection, Fragalysis export |
| `sets/interaction.py` | InteractionSet: fingerprints, h-index, aggregation |
| `sets/ingredient.py` | IngredientSet: DataFrame-backed compound amounts |
| `sets/reaction.py` | ReactionSet: product/reactant composition queries |
| `sets/route.py` | RouteSet: pruning, supplier availability |

### Data Layer

| File | Purpose |
|------|---------|
| `hippo/designdb/models.py` | Django ORM models for all domain entities |
| `hippo/designdb/managers.py` | Custom queryset managers |
| `hippo/designdb/recipe.py` | Recipe/Route aggregates and scoring |
| `images/xchem-designdb/01_schema.sql` | PostgreSQL schema with RDKit extensions |

### Utilities

| File | Purpose |
|------|---------|
| `hippo/designdb/utils.py` | SMILES sanitization, InChI-Key, tautomer hashing |
| `hippo/designdb/utils_chem.py` | Reaction stoichiometry validation |
| `hippo/designdb/utils_frag.py` | Fragalysis longcode parsing, URL resolution |
| `hippo/designdb/utils_xca.py` | XChemAlign longcode decoding |
| `hippo/designdb/interactions.py` | Interaction type constants and cutoffs |
| `hippo/designdb/plotting.py` | Plotly interaction punch-card visualizations |
| `hippo/designdb/settings.py` | Configuration settings |

### Tests

| File | Purpose |
|------|---------|
| `tests/test_01_fragalysis_download.py` | Fragalysis download integration |
| `tests/test_02_setup_animal.py` | Animal instantiation |
| `tests/test_03_add_hits.py` | Hit loading pipeline |
| `tests/test_04_interactions.py` | Interaction calculation |
| `tests/test_05_scaffolds.py` | Scaffold computation |
| `tests/test_compound.py` | Compound property smoke tests |
| `tests/test_pose.py` | Pose property smoke tests |

---

## Complexity Hotspots

These files have the highest complexity and warrant careful study:

| File | Why It's Complex |
|------|-----------------|
| `hippo/designdb/animal.py` | God object with 12+ imports; coordinates all services |
| `hippo/designdb/models.py` | 16 incoming edges; entire schema in one file |
| `hippo/designdb/services/ingestion.py` | Multi-format parser (Fragalysis, SDF, Syndirella) |
| `hippo/designdb/services/interaction.py` | 3D geometric calculations with distance/angle cutoffs |
| `hippo/designdb/sets/pose.py` | Rich querying, scoring, community detection, export |
| `hippo/designdb/components/compound.py` | SMILES manipulation, tautomers, molecular descriptors |
| `hippo/designdb/recipe.py` | Recipe/Route aggregates, scoring, optimization |
| `hippo/designdb/services/generation.py` | Random sampling within budget constraints |
| `hippo/designdb/services/recipe_score.py` | Multi-attribute weighted scoring |
| `images/postgres/Dockerfile` | Multi-stage RDKit C++ compilation for PostgreSQL |

---

## Getting Started

1. **Clone and set up Docker**: `docker compose up` brings up PostgreSQL with RDKit cartridge
2. **Install the package**: `pip install -e .` (or use `uv`)
3. **Run tests**: See `README-test.md` for Docker-based test setup
4. **Try the API**:
   ```python
   from hippo import load_hippo
   h = load_hippo("my_target")
   ```
5. **Follow the tour above** starting from the entry point

---

*Generated from knowledge graph analysis on 2026-07-02.*
