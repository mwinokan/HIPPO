-- =========================================================
-- designdb Database Schema
-- =========================================================

-- =========================================================
DROP SCHEMA IF EXISTS designdb CASCADE;
-- =========================================================

-- =========================================================
-- PREREQUISITES & EXTENSIONS
-- =========================================================

CREATE SCHEMA IF NOT EXISTS designdb;
CREATE SCHEMA IF NOT EXISTS rdkit;

CREATE EXTENSION IF NOT EXISTS rdkit WITH SCHEMA rdkit;
CREATE EXTENSION IF NOT EXISTS pgcrypto WITH SCHEMA designdb;

SET search_path TO designdb, rdkit, public;

REVOKE CREATE ON SCHEMA public FROM PUBLIC;

-- =========================================================
-- TABLES (ordered by FK dependencies)
-- =========================================================

CREATE TABLE IF NOT EXISTS designdb.targets (
    id BIGSERIAL PRIMARY KEY, --Internal ID inserted when registering target via Fragalysis
    external_target_id BIGINT, -- ID of this target in the external database (Scarab link)
    target_name TEXT NOT NULL, --Insert from HIPPO codebase. Must be a link to Scarab protein production target
    target_metadata TEXT, -- Not populated by code
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_target UNIQUE (target_name)
);

CREATE TABLE IF NOT EXISTS designdb.compounds (
    id BIGSERIAL PRIMARY KEY,
    compound_inchikey TEXT, -- Populated by RDKit cartridge trigger from compound_smiles (do not insert by code). Maybe insert by the codebase or jupyter. Do we need to write this by code or can be calculated by the cartridge?
    compound_alias TEXT, -- Maybe insert by the codebase.
    compound_smiles TEXT, -- Inserted by the codebase. Trigger populates compound_mol and compound_inchikey. 2D flat SMILES. Is this 2d flat smiles without any stereochemistry? LR - Yes 2D. Looks like designdb function sanitise_smiles does remove stereochemistry - will this be a problem when a user wants to register a design with defined stereochemistry?
    compound_hash TEXT NOT NULL, -- Canonical identity hash from application pipeline (e.g. RDKit RegistrationHash after SuperParent); pair with rdkit_version for reproducibility
    base_compound_id BIGINT REFERENCES designdb.compounds (id) ON DELETE SET NULL, -- Not populated by code
    -- compound_mol rdkit.mol, -- Replaced by TEXT CTAB below: JDBC showed SMILES text; mol_to_ctab gives a molfile string that Scarab will easily convert to structure.
    compound_mol TEXT, -- V2000 CTAB (mol block) from rdkit.mol_to_ctab(mol_from_smiles(...))
    compound_pattern_bfp bit(2048), -- Postgres RDkit cartridge can calc this, Chemicalite does, not sure if insert by codebase. Currently seems broken
    compound_morgan_bfp bit(2048), -- Postgres cartridge can't calc this. Must be inserted by codebase, but currently its broken
    compound_metadata TEXT, -- currently Null
    note TEXT,  -- New column
    rdkit_version TEXT, -- RDKit version string used when computing compound_hash (application-set; cartridge may also populate)
    inchi_version TEXT NOT NULL, -- InChI software version (rdkit.Chem.inchi.GetInchiVersion)
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_compound_alias UNIQUE (compound_alias),
    CONSTRAINT uc_compound_inchikey UNIQUE (compound_inchikey),
    CONSTRAINT uc_compound_smiles UNIQUE (compound_smiles)
);

CREATE TABLE IF NOT EXISTS designdb.subsites (
    id BIGSERIAL PRIMARY KEY,
    target_id BIGINT NOT NULL REFERENCES designdb.targets (id) ON DELETE RESTRICT, -- Insert by codebase/notebook
    subsite_name TEXT NOT NULL, -- Insert by codebase/notebook
    subsite_metadata TEXT, -- Not populated by code
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_subsite UNIQUE (target_id, subsite_name)
);

CREATE TABLE IF NOT EXISTS designdb.poses (
    id BIGSERIAL PRIMARY KEY,
    pose_inchikey TEXT, -- Populated by RDKit cartridge trigger from pose_mol (do not insert by code). Originally, nserted by codebase when registering poses? Might be done from pose.mol?
    pose_alias TEXT,
    pose_smiles TEXT, -- Populated by RDKit cartridge trigger from pose_mol (do not insert by code). LR - necessary because will contain defined stereochemistry - should these be canonicalised? Is it done by codebase from pose.mol? Could be done by RDkit cartridge.
    pose_reference INTEGER,
    pose_path TEXT,
    compound_id BIGINT NOT NULL REFERENCES designdb.compounds (id) ON DELETE RESTRICT,
    target_id BIGINT NOT NULL REFERENCES designdb.targets (id) ON DELETE RESTRICT,
    pose_mol rdkit.mol, -- Insert by codebase. Trigger populates pose_inchikey and pose_smiles via cartridge.
    pose_fingerprint INTEGER, --Not sure if it null or actually calcualated somewhere.
    --pose_energy_score REAL, -- LR - redundant; use designdb.score_values
    --pose_distance_score REAL, -- LR - redundant; use designdb.score_values
    --pose_inspiration_score REAL, -- LR - redundant; use designdb.score_values
    pose_metadata TEXT,
    note TEXT,  -- New column
    rdkit_version TEXT, --Can be done by RDkit cartridge
    inchi_version TEXT,  -- Must be done by codebase
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now()
    -- CONSTRAINT uc_pose_alias UNIQUE (pose_alias), -- Removed
    -- CONSTRAINT uc_pose_path UNIQUE (pose_path) -- Removed
);

CREATE TABLE IF NOT EXISTS designdb.subsite_tags (
    id BIGSERIAL PRIMARY KEY,
    subsite_id BIGINT NOT NULL REFERENCES designdb.subsites (id) ON DELETE RESTRICT, -- Insert by codebase
    pose_id BIGINT NOT NULL REFERENCES designdb.poses (id) ON DELETE RESTRICT, -- Insert by codebase
    subsite_tag_metadata TEXT, -- Not populated by code
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_subsite_tag UNIQUE (subsite_id, pose_id)
);

-- New table
CREATE TABLE IF NOT EXISTS designdb.pose_methods (
    id BIGSERIAL PRIMARY KEY,
    pose_method_name TEXT,
    pose_method_description TEXT,
    pose_method_version TEXT,
    pose_method_organization TEXT,
    pose_method_link TEXT,
    pose_method_note TEXT,
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_pose_method UNIQUE NULLS NOT DISTINCT (pose_method_name, pose_method_version)
);

-- New table
CREATE TABLE IF NOT EXISTS designdb.has_pose_methods (
    pose_id BIGINT NOT NULL REFERENCES designdb.poses (id) ON DELETE CASCADE,
    pose_method_id BIGINT NOT NULL REFERENCES designdb.pose_methods (id) ON DELETE CASCADE,
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    PRIMARY KEY (pose_id, pose_method_id)
);

-- New table
CREATE TABLE IF NOT EXISTS designdb.pose_tags (
    id BIGSERIAL PRIMARY KEY,
    pose_tag_name TEXT NOT NULL,
    pose_tag_description TEXT,
    pose_tag_note TEXT,
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_pose_tag_name UNIQUE (pose_tag_name)
);

-- New table
CREATE TABLE IF NOT EXISTS designdb.has_pose_tags (
    pose_id BIGINT NOT NULL REFERENCES designdb.poses (id) ON DELETE CASCADE,
    pose_tag_id BIGINT NOT NULL REFERENCES designdb.pose_tags (id) ON DELETE CASCADE,
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    PRIMARY KEY (pose_id, pose_tag_id)
);

CREATE TABLE IF NOT EXISTS designdb.inspirations (
    id BIGSERIAL PRIMARY KEY,
    original_pose_id BIGINT REFERENCES designdb.poses (id) ON DELETE SET NULL, -- Insert by codebase
    derivative_pose_id BIGINT REFERENCES designdb.poses (id) ON DELETE SET NULL, -- Insert by codebase
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_inspiration UNIQUE (original_pose_id, derivative_pose_id)
);

CREATE TABLE IF NOT EXISTS designdb.features (
    id BIGSERIAL PRIMARY KEY,
    feature_family TEXT, -- Insert by codebase
    target_id BIGINT NOT NULL REFERENCES designdb.targets (id) ON DELETE RESTRICT, -- Insert by codebase
    feature_chain_name TEXT, -- Insert by codebase
    feature_residue_name TEXT, -- Insert by codebase
    feature_residue_number INTEGER, -- Insert by codebase
    feature_atom_name TEXT, -- Insert by codebase
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_feature UNIQUE (feature_family, target_id, feature_chain_name, feature_residue_number, feature_residue_name, feature_atom_name)
);

CREATE TABLE IF NOT EXISTS designdb.interactions (
    id BIGSERIAL PRIMARY KEY,
    feature_id BIGINT NOT NULL REFERENCES designdb.features (id) ON DELETE RESTRICT, -- Insert by codebase
    pose_id BIGINT NOT NULL REFERENCES designdb.poses (id) ON DELETE RESTRICT, -- Insert by codebase
    interaction_type TEXT NOT NULL, -- Insert by codebase
    interaction_family TEXT NOT NULL, -- Insert by codebase
    interaction_atom_id TEXT NOT NULL, -- Insert by codebase
    interaction_prot_coord TEXT NOT NULL, -- Insert by codebase. Not populated by ProLIF, needs reviewing
    interaction_lig_coord TEXT NOT NULL, -- Insert by codebase. Not populated by ProLIF, needs reviewing
    interaction_distance REAL NOT NULL, -- Insert by codebase
    interaction_angle REAL, -- Insert by codebase
    interaction_energy REAL, -- Insert by codebase. Not populated by ProLIF, needs reviewing
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_interaction UNIQUE (feature_id, pose_id, interaction_type, interaction_family, interaction_atom_id)
);

-- New table
CREATE TABLE IF NOT EXISTS designdb.compound_tags (
    id BIGSERIAL PRIMARY KEY,
    compound_tag_name TEXT NOT NULL,
    compound_tag_description TEXT,
    compound_tag_note TEXT,
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_compound_tag_name UNIQUE (compound_tag_name)
);

-- New table
CREATE TABLE IF NOT EXISTS designdb.has_compound_tags (
    compound_id BIGINT NOT NULL REFERENCES designdb.compounds (id) ON DELETE CASCADE,
    compound_tag_id BIGINT NOT NULL REFERENCES designdb.compound_tags (id) ON DELETE CASCADE,
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    PRIMARY KEY (compound_id, compound_tag_id)
);

-- New table
CREATE TABLE IF NOT EXISTS designdb.enumeration_methods (
    id BIGSERIAL PRIMARY KEY,
    enum_name TEXT,
    enum_description TEXT,
    enum_version TEXT,
    enum_organization TEXT,
    enum_link TEXT,
    enum_note TEXT,
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_enumeration_method UNIQUE NULLS NOT DISTINCT (enum_name, enum_version)
);

-- New table
CREATE TABLE IF NOT EXISTS designdb.has_enumeration_methods (
    compound_id BIGINT NOT NULL REFERENCES designdb.compounds (id) ON DELETE CASCADE,
    enumeration_method_id BIGINT NOT NULL REFERENCES designdb.enumeration_methods (id) ON DELETE CASCADE,
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    PRIMARY KEY (compound_id, enumeration_method_id)
);

-- New table
CREATE TABLE IF NOT EXISTS designdb.scoring_methods (
    id BIGSERIAL PRIMARY KEY,
    method_name TEXT,
    method_description TEXT,
    method_version TEXT,
    method_organization TEXT,
    method_link TEXT,
    note TEXT,
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_scoring_method UNIQUE NULLS NOT DISTINCT (method_name, method_version)
);

-- One row per (pose_id, compound_id, scoring_method_id). score JSONB stores {"score": value} (numeric or text).
CREATE TABLE IF NOT EXISTS designdb.score_values (
    pose_id BIGINT NOT NULL REFERENCES designdb.poses (id) ON DELETE RESTRICT,
    compound_id BIGINT NOT NULL REFERENCES designdb.compounds (id) ON DELETE RESTRICT,
    scoring_method_id BIGINT NOT NULL REFERENCES designdb.scoring_methods (id) ON DELETE RESTRICT,
    score JSONB NOT NULL,
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT pk_score_values PRIMARY KEY (pose_id, compound_id, scoring_method_id)
);

CREATE TABLE IF NOT EXISTS designdb.reactions (
    id BIGSERIAL PRIMARY KEY,
    reaction_type TEXT, -- Insert by codebase/notebook, Synderilla
    product_compound_id BIGINT NOT NULL REFERENCES designdb.compounds (id) ON DELETE RESTRICT, -- Insert by codebase/notebook, Synderilla
    reaction_product_yield REAL, -- Insert by codebase/notebook, Synderilla
    reaction_metadata TEXT, -- Not populated by code
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now()
);

CREATE TABLE IF NOT EXISTS designdb.reactants (
    id BIGSERIAL PRIMARY KEY,
    reactant_amount REAL, -- Insert by codebase/notebook, Synderilla
    reaction_id BIGINT NOT NULL REFERENCES designdb.reactions (id) ON DELETE CASCADE, -- Insert by codebase/notebook, Synderilla
    compound_id BIGINT NOT NULL REFERENCES designdb.compounds (id) ON DELETE RESTRICT, -- Insert by codebase/notebook, Synderilla
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_reactant UNIQUE (reaction_id, compound_id)
);

-- Registration identity: one row per distinct catalogue_smiles (unique). catalogue_inchikey is NOT unique —
-- different SMILES strings can map to the same Standard InChIKey after cartridge normalization (trigger from SMILES).
-- catalogue_hash: application / loader pipeline (same algorithm as compounds.compound_hash); NOT unique — multiple
-- rows may share a hash when SuperParent/registration hash collapses stereoisomers differently than stored SMILES and that's the correct behaviour
CREATE TABLE IF NOT EXISTS designdb.catalogue_compounds (
    id BIGSERIAL PRIMARY KEY,
    catalogue_smiles TEXT NOT NULL,
    catalogue_inchikey TEXT NOT NULL, -- Populated by RDKit cartridge trigger from catalogue_smiles
    catalogue_hash TEXT NOT NULL, -- Set by Enamine parsing script on insert/update; links via designdb.compound_catalogue_map
    rdkit_version TEXT NOT NULL,
    inchi_version TEXT NOT NULL,
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uq_catalogue_compounds_smiles UNIQUE (catalogue_smiles),
    CONSTRAINT ck_catalogue_compounds_hash_nonempty CHECK (length(trim(catalogue_hash)) > 0)
);

-- Former quotes rows split out. supplier = old quote_catalogue; supplier_id = old quote_entry; vendor = old quote_supplier.
CREATE TABLE IF NOT EXISTS designdb.catalogue_prices (
    id BIGSERIAL PRIMARY KEY,
    catalogue_id BIGINT NOT NULL REFERENCES designdb.catalogue_compounds (id) ON DELETE CASCADE,
    vendor TEXT NOT NULL,
    supplier TEXT,
    supplier_id TEXT NOT NULL,
    amount REAL NOT NULL,
    price REAL,
    currency TEXT,
    purity REAL,
    lead_time INTEGER,
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_catalogue_price UNIQUE (catalogue_id, vendor, supplier, supplier_id, amount)
);

-- Many-to-many: compounds - catalogue price lines matched on shared identity hash (compound_hash = catalogue_hash).
CREATE TABLE IF NOT EXISTS designdb.compound_catalogue_map (
    compound_id BIGINT NOT NULL REFERENCES designdb.compounds (id) ON DELETE CASCADE,
    catalogue_price_id BIGINT NOT NULL REFERENCES designdb.catalogue_prices (id) ON DELETE CASCADE,
    match_hash TEXT NOT NULL,
    -- Will be deleted automatically from db from here:
    -- HIPPO temporary columns, remove when HIPPO can handle catalogue_prices and catalogue_compounds:
    -- catalogue_inchikey TEXT, --Needs to be removed when HIPPO codebase ready to handle prices properly.
    -- supplier TEXT, --Needs to be removed when HIPPO codebase ready to handle prices properly.
    -- amount REAL, --Needs to be removed when HIPPO codebase ready to handle prices properly.
    -- price REAL, --Needs to be removed when HIPPO codebase ready to handle prices properly.
    -- lead_time INTEGER, --Needs to be removed when HIPPO codebase ready to handle prices properly.
    -- Will be deleted automatically until here
    created_on TIMESTAMPTZ NOT NULL DEFAULT now(),
    updated_on TIMESTAMPTZ NOT NULL DEFAULT now(),
    PRIMARY KEY (compound_id, catalogue_price_id),
    CONSTRAINT ck_compound_catalogue_map_match_hash_nonempty CHECK (length(trim(match_hash)) > 0)
);

CREATE TABLE IF NOT EXISTS designdb.scaffolds (
    id BIGSERIAL PRIMARY KEY,
    base_compound_id BIGINT REFERENCES designdb.compounds (id) ON DELETE SET NULL, -- Insert by codebase
    superstructure_compound_id BIGINT REFERENCES designdb.compounds (id) ON DELETE SET NULL, -- Insert by codebase
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_scaffold UNIQUE (base_compound_id, superstructure_compound_id)
);

CREATE TABLE IF NOT EXISTS designdb.routes (
    id BIGSERIAL PRIMARY KEY,
    product_compound_id BIGINT NOT NULL REFERENCES designdb.compounds (id) ON DELETE RESTRICT, -- Insert by codebase/notebook, Synderilla
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now()
);

CREATE TABLE IF NOT EXISTS designdb.components (
    id BIGSERIAL PRIMARY KEY,
    route_id BIGINT NOT NULL REFERENCES designdb.routes (id) ON DELETE RESTRICT, -- Insert by codebase/notebook, Synderilla
    component_type INTEGER, -- Insert by codebase/notebook, Synderilla--
    component_ref INTEGER, -- Insert by codebase/notebook, Synderilla
    component_amount REAL, -- Insert by codebase/notebook, Synderilla
    created_on TIMESTAMPTZ DEFAULT now(),
    updated_on TIMESTAMPTZ DEFAULT now(),
    CONSTRAINT uc_component UNIQUE (route_id, component_ref, component_type)
);

-- Replaced table with individual tag tables
-- CREATE TABLE IF NOT EXISTS designdb.tags (
--     id BIGSERIAL PRIMARY KEY,
--     tag_name TEXT, -- Insert by codebase
--     tag_description TEXT, -- New column
--     note TEXT, -- New column
--     -- compound_id BIGINT REFERENCES designdb.compounds (id) ON DELETE SET NULL, -- Insert by codebase, need to be removed, change in code needed
--     -- pose_id BIGINT REFERENCES designdb.poses (id) ON DELETE SET NULL, -- Insert by codebase, need to be removed, change in code needed
--     created_on TIMESTAMPTZ DEFAULT now(),
--     updated_on TIMESTAMPTZ DEFAULT now()
--     -- CONSTRAINT uc_tag_compound UNIQUE (tag_name, compound_id),
--     -- CONSTRAINT uc_tag_pose UNIQUE (tag_name, pose_id)
-- );


-- =========================================================
-- AUDIT TABLES
-- =========================================================

-- Event audit for catalogue_compounds (chemistry / identity rows)
CREATE TABLE IF NOT EXISTS designdb.catalogue_compounds_event_audit (
    audit_pk BIGSERIAL PRIMARY KEY,
    id BIGINT NOT NULL,
    operation CHAR(1) NOT NULL CHECK (operation IN ('I','U','D')),
    old_values JSONB,
    new_values JSONB,
    changed_by TEXT,
    changed_at TIMESTAMPTZ DEFAULT now()
);

-- Event audit for catalogue_prices (vendor / pricing lines)
CREATE TABLE IF NOT EXISTS designdb.catalogue_prices_event_audit (
    audit_pk BIGSERIAL PRIMARY KEY,
    id BIGINT NOT NULL,
    operation CHAR(1) NOT NULL CHECK (operation IN ('I','U','D')),
    old_values JSONB,
    new_values JSONB,
    changed_by TEXT,
    changed_at TIMESTAMPTZ DEFAULT now()
);

-- Event audit for pose_tags
CREATE TABLE IF NOT EXISTS designdb.pose_tags_event_audit (
    audit_pk BIGSERIAL PRIMARY KEY,
    id BIGINT NOT NULL,
    operation CHAR(1) NOT NULL CHECK (operation IN ('I','U','D')),
    old_values JSONB,
    new_values JSONB,
    changed_by TEXT,
    changed_at TIMESTAMPTZ DEFAULT now()
);

-- Event audit for compound_tags
CREATE TABLE IF NOT EXISTS designdb.compound_tags_event_audit (
    audit_pk BIGSERIAL PRIMARY KEY,
    id BIGINT NOT NULL,
    operation CHAR(1) NOT NULL CHECK (operation IN ('I','U','D')),
    old_values JSONB,
    new_values JSONB,
    changed_by TEXT,
    changed_at TIMESTAMPTZ DEFAULT now()
);

-- Event audit for pose_methods
CREATE TABLE IF NOT EXISTS designdb.pose_methods_event_audit (
    audit_pk BIGSERIAL PRIMARY KEY,
    id BIGINT NOT NULL,
    operation CHAR(1) NOT NULL CHECK (operation IN ('I','U','D')),
    old_values JSONB,
    new_values JSONB,
    changed_by TEXT,
    changed_at TIMESTAMPTZ DEFAULT now()
);

-- Event audit for enumeration_methods
CREATE TABLE IF NOT EXISTS designdb.enumeration_methods_event_audit (
    audit_pk BIGSERIAL PRIMARY KEY,
    id BIGINT NOT NULL,
    operation CHAR(1) NOT NULL CHECK (operation IN ('I','U','D')),
    old_values JSONB,
    new_values JSONB,
    changed_by TEXT,
    changed_at TIMESTAMPTZ DEFAULT now()
);

-- Event audit for scoring_methods
CREATE TABLE IF NOT EXISTS designdb.scoring_methods_event_audit (
    audit_pk BIGSERIAL PRIMARY KEY,
    id BIGINT NOT NULL,
    operation CHAR(1) NOT NULL CHECK (operation IN ('I','U','D')),
    old_values JSONB,
    new_values JSONB,
    changed_by TEXT,
    changed_at TIMESTAMPTZ DEFAULT now()
);

-- =========================================================
-- INDEXES
-- =========================================================

CREATE INDEX IF NOT EXISTS idx_target_name ON designdb.targets(target_name);
CREATE INDEX IF NOT EXISTS idx_target_created ON designdb.targets(created_on);

CREATE INDEX IF NOT EXISTS idx_scoring_method_name ON designdb.scoring_methods(method_name);
CREATE INDEX IF NOT EXISTS idx_scoring_method_created ON designdb.scoring_methods(created_on);

CREATE INDEX IF NOT EXISTS idx_enumeration_method_name ON designdb.enumeration_methods(enum_name);
CREATE INDEX IF NOT EXISTS idx_enumeration_method_created ON designdb.enumeration_methods(created_on);

CREATE INDEX IF NOT EXISTS idx_pose_method_name ON designdb.pose_methods(pose_method_name);
CREATE INDEX IF NOT EXISTS idx_pose_method_created ON designdb.pose_methods(created_on);

CREATE INDEX IF NOT EXISTS idx_compound_base_compound_id ON designdb.compounds(base_compound_id);
CREATE INDEX IF NOT EXISTS idx_compound_inchikey ON designdb.compounds(compound_inchikey);
CREATE INDEX IF NOT EXISTS idx_compound_smiles ON designdb.compounds(compound_smiles);
CREATE INDEX IF NOT EXISTS idx_compound_compound_hash ON designdb.compounds(compound_hash);
CREATE INDEX IF NOT EXISTS idx_compound_created ON designdb.compounds(created_on);

CREATE INDEX IF NOT EXISTS idx_feature_target_id ON designdb.features(target_id);
CREATE INDEX IF NOT EXISTS idx_feature_created ON designdb.features(created_on);

CREATE INDEX IF NOT EXISTS idx_route_product_compound_id ON designdb.routes(product_compound_id);
CREATE INDEX IF NOT EXISTS idx_route_created ON designdb.routes(created_on);

CREATE INDEX IF NOT EXISTS idx_reaction_product_compound_id ON designdb.reactions(product_compound_id);
CREATE INDEX IF NOT EXISTS idx_reaction_created ON designdb.reactions(created_on);

CREATE INDEX IF NOT EXISTS idx_pose_compound_id ON designdb.poses(compound_id);
CREATE INDEX IF NOT EXISTS idx_pose_target_id ON designdb.poses(target_id);
CREATE INDEX IF NOT EXISTS idx_pose_path ON designdb.poses(pose_path);
CREATE INDEX IF NOT EXISTS idx_pose_created ON designdb.poses(created_on);

CREATE INDEX IF NOT EXISTS idx_score_values_pose_id ON designdb.score_values(pose_id);
CREATE INDEX IF NOT EXISTS idx_score_values_compound_id ON designdb.score_values(compound_id);
CREATE INDEX IF NOT EXISTS idx_score_values_scoring_method_id ON designdb.score_values(scoring_method_id);
CREATE INDEX IF NOT EXISTS idx_score_values_created ON designdb.score_values(created_on);

CREATE INDEX IF NOT EXISTS idx_subsite_target_id ON designdb.subsites(target_id);
CREATE INDEX IF NOT EXISTS idx_subsite_created ON designdb.subsites(created_on);

CREATE INDEX IF NOT EXISTS idx_component_route_id ON designdb.components(route_id);
CREATE INDEX IF NOT EXISTS idx_component_created ON designdb.components(created_on);

CREATE INDEX IF NOT EXISTS idx_inspiration_original_pose_id ON designdb.inspirations(original_pose_id);
CREATE INDEX IF NOT EXISTS idx_inspiration_derivative_pose_id ON designdb.inspirations(derivative_pose_id);
CREATE INDEX IF NOT EXISTS idx_inspiration_created ON designdb.inspirations(created_on);

CREATE INDEX IF NOT EXISTS idx_interaction_feature_id ON designdb.interactions(feature_id);
CREATE INDEX IF NOT EXISTS idx_interaction_pose_id ON designdb.interactions(pose_id);
CREATE INDEX IF NOT EXISTS idx_interaction_created ON designdb.interactions(created_on);

-- UNIQUE(catalogue_smiles) supplies btree on catalogue_smiles; btree on catalogue_inchikey (non-unique) and catalogue_hash for lookups
CREATE INDEX IF NOT EXISTS idx_catalogue_compounds_created ON designdb.catalogue_compounds(created_on);
CREATE INDEX IF NOT EXISTS idx_catalogue_compounds_inchikey ON designdb.catalogue_compounds(catalogue_inchikey);
CREATE INDEX IF NOT EXISTS idx_catalogue_compounds_hash ON designdb.catalogue_compounds(catalogue_hash);

CREATE INDEX IF NOT EXISTS idx_catalogue_prices_catalogue_id ON designdb.catalogue_prices(catalogue_id);
CREATE INDEX IF NOT EXISTS idx_catalogue_prices_created ON designdb.catalogue_prices(created_on);

CREATE INDEX IF NOT EXISTS idx_compound_catalogue_map_match_hash ON designdb.compound_catalogue_map(match_hash);
CREATE INDEX IF NOT EXISTS idx_compound_catalogue_map_catalogue_price_id ON designdb.compound_catalogue_map(catalogue_price_id);
CREATE INDEX IF NOT EXISTS idx_compound_catalogue_map_created ON designdb.compound_catalogue_map(created_on);

-- =========================================================
-- AUDIT INDEXES
-- =========================================================

CREATE INDEX IF NOT EXISTS idx_catalogue_compounds_event_audit_id ON designdb.catalogue_compounds_event_audit(id);
CREATE INDEX IF NOT EXISTS idx_catalogue_compounds_event_audit_operation ON designdb.catalogue_compounds_event_audit(operation);
CREATE INDEX IF NOT EXISTS idx_catalogue_compounds_event_audit_changed_at ON designdb.catalogue_compounds_event_audit(changed_at);
CREATE INDEX IF NOT EXISTS idx_catalogue_compounds_event_audit_changed_by ON designdb.catalogue_compounds_event_audit(changed_by);
CREATE INDEX IF NOT EXISTS idx_catalogue_compounds_event_audit_old_gin ON designdb.catalogue_compounds_event_audit USING GIN (old_values);
CREATE INDEX IF NOT EXISTS idx_catalogue_compounds_event_audit_new_gin ON designdb.catalogue_compounds_event_audit USING GIN (new_values);
CREATE INDEX IF NOT EXISTS idx_catalogue_compounds_event_audit_id_changed ON designdb.catalogue_compounds_event_audit(id, changed_at DESC);

CREATE INDEX IF NOT EXISTS idx_catalogue_prices_event_audit_id ON designdb.catalogue_prices_event_audit(id);
CREATE INDEX IF NOT EXISTS idx_catalogue_prices_event_audit_operation ON designdb.catalogue_prices_event_audit(operation);
CREATE INDEX IF NOT EXISTS idx_catalogue_prices_event_audit_changed_at ON designdb.catalogue_prices_event_audit(changed_at);
CREATE INDEX IF NOT EXISTS idx_catalogue_prices_event_audit_changed_by ON designdb.catalogue_prices_event_audit(changed_by);
CREATE INDEX IF NOT EXISTS idx_catalogue_prices_event_audit_old_gin ON designdb.catalogue_prices_event_audit USING GIN (old_values);
CREATE INDEX IF NOT EXISTS idx_catalogue_prices_event_audit_new_gin ON designdb.catalogue_prices_event_audit USING GIN (new_values);
CREATE INDEX IF NOT EXISTS idx_catalogue_prices_event_audit_id_changed ON designdb.catalogue_prices_event_audit(id, changed_at DESC);

CREATE INDEX IF NOT EXISTS idx_pose_tags_event_audit_id ON designdb.pose_tags_event_audit(id);
CREATE INDEX IF NOT EXISTS idx_pose_tags_event_audit_operation ON designdb.pose_tags_event_audit(operation);
CREATE INDEX IF NOT EXISTS idx_pose_tags_event_audit_changed_at ON designdb.pose_tags_event_audit(changed_at);
CREATE INDEX IF NOT EXISTS idx_pose_tags_event_audit_old_gin ON designdb.pose_tags_event_audit USING GIN (old_values);
CREATE INDEX IF NOT EXISTS idx_pose_tags_event_audit_new_gin ON designdb.pose_tags_event_audit USING GIN (new_values);
CREATE INDEX IF NOT EXISTS idx_pose_tags_event_audit_id_changed ON designdb.pose_tags_event_audit(id, changed_at DESC);

CREATE INDEX IF NOT EXISTS idx_compound_tags_event_audit_id ON designdb.compound_tags_event_audit(id);
CREATE INDEX IF NOT EXISTS idx_compound_tags_event_audit_operation ON designdb.compound_tags_event_audit(operation);
CREATE INDEX IF NOT EXISTS idx_compound_tags_event_audit_changed_at ON designdb.compound_tags_event_audit(changed_at);
CREATE INDEX IF NOT EXISTS idx_compound_tags_event_audit_old_gin ON designdb.compound_tags_event_audit USING GIN (old_values);
CREATE INDEX IF NOT EXISTS idx_compound_tags_event_audit_new_gin ON designdb.compound_tags_event_audit USING GIN (new_values);
CREATE INDEX IF NOT EXISTS idx_compound_tags_event_audit_id_changed ON designdb.compound_tags_event_audit(id, changed_at DESC);

CREATE INDEX IF NOT EXISTS idx_pose_methods_event_audit_id ON designdb.pose_methods_event_audit(id);
CREATE INDEX IF NOT EXISTS idx_pose_methods_event_audit_operation ON designdb.pose_methods_event_audit(operation);
CREATE INDEX IF NOT EXISTS idx_pose_methods_event_audit_changed_at ON designdb.pose_methods_event_audit(changed_at);
CREATE INDEX IF NOT EXISTS idx_pose_methods_event_audit_old_gin ON designdb.pose_methods_event_audit USING GIN (old_values);
CREATE INDEX IF NOT EXISTS idx_pose_methods_event_audit_new_gin ON designdb.pose_methods_event_audit USING GIN (new_values);
CREATE INDEX IF NOT EXISTS idx_pose_methods_event_audit_id_changed ON designdb.pose_methods_event_audit(id, changed_at DESC);

CREATE INDEX IF NOT EXISTS idx_enumeration_methods_event_audit_id ON designdb.enumeration_methods_event_audit(id);
CREATE INDEX IF NOT EXISTS idx_enumeration_methods_event_audit_operation ON designdb.enumeration_methods_event_audit(operation);
CREATE INDEX IF NOT EXISTS idx_enumeration_methods_event_audit_changed_at ON designdb.enumeration_methods_event_audit(changed_at);
CREATE INDEX IF NOT EXISTS idx_enumeration_methods_event_audit_old_gin ON designdb.enumeration_methods_event_audit USING GIN (old_values);
CREATE INDEX IF NOT EXISTS idx_enumeration_methods_event_audit_new_gin ON designdb.enumeration_methods_event_audit USING GIN (new_values);
CREATE INDEX IF NOT EXISTS idx_enumeration_methods_event_audit_id_changed ON designdb.enumeration_methods_event_audit(id, changed_at DESC);

CREATE INDEX IF NOT EXISTS idx_scoring_methods_event_audit_id ON designdb.scoring_methods_event_audit(id);
CREATE INDEX IF NOT EXISTS idx_scoring_methods_event_audit_operation ON designdb.scoring_methods_event_audit(operation);
CREATE INDEX IF NOT EXISTS idx_scoring_methods_event_audit_changed_at ON designdb.scoring_methods_event_audit(changed_at);
CREATE INDEX IF NOT EXISTS idx_scoring_methods_event_audit_old_gin ON designdb.scoring_methods_event_audit USING GIN (old_values);
CREATE INDEX IF NOT EXISTS idx_scoring_methods_event_audit_new_gin ON designdb.scoring_methods_event_audit USING GIN (new_values);
CREATE INDEX IF NOT EXISTS idx_scoring_methods_event_audit_id_changed ON designdb.scoring_methods_event_audit(id, changed_at DESC);

CREATE INDEX IF NOT EXISTS idx_reactant_reaction_id ON designdb.reactants(reaction_id);
CREATE INDEX IF NOT EXISTS idx_reactant_compound_id ON designdb.reactants(compound_id);
CREATE INDEX IF NOT EXISTS idx_reactant_created ON designdb.reactants(created_on);

CREATE INDEX IF NOT EXISTS idx_scaffold_base_compound_id ON designdb.scaffolds(base_compound_id);
CREATE INDEX IF NOT EXISTS idx_scaffold_superstructure_compound_id ON designdb.scaffolds(superstructure_compound_id);
CREATE INDEX IF NOT EXISTS idx_scaffold_created ON designdb.scaffolds(created_on);

CREATE INDEX IF NOT EXISTS idx_subsite_tag_subsite_id ON designdb.subsite_tags(subsite_id);
CREATE INDEX IF NOT EXISTS idx_subsite_tag_pose_id ON designdb.subsite_tags(pose_id);
CREATE INDEX IF NOT EXISTS idx_subsite_tag_created ON designdb.subsite_tags(created_on);

-- Removed due to replaced tables
-- CREATE INDEX IF NOT EXISTS idx_tag_compound ON designdb.tags(tag_compound);
-- CREATE INDEX IF NOT EXISTS idx_tag_pose ON designdb.tags(tag_pose);
-- CREATE INDEX IF NOT EXISTS idx_tag_created ON designdb.tags(created_on);

CREATE INDEX IF NOT EXISTS idx_pose_tag_created ON designdb.pose_tags(created_on);
CREATE INDEX IF NOT EXISTS idx_compound_tag_created ON designdb.compound_tags(created_on);

CREATE INDEX IF NOT EXISTS idx_has_pose_methods_pose_method_id ON designdb.has_pose_methods(pose_method_id);
CREATE INDEX IF NOT EXISTS idx_has_pose_methods_created ON designdb.has_pose_methods(created_on);

CREATE INDEX IF NOT EXISTS idx_has_pose_tag_pose_tag_id ON designdb.has_pose_tags(pose_tag_id);
CREATE INDEX IF NOT EXISTS idx_has_pose_tag_created ON designdb.has_pose_tags(created_on);

CREATE INDEX IF NOT EXISTS idx_has_enumeration_methods_enumeration_method_id ON designdb.has_enumeration_methods(enumeration_method_id);
CREATE INDEX IF NOT EXISTS idx_has_enumeration_methods_created ON designdb.has_enumeration_methods(created_on);

CREATE INDEX IF NOT EXISTS idx_has_compound_tag_compound_tag_id ON designdb.has_compound_tags(compound_tag_id);
CREATE INDEX IF NOT EXISTS idx_has_compound_tag_created ON designdb.has_compound_tags(created_on);

-- =========================================================
-- MATERIALIZED VIEWS
-- =========================================================
-- designdb.scores_per_pose_pivoted_mv: pose_id, compound_id + one column per (method_name, method_version).
-- Pivoted from score_values joined with scoring_methods. Dynamically re-generated when new method added.

-- =========================================================
-- VIEWS
-- =========================================================

-- Price-line UPDATEs from catalogue_prices_event_audit (JSON keys match catalogue_prices column names)
CREATE OR REPLACE VIEW designdb.catalogue_prices_price_changes_v AS
SELECT
  a.id AS catalogue_price_id,
  (NULLIF(COALESCE(a.new_values->>'catalogue_id', a.old_values->>'catalogue_id'), ''))::BIGINT AS catalogue_id,
  COALESCE(a.new_values->>'vendor', a.old_values->>'vendor') AS vendor,
  COALESCE(a.new_values->>'supplier', a.old_values->>'supplier') AS supplier,
  COALESCE(a.new_values->>'supplier_id', a.old_values->>'supplier_id') AS supplier_id,
  (NULLIF(COALESCE(a.new_values->>'amount', a.old_values->>'amount'), ''))::DOUBLE PRECISION AS amount,
  (NULLIF(a.old_values->>'price', ''))::DOUBLE PRECISION AS price_old,
  (NULLIF(a.new_values->>'price', ''))::DOUBLE PRECISION AS price_new,
  COALESCE(a.new_values->>'currency', a.old_values->>'currency') AS currency,
  (NULLIF(COALESCE(a.new_values->>'purity', a.old_values->>'purity'), ''))::DOUBLE PRECISION AS purity,
  (NULLIF(COALESCE(a.new_values->>'lead_time', a.old_values->>'lead_time'), ''))::INTEGER AS lead_time,
  a.changed_at
FROM designdb.catalogue_prices_event_audit a
WHERE a.operation = 'U';

-- Pose tags: UPDATE events with old/new name, description, note.
CREATE OR REPLACE VIEW designdb.pose_tags_changes_v AS
SELECT
  a.id AS pose_tag_id,
  a.old_values->>'pose_tag_name' AS pose_tag_name_old,
  a.new_values->>'pose_tag_name' AS pose_tag_name_new,
  a.old_values->>'pose_tag_description' AS pose_tag_description_old,
  a.new_values->>'pose_tag_description' AS pose_tag_description_new,
  a.old_values->>'pose_tag_note' AS pose_tag_note_old,
  a.new_values->>'pose_tag_note' AS pose_tag_note_new,
  a.changed_by,
  a.changed_at
FROM designdb.pose_tags_event_audit a
WHERE a.operation = 'U';

-- Compound tags: UPDATE events with old/new name, description, note.
CREATE OR REPLACE VIEW designdb.compound_tags_changes_v AS
SELECT
  a.id AS compound_tag_id,
  a.old_values->>'compound_tag_name' AS compound_tag_name_old,
  a.new_values->>'compound_tag_name' AS compound_tag_name_new,
  a.old_values->>'compound_tag_description' AS compound_tag_description_old,
  a.new_values->>'compound_tag_description' AS compound_tag_description_new,
  a.old_values->>'compound_tag_note' AS compound_tag_note_old,
  a.new_values->>'compound_tag_note' AS compound_tag_note_new,
  a.changed_by,
  a.changed_at
FROM designdb.compound_tags_event_audit a
WHERE a.operation = 'U';

-- Pose methods: UPDATE events with old/new name, description, version, etc.
CREATE OR REPLACE VIEW designdb.pose_methods_changes_v AS
SELECT
  a.id AS pose_method_id,
  a.old_values->>'pose_method_name' AS pose_method_name_old,
  a.new_values->>'pose_method_name' AS pose_method_name_new,
  a.old_values->>'pose_method_description' AS pose_method_description_old,
  a.new_values->>'pose_method_description' AS pose_method_description_new,
  a.old_values->>'pose_method_version' AS pose_method_version_old,
  a.new_values->>'pose_method_version' AS pose_method_version_new,
  a.changed_by,
  a.changed_at
FROM designdb.pose_methods_event_audit a
WHERE a.operation = 'U';

-- Enumeration methods: UPDATE events with old/new name, description, version, etc.
CREATE OR REPLACE VIEW designdb.enumeration_methods_changes_v AS
SELECT
  a.id AS enumeration_method_id,
  a.old_values->>'enum_name' AS enum_name_old,
  a.new_values->>'enum_name' AS enum_name_new,
  a.old_values->>'enum_description' AS enum_description_old,
  a.new_values->>'enum_description' AS enum_description_new,
  a.old_values->>'enum_version' AS enum_version_old,
  a.new_values->>'enum_version' AS enum_version_new,
  a.changed_by,
  a.changed_at
FROM designdb.enumeration_methods_event_audit a
WHERE a.operation = 'U';

-- Scoring methods: UPDATE events with old/new name, description, version, etc.
CREATE OR REPLACE VIEW designdb.scoring_methods_changes_v AS
SELECT
  a.id AS scoring_method_id,
  a.old_values->>'method_name' AS method_name_old,
  a.new_values->>'method_name' AS method_name_new,
  a.old_values->>'method_description' AS method_description_old,
  a.new_values->>'method_description' AS method_description_new,
  a.old_values->>'method_version' AS method_version_old,
  a.new_values->>'method_version' AS method_version_new,
  a.changed_by,
  a.changed_at
FROM designdb.scoring_methods_event_audit a
WHERE a.operation = 'U';

-- =========================================================
-- FUNCTIONS
-- =========================================================

-- =========================================================
-- RDKIT CARTRIDGE – COMPOUND WRAPPERS
-- =========================================================
-- Schema-qualified wrappers for RDKit mol_from_smiles, mol_to_smiles, mol_inchikey, mol_to_ctab (compound, pose, catalogue_compounds triggers).

CREATE OR REPLACE FUNCTION designdb.mol_from_smiles(smiles TEXT) RETURNS rdkit.mol
  LANGUAGE SQL AS $$ SELECT rdkit.mol_from_smiles(smiles::cstring); $$;

CREATE OR REPLACE FUNCTION designdb.mol_to_smiles(m rdkit.mol) RETURNS text
  LANGUAGE SQL AS $$ SELECT rdkit.mol_to_smiles(m); $$;

CREATE OR REPLACE FUNCTION designdb.mol_to_inchikey(m rdkit.mol) RETURNS text
  LANGUAGE SQL AS $$ SELECT rdkit.mol_inchikey(m); $$;

CREATE OR REPLACE FUNCTION designdb.mol_to_ctab(m rdkit.mol) RETURNS text
  LANGUAGE SQL AS $$ SELECT rdkit.mol_to_ctab(m); $$;

-- =========================================================
-- RDKIT CARTRIDGE – COMPOUND TRIGGER
-- =========================================================
-- Input: compound_smiles (inserted by application). Populates compound_mol (CTAB) and compound_inchikey.

CREATE OR REPLACE FUNCTION designdb.populate_compound_cartridge_from_smiles()
RETURNS trigger
LANGUAGE plpgsql
AS $$
DECLARE
  v_mol rdkit.mol;
BEGIN
  IF NEW.compound_smiles IS NOT NULL THEN
    BEGIN
      v_mol := designdb.mol_from_smiles(NEW.compound_smiles);
      IF v_mol IS NOT NULL THEN
        -- NEW.compound_mol := v_mol; -- store rdkit.mol (was default text form ~ SMILES over JDBC)
        NEW.compound_mol := designdb.mol_to_ctab(v_mol);
        NEW.compound_inchikey := designdb.mol_to_inchikey(v_mol);
      END IF;
    EXCEPTION WHEN OTHERS THEN
      NULL;
    END;
  END IF;
  RETURN NEW;
END;
$$;

DROP TRIGGER IF EXISTS trg_populate_compound_cartridge_from_smiles ON designdb.compounds;
CREATE TRIGGER trg_populate_compound_cartridge_from_smiles
  BEFORE INSERT OR UPDATE OF compound_smiles ON designdb.compounds
  FOR EACH ROW
  EXECUTE FUNCTION designdb.populate_compound_cartridge_from_smiles();

-- =========================================================
-- RDKIT CARTRIDGE – CATALOGUE_COMPOUNDS TRIGGER
-- =========================================================
-- Input: catalogue_smiles (inserted by application). Populates catalogue_inchikey (NOT NULL column).

CREATE OR REPLACE FUNCTION designdb.populate_catalogue_cartridge_from_smiles()
RETURNS trigger
LANGUAGE plpgsql
AS $$
DECLARE
  v_mol rdkit.mol;
  v_ik TEXT;
BEGIN
  IF NEW.catalogue_smiles IS NULL OR btrim(NEW.catalogue_smiles) = '' THEN
    RAISE EXCEPTION 'designdb.catalogue_compounds: catalogue_smiles is required';
  END IF;
  v_mol := designdb.mol_from_smiles(NEW.catalogue_smiles);
  IF v_mol IS NULL THEN
    RAISE EXCEPTION 'designdb.catalogue_compounds: mol_from_smiles returned NULL for catalogue_smiles';
  END IF;
  v_ik := designdb.mol_to_inchikey(v_mol);
  IF v_ik IS NULL OR btrim(v_ik) = '' THEN
    RAISE EXCEPTION 'designdb.catalogue_compounds: mol_to_inchikey returned empty for catalogue_smiles';
  END IF;
  NEW.catalogue_inchikey := v_ik;
  RETURN NEW;
END;
$$;

DROP TRIGGER IF EXISTS trg_populate_catalogue_cartridge_from_smiles ON designdb.catalogue_compounds;
CREATE TRIGGER trg_populate_catalogue_cartridge_from_smiles
  BEFORE INSERT OR UPDATE OF catalogue_smiles ON designdb.catalogue_compounds
  FOR EACH ROW
  EXECUTE FUNCTION designdb.populate_catalogue_cartridge_from_smiles();

-- =========================================================
-- RDKIT CARTRIDGE – POSE TRIGGER
-- =========================================================
-- Input: pose_mol (inserted by application). Populates pose_inchikey and pose_smiles.

CREATE OR REPLACE FUNCTION designdb.populate_pose_cartridge_from_mol()
RETURNS trigger
LANGUAGE plpgsql
AS $$
BEGIN
  IF NEW.pose_mol IS NOT NULL THEN
    BEGIN
      NEW.pose_smiles := designdb.mol_to_smiles(NEW.pose_mol);
      NEW.pose_inchikey := designdb.mol_to_inchikey(NEW.pose_mol);
    EXCEPTION WHEN OTHERS THEN
      NULL;
    END;
  END IF;
  RETURN NEW;
END;
$$;

DROP TRIGGER IF EXISTS trg_populate_pose_cartridge_from_mol ON designdb.poses;
CREATE TRIGGER trg_populate_pose_cartridge_from_mol
  BEFORE INSERT OR UPDATE OF pose_mol ON designdb.poses
  FOR EACH ROW
  EXECUTE FUNCTION designdb.populate_pose_cartridge_from_mol();

-- =========================================================
-- SCORE_VALUES – ENFORCE compound_id MATCHES pose's compound_id
-- =========================================================
CREATE OR REPLACE FUNCTION designdb.check_score_values_compound_matches_pose()
RETURNS trigger
LANGUAGE plpgsql
AS $$
DECLARE
  v_pose_compound_id BIGINT;
BEGIN
  SELECT compound_id INTO v_pose_compound_id
  FROM designdb.poses
  WHERE id = NEW.pose_id;
  IF v_pose_compound_id IS NULL THEN
    RAISE EXCEPTION 'pose_id % does not exist', NEW.pose_id;
  END IF;
  IF NEW.compound_id IS DISTINCT FROM v_pose_compound_id THEN
    RAISE EXCEPTION 'score_values.compound_id (%) must match poses.compound_id (%) for pose_id %',
      NEW.compound_id, v_pose_compound_id, NEW.pose_id;
  END IF;
  RETURN NEW;
END;
$$;

DROP TRIGGER IF EXISTS trg_check_score_values_compound_matches_pose ON designdb.score_values;
CREATE TRIGGER trg_check_score_values_compound_matches_pose
  BEFORE INSERT OR UPDATE OF pose_id, compound_id ON designdb.score_values
  FOR EACH ROW
  EXECUTE FUNCTION designdb.check_score_values_compound_matches_pose();

-- =========================================================
-- (Re)creates materialized view designdb.scores_per_pose_pivoted_mv with columns from scoring_method.
-- Columns: pose_id, compound_id, then one JSONB column per (method_name, method_version).
-- Column names use suffix _m{scoring_method_id} to avoid collisions (e.g. vina_1_0_m1).
-- Value in column: if score is numeric then JSONB number, if text then JSONB string.
CREATE OR REPLACE FUNCTION designdb.create_scores_per_pose_pivoted_mv()
RETURNS void
LANGUAGE plpgsql
AS $$
DECLARE
  select_qry text;
  col text;
  method_rec record;
  score_txt text;
  value_expr text;
  numeric_pat text := '^\-?[0-9]*\.?[0-9]+([eE][\-+]?[0-9]+)?$';
BEGIN
  select_qry := 'SELECT sv.pose_id, sv.compound_id';
  FOR method_rec IN
    SELECT m.id, m.method_name, m.method_version
    FROM designdb.scoring_methods m
    ORDER BY m.id
  LOOP
    col := regexp_replace(
      trim(method_rec.method_name) || '_' || coalesce(
        replace(replace(trim(coalesce(method_rec.method_version, '')), ' ', '_'), '.', '_'),
        ''
      ),
      '[^a-zA-Z0-9_]', '_', 'g'
    ) || '_m' || method_rec.id;
    IF col <> '' AND col <> '_' THEN
      col := quote_ident(col);
      score_txt := '(sv.score->>' || quote_literal('score') || ')';
      value_expr := '(CASE WHEN ' || score_txt || ' IS NOT NULL AND ' || score_txt || ' ~ ' || quote_literal(numeric_pat)
        || ' THEN to_jsonb((' || score_txt || ')::numeric) ELSE to_jsonb(' || score_txt || ') END)';
      select_qry := select_qry || ', MAX(CASE WHEN sv.scoring_method_id = ' || method_rec.id
        || ' THEN ' || value_expr || ' END) AS ' || col;
    END IF;
  END LOOP;
  select_qry := select_qry || ' FROM designdb.score_values sv GROUP BY sv.pose_id, sv.compound_id';
  EXECUTE 'DROP MATERIALIZED VIEW IF EXISTS designdb.scores_per_pose_pivoted_mv CASCADE';
  EXECUTE 'CREATE MATERIALIZED VIEW designdb.scores_per_pose_pivoted_mv AS ' || select_qry;
  EXECUTE 'CREATE UNIQUE INDEX ON designdb.scores_per_pose_pivoted_mv (pose_id, compound_id)';
END;
$$;

CREATE OR REPLACE FUNCTION designdb.trg_recreate_scores_pivoted_mv()
RETURNS trigger
LANGUAGE plpgsql
AS $$
BEGIN
  PERFORM designdb.create_scores_per_pose_pivoted_mv();
  RETURN NULL;
END;
$$;

CREATE OR REPLACE FUNCTION designdb.refresh_scores_per_pose_pivoted_mv()
RETURNS trigger
LANGUAGE plpgsql
AS $$
BEGIN
  REFRESH MATERIALIZED VIEW CONCURRENTLY designdb.scores_per_pose_pivoted_mv;
  RETURN NULL;
END;
$$;

CREATE OR REPLACE FUNCTION designdb.update_updated_on()
RETURNS trigger AS $$
BEGIN
    NEW.updated_on = now();
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

/*
-- Will be deleted automatically from db from here:
-- === HIPPO TEMPORARY: denormalized map columns; remove when HIPPO ready ===

CREATE OR REPLACE FUNCTION designdb.trg_compound_catalogue_map_from_compound_temp()
RETURNS trigger
LANGUAGE plpgsql
AS $$
BEGIN
    IF TG_OP = 'UPDATE' THEN
        IF OLD.compound_hash IS NOT DISTINCT FROM NEW.compound_hash THEN
            RETURN NEW;
        END IF;
        DELETE FROM designdb.compound_catalogue_map WHERE compound_id = NEW.id;
    END IF;

    INSERT INTO designdb.compound_catalogue_map (
        compound_id, catalogue_price_id, match_hash,
        catalogue_inchikey, supplier, amount, price, lead_time
    )
    SELECT NEW.id, cp.id, NEW.compound_hash,
           cat.catalogue_inchikey, cp.supplier, cp.amount, cp.price, cp.lead_time
    FROM designdb.catalogue_prices cp
    JOIN designdb.catalogue_compounds cat ON cat.id = cp.catalogue_id
    WHERE cat.catalogue_hash = NEW.compound_hash
    ON CONFLICT (compound_id, catalogue_price_id) DO NOTHING;

    RETURN NEW;
END;
$$;

CREATE OR REPLACE FUNCTION designdb.trg_compound_catalogue_map_from_catalogue_temp()
RETURNS trigger
LANGUAGE plpgsql
AS $$
BEGIN
    IF TG_OP = 'UPDATE' THEN
        IF OLD.catalogue_hash IS NOT DISTINCT FROM NEW.catalogue_hash THEN
            RETURN NEW;
        END IF;
        DELETE FROM designdb.compound_catalogue_map
        WHERE catalogue_price_id IN (
            SELECT id FROM designdb.catalogue_prices WHERE catalogue_id = NEW.id
        );
    END IF;

    INSERT INTO designdb.compound_catalogue_map (
        compound_id, catalogue_price_id, match_hash,
        catalogue_inchikey, supplier, amount, price, lead_time
    )
    SELECT c.id, cp.id, c.compound_hash,
           NEW.catalogue_inchikey, cp.supplier, cp.amount, cp.price, cp.lead_time
    FROM designdb.compounds c
    CROSS JOIN designdb.catalogue_prices cp
    WHERE cp.catalogue_id = NEW.id AND c.compound_hash = NEW.catalogue_hash
    ON CONFLICT (compound_id, catalogue_price_id) DO NOTHING;

    RETURN NEW;
END;
$$;

CREATE OR REPLACE FUNCTION designdb.trg_compound_catalogue_map_from_catalogue_price_temp()
RETURNS trigger
LANGUAGE plpgsql
AS $$
BEGIN
    IF TG_OP = 'UPDATE' AND OLD.catalogue_id IS NOT DISTINCT FROM NEW.catalogue_id THEN
        RETURN NEW;
    END IF;
    IF TG_OP = 'UPDATE' THEN
        DELETE FROM designdb.compound_catalogue_map WHERE catalogue_price_id = NEW.id;
    END IF;

    INSERT INTO designdb.compound_catalogue_map (
        compound_id, catalogue_price_id, match_hash,
        catalogue_inchikey, supplier, amount, price, lead_time
    )
    SELECT c.id, NEW.id, c.compound_hash,
           cat.catalogue_inchikey, NEW.supplier, NEW.amount, NEW.price, NEW.lead_time
    FROM designdb.compounds c
    JOIN designdb.catalogue_compounds cat ON cat.id = NEW.catalogue_id
    WHERE c.compound_hash = cat.catalogue_hash
    ON CONFLICT (compound_id, catalogue_price_id) DO NOTHING;

    RETURN NEW;
END;
$$;

CREATE OR REPLACE FUNCTION designdb.trg_compound_catalogue_map_refresh_denorm_from_catalogue_price_temp()
RETURNS trigger
LANGUAGE plpgsql
AS $$
BEGIN
    UPDATE designdb.compound_catalogue_map m
    SET supplier = NEW.supplier,
        amount = NEW.amount,
        price = NEW.price,
        lead_time = NEW.lead_time,
        updated_on = now()
    WHERE m.catalogue_price_id = NEW.id;
    RETURN NEW;
END;
$$;

CREATE OR REPLACE FUNCTION designdb.trg_compound_catalogue_map_refresh_denorm_from_catalogue_compound_temp()
RETURNS trigger
LANGUAGE plpgsql
AS $$
BEGIN
    IF TG_OP = 'UPDATE' THEN
        IF OLD.catalogue_inchikey IS NOT DISTINCT FROM NEW.catalogue_inchikey THEN
            RETURN NEW;
        END IF;
    END IF;
    UPDATE designdb.compound_catalogue_map m
    SET catalogue_inchikey = NEW.catalogue_inchikey,
        updated_on = now()
    WHERE m.catalogue_price_id IN (
        SELECT id FROM designdb.catalogue_prices WHERE catalogue_id = NEW.id
    );
    RETURN NEW;
END;
$$;
-- Will be deleted automatically until here
*/

-- =========================================================
-- AUDIT FUNCTIONS
-- =========================================================

-- Event audit trigger: records INSERT/UPDATE/DELETE to an audit table with JSONB old/new values.
-- TG_ARGV[0]=audit_table, [1]=pk_column, [2]=excluded_columns (comma-sep), [3]=binary_columns (comma-sep, hashed as sha256:hex).
CREATE OR REPLACE FUNCTION designdb.event_audit_trigger()
RETURNS trigger AS
$$
DECLARE
    v_audit_table TEXT := TG_ARGV[0];
    v_pk_col TEXT := TG_ARGV[1];
    v_excluded TEXT := COALESCE(TG_ARGV[2], '');
    v_bincols TEXT := COALESCE(TG_ARGV[3], '');
    v_excluded_arr TEXT[];
    v_bincols_arr TEXT[];
    v_changed_by TEXT := COALESCE(current_setting('app.current_user', true), current_user);
    v_old_json JSONB;
    v_new_json JSONB;
    v_record_pk TEXT;
    v_bin_hash TEXT;
    v_col TEXT;
BEGIN
    IF v_excluded = '' THEN
        v_excluded_arr := ARRAY[]::text[];
    ELSE
        v_excluded_arr := ARRAY(SELECT trim(x) FROM regexp_split_to_table(v_excluded, ',') AS x);
    END IF;

    IF v_bincols = '' THEN
        v_bincols_arr := ARRAY[]::text[];
    ELSE
        v_bincols_arr := ARRAY(SELECT trim(x) FROM regexp_split_to_table(v_bincols, ',') AS x);
    END IF;

    IF TG_OP = 'INSERT' THEN
        v_old_json := NULL;
        v_new_json := to_jsonb(NEW);
    ELSIF TG_OP = 'DELETE' THEN
        v_old_json := to_jsonb(OLD);
        v_new_json := NULL;
    ELSE
        IF to_jsonb(OLD) IS NOT DISTINCT FROM to_jsonb(NEW) THEN
            RETURN NEW;
        END IF;
        v_old_json := to_jsonb(OLD);
        v_new_json := to_jsonb(NEW);
    END IF;

    FOREACH v_col IN ARRAY v_excluded_arr LOOP
        IF v_old_json IS NOT NULL THEN v_old_json := v_old_json - v_col; END IF;
        IF v_new_json IS NOT NULL THEN v_new_json := v_new_json - v_col; END IF;
    END LOOP;

    FOREACH v_col IN ARRAY v_bincols_arr LOOP
        IF v_old_json IS NOT NULL AND v_old_json ? v_col THEN
            BEGIN
                EXECUTE format('SELECT CASE WHEN ($1).%I IS NULL THEN NULL ELSE encode(digest(($1).%I::bytea, ''sha256''), ''hex'') END', v_col, v_col)
                USING OLD INTO v_bin_hash;
                IF v_bin_hash IS NOT NULL THEN
                    v_old_json := jsonb_set(v_old_json, ARRAY[v_col], to_jsonb('sha256:' || v_bin_hash));
                ELSE
                    v_old_json := v_old_json - v_col;
                END IF;
            EXCEPTION WHEN others THEN
                v_old_json := v_old_json - v_col;
            END;
        END IF;
        IF v_new_json IS NOT NULL AND v_new_json ? v_col THEN
            BEGIN
                EXECUTE format('SELECT CASE WHEN ($1).%I IS NULL THEN NULL ELSE encode(digest(($1).%I::bytea, ''sha256''), ''hex'') END', v_col, v_col)
                USING NEW INTO v_bin_hash;
                IF v_bin_hash IS NOT NULL THEN
                    v_new_json := jsonb_set(v_new_json, ARRAY[v_col], to_jsonb('sha256:' || v_bin_hash));
                ELSE
                    v_new_json := v_new_json - v_col;
                END IF;
            EXCEPTION WHEN others THEN
                v_new_json := v_new_json - v_col;
            END;
        END IF;
    END LOOP;

    IF TG_OP = 'DELETE' THEN
        EXECUTE format('SELECT ($1).%I::text', v_pk_col) USING OLD INTO v_record_pk;
    ELSE
        EXECUTE format('SELECT ($1).%I::text', v_pk_col) USING NEW INTO v_record_pk;
    END IF;

    EXECUTE format('INSERT INTO %s (%s, operation, old_values, new_values, changed_by, changed_at) VALUES ($1,$2,$3,$4,$5,NOW())', v_audit_table, v_pk_col)
    USING v_record_pk::BIGINT, TG_OP::CHAR, v_old_json, v_new_json, v_changed_by;

    RETURN COALESCE(NEW, OLD);
END;
$$ LANGUAGE plpgsql VOLATILE;

-- =========================================================
-- TRIGGERS (updated_on)
-- =========================================================

DROP TRIGGER IF EXISTS trg_target_updated_on ON designdb.targets;
CREATE TRIGGER trg_target_updated_on BEFORE UPDATE ON designdb.targets FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_scoring_method_updated_on ON designdb.scoring_methods;
CREATE TRIGGER trg_scoring_method_updated_on BEFORE UPDATE ON designdb.scoring_methods FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_scoring_method_recreate_pivoted_mv ON designdb.scoring_methods;
CREATE TRIGGER trg_scoring_method_recreate_pivoted_mv
  AFTER INSERT OR UPDATE OR DELETE ON designdb.scoring_methods
  FOR EACH STATEMENT EXECUTE FUNCTION designdb.trg_recreate_scores_pivoted_mv();

DROP TRIGGER IF EXISTS trg_enumeration_method_updated_on ON designdb.enumeration_methods;
CREATE TRIGGER trg_enumeration_method_updated_on BEFORE UPDATE ON designdb.enumeration_methods FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_pose_method_updated_on ON designdb.pose_methods;
CREATE TRIGGER trg_pose_method_updated_on BEFORE UPDATE ON designdb.pose_methods FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_compound_updated_on ON designdb.compounds;
CREATE TRIGGER trg_compound_updated_on BEFORE UPDATE ON designdb.compounds FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

/*
-- Will be deleted automatically from db from here:
-- Production trigger names reserved (see /* */ block at end of file); active temp triggers:
DROP TRIGGER IF EXISTS trg_compound_catalogue_map_sync_compound ON designdb.compounds;
DROP TRIGGER IF EXISTS trg_compound_catalogue_map_sync_compound_temp ON designdb.compounds;
CREATE TRIGGER trg_compound_catalogue_map_sync_compound_temp
    AFTER INSERT OR UPDATE OF compound_hash ON designdb.compounds
    FOR EACH ROW
    EXECUTE FUNCTION designdb.trg_compound_catalogue_map_from_compound_temp();
-- Will be deleted automatically until here
*/

DROP TRIGGER IF EXISTS trg_feature_updated_on ON designdb.features;
CREATE TRIGGER trg_feature_updated_on BEFORE UPDATE ON designdb.features FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_route_updated_on ON designdb.routes;
CREATE TRIGGER trg_route_updated_on BEFORE UPDATE ON designdb.routes FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_reaction_updated_on ON designdb.reactions;
CREATE TRIGGER trg_reaction_updated_on BEFORE UPDATE ON designdb.reactions FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_pose_updated_on ON designdb.poses;
CREATE TRIGGER trg_pose_updated_on BEFORE UPDATE ON designdb.poses FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_score_values_updated_on ON designdb.score_values;
CREATE TRIGGER trg_score_values_updated_on BEFORE UPDATE ON designdb.score_values FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_score_values_refresh_pivoted_mv ON designdb.score_values;
CREATE TRIGGER trg_score_values_refresh_pivoted_mv
  AFTER INSERT OR UPDATE OR DELETE ON designdb.score_values
  FOR EACH STATEMENT EXECUTE FUNCTION designdb.refresh_scores_per_pose_pivoted_mv();

DROP TRIGGER IF EXISTS trg_subsite_updated_on ON designdb.subsites;
CREATE TRIGGER trg_subsite_updated_on BEFORE UPDATE ON designdb.subsites FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_component_updated_on ON designdb.components;
CREATE TRIGGER trg_component_updated_on BEFORE UPDATE ON designdb.components FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_inspiration_updated_on ON designdb.inspirations;
CREATE TRIGGER trg_inspiration_updated_on BEFORE UPDATE ON designdb.inspirations FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_interaction_updated_on ON designdb.interactions;
CREATE TRIGGER trg_interaction_updated_on BEFORE UPDATE ON designdb.interactions FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_catalogue_updated_on ON designdb.catalogue_compounds;
CREATE TRIGGER trg_catalogue_updated_on BEFORE UPDATE ON designdb.catalogue_compounds FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_catalogue_price_updated_on ON designdb.catalogue_prices;
CREATE TRIGGER trg_catalogue_price_updated_on BEFORE UPDATE ON designdb.catalogue_prices FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

/*
-- Will be deleted automatically from db from here:
DROP TRIGGER IF EXISTS trg_compound_catalogue_map_sync_catalogue ON designdb.catalogue_compounds;
DROP TRIGGER IF EXISTS trg_compound_catalogue_map_sync_catalogue_temp ON designdb.catalogue_compounds;
CREATE TRIGGER trg_compound_catalogue_map_sync_catalogue_temp
    AFTER INSERT OR UPDATE OF catalogue_hash ON designdb.catalogue_compounds
    FOR EACH ROW
    EXECUTE FUNCTION designdb.trg_compound_catalogue_map_from_catalogue_temp();

DROP TRIGGER IF EXISTS trg_compound_catalogue_map_sync_catalogue_price ON designdb.catalogue_prices;
DROP TRIGGER IF EXISTS trg_compound_catalogue_map_sync_catalogue_price_temp ON designdb.catalogue_prices;
CREATE TRIGGER trg_compound_catalogue_map_sync_catalogue_price_temp
    AFTER INSERT OR UPDATE OF catalogue_id ON designdb.catalogue_prices
    FOR EACH ROW
    EXECUTE FUNCTION designdb.trg_compound_catalogue_map_from_catalogue_price_temp();

DROP TRIGGER IF EXISTS trg_compound_catalogue_map_refresh_denorm_from_catalogue_price ON designdb.catalogue_prices;
DROP TRIGGER IF EXISTS trg_compound_catalogue_map_refresh_denorm_from_catalogue_price_temp ON designdb.catalogue_prices;
CREATE TRIGGER trg_compound_catalogue_map_refresh_denorm_from_catalogue_price_temp
    AFTER UPDATE OF supplier, amount, price, lead_time ON designdb.catalogue_prices
    FOR EACH ROW
    EXECUTE FUNCTION designdb.trg_compound_catalogue_map_refresh_denorm_from_catalogue_price_temp();

DROP TRIGGER IF EXISTS trg_compound_catalogue_map_refresh_denorm_from_catalogue_compound ON designdb.catalogue_compounds;
DROP TRIGGER IF EXISTS trg_compound_catalogue_map_refresh_denorm_from_catalogue_compound_temp ON designdb.catalogue_compounds;
CREATE TRIGGER trg_compound_catalogue_map_refresh_denorm_from_catalogue_compound_temp
    AFTER UPDATE OF catalogue_inchikey ON designdb.catalogue_compounds
    FOR EACH ROW
    EXECUTE FUNCTION designdb.trg_compound_catalogue_map_refresh_denorm_from_catalogue_compound_temp();
-- Will be deleted automatically until here
*/

DROP TRIGGER IF EXISTS trg_compound_catalogue_map_updated_on ON designdb.compound_catalogue_map;
CREATE TRIGGER trg_compound_catalogue_map_updated_on BEFORE UPDATE ON designdb.compound_catalogue_map FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

-- =========================================================
-- AUDIT TRIGGERS
-- =========================================================

DROP TRIGGER IF EXISTS trg_catalogue_compounds_event_audit ON designdb.catalogue_compounds;
CREATE TRIGGER trg_catalogue_compounds_event_audit
    AFTER INSERT OR UPDATE OR DELETE ON designdb.catalogue_compounds
    FOR EACH ROW
    EXECUTE FUNCTION designdb.event_audit_trigger(
        'designdb.catalogue_compounds_event_audit',
        'id',
        'created_on,updated_on',
        ''
    );

DROP TRIGGER IF EXISTS trg_catalogue_prices_event_audit ON designdb.catalogue_prices;
CREATE TRIGGER trg_catalogue_prices_event_audit
    AFTER INSERT OR UPDATE OR DELETE ON designdb.catalogue_prices
    FOR EACH ROW
    EXECUTE FUNCTION designdb.event_audit_trigger(
        'designdb.catalogue_prices_event_audit',
        'id',
        'created_on,updated_on',
        ''
    );

DROP TRIGGER IF EXISTS trg_pose_tags_event_audit ON designdb.pose_tags;
CREATE TRIGGER trg_pose_tags_event_audit
    AFTER INSERT OR UPDATE OR DELETE ON designdb.pose_tags
    FOR EACH ROW
    EXECUTE FUNCTION designdb.event_audit_trigger(
        'designdb.pose_tags_event_audit',
        'id',
        'created_on,updated_on',
        ''
    );

DROP TRIGGER IF EXISTS trg_compound_tags_event_audit ON designdb.compound_tags;
CREATE TRIGGER trg_compound_tags_event_audit
    AFTER INSERT OR UPDATE OR DELETE ON designdb.compound_tags
    FOR EACH ROW
    EXECUTE FUNCTION designdb.event_audit_trigger(
        'designdb.compound_tags_event_audit',
        'id',
        'created_on,updated_on',
        ''
    );

DROP TRIGGER IF EXISTS trg_pose_methods_event_audit ON designdb.pose_methods;
CREATE TRIGGER trg_pose_methods_event_audit
    AFTER INSERT OR UPDATE OR DELETE ON designdb.pose_methods
    FOR EACH ROW
    EXECUTE FUNCTION designdb.event_audit_trigger(
        'designdb.pose_methods_event_audit',
        'id',
        'created_on,updated_on',
        ''
    );

DROP TRIGGER IF EXISTS trg_enumeration_methods_event_audit ON designdb.enumeration_methods;
CREATE TRIGGER trg_enumeration_methods_event_audit
    AFTER INSERT OR UPDATE OR DELETE ON designdb.enumeration_methods
    FOR EACH ROW
    EXECUTE FUNCTION designdb.event_audit_trigger(
        'designdb.enumeration_methods_event_audit',
        'id',
        'created_on,updated_on',
        ''
    );

DROP TRIGGER IF EXISTS trg_scoring_methods_event_audit ON designdb.scoring_methods;
CREATE TRIGGER trg_scoring_methods_event_audit
    AFTER INSERT OR UPDATE OR DELETE ON designdb.scoring_methods
    FOR EACH ROW
    EXECUTE FUNCTION designdb.event_audit_trigger(
        'designdb.scoring_methods_event_audit',
        'id',
        'created_on,updated_on',
        ''
    );

DROP TRIGGER IF EXISTS trg_reactant_updated_on ON designdb.reactants;
CREATE TRIGGER trg_reactant_updated_on BEFORE UPDATE ON designdb.reactants FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_scaffold_updated_on ON designdb.scaffolds;
CREATE TRIGGER trg_scaffold_updated_on BEFORE UPDATE ON designdb.scaffolds FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_subsite_tag_updated_on ON designdb.subsite_tags;
CREATE TRIGGER trg_subsite_tag_updated_on BEFORE UPDATE ON designdb.subsite_tags FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

-- Removed due to replaced tables
-- DROP TRIGGER IF EXISTS trg_tag_updated_on ON designdb.tags;
-- CREATE TRIGGER trg_tag_updated_on BEFORE UPDATE ON designdb.tags FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_pose_tag_updated_on ON designdb.pose_tags;
CREATE TRIGGER trg_pose_tag_updated_on BEFORE UPDATE ON designdb.pose_tags FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_compound_tag_updated_on ON designdb.compound_tags;
CREATE TRIGGER trg_compound_tag_updated_on BEFORE UPDATE ON designdb.compound_tags FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_has_pose_tag_updated_on ON designdb.has_pose_tags;
CREATE TRIGGER trg_has_pose_tag_updated_on BEFORE UPDATE ON designdb.has_pose_tags FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_has_pose_methods_updated_on ON designdb.has_pose_methods;
CREATE TRIGGER trg_has_pose_methods_updated_on BEFORE UPDATE ON designdb.has_pose_methods FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_has_compound_tag_updated_on ON designdb.has_compound_tags;
CREATE TRIGGER trg_has_compound_tag_updated_on BEFORE UPDATE ON designdb.has_compound_tags FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

DROP TRIGGER IF EXISTS trg_has_enumeration_methods_updated_on ON designdb.has_enumeration_methods;
CREATE TRIGGER trg_has_enumeration_methods_updated_on BEFORE UPDATE ON designdb.has_enumeration_methods FOR EACH ROW EXECUTE FUNCTION designdb.update_updated_on();

-- Pivoted materialized view once at schema load, this will be mapped in Scarab to do the filtering based on any type of scores/methods
SELECT designdb.create_scores_per_pose_pivoted_mv();

-- =============================================================================
-- HIPPO → production revert (manual step on live DB): strip leading "-- " from the
-- DROP/ALTER lines below and run in order; then remove the "/*" and "*/" around the
-- production block so it can be manually executed.
--
-- DROP TRIGGER IF EXISTS trg_compound_catalogue_map_refresh_denorm_from_catalogue_compound_temp ON designdb.catalogue_compounds;
-- DROP TRIGGER IF EXISTS trg_compound_catalogue_map_sync_catalogue_temp ON designdb.catalogue_compounds;
-- DROP TRIGGER IF EXISTS trg_compound_catalogue_map_sync_compound_temp ON designdb.compounds;
-- DROP TRIGGER IF EXISTS trg_compound_catalogue_map_refresh_denorm_from_catalogue_price_temp ON designdb.catalogue_prices;
-- DROP TRIGGER IF EXISTS trg_compound_catalogue_map_sync_catalogue_price_temp ON designdb.catalogue_prices;
--
-- DROP FUNCTION IF EXISTS designdb.trg_compound_catalogue_map_refresh_denorm_from_catalogue_compound_temp();
-- DROP FUNCTION IF EXISTS designdb.trg_compound_catalogue_map_from_catalogue_temp();
-- DROP FUNCTION IF EXISTS designdb.trg_compound_catalogue_map_refresh_denorm_from_catalogue_price_temp();
-- DROP FUNCTION IF EXISTS designdb.trg_compound_catalogue_map_from_catalogue_price_temp();
-- DROP FUNCTION IF EXISTS designdb.trg_compound_catalogue_map_from_compound_temp();
--
-- ALTER TABLE designdb.compound_catalogue_map DROP COLUMN IF EXISTS catalogue_inchikey;
-- ALTER TABLE designdb.compound_catalogue_map DROP COLUMN IF EXISTS supplier;
-- ALTER TABLE designdb.compound_catalogue_map DROP COLUMN IF EXISTS amount;
-- ALTER TABLE designdb.compound_catalogue_map DROP COLUMN IF EXISTS price;
-- ALTER TABLE designdb.compound_catalogue_map DROP COLUMN IF EXISTS lead_time;
-- =============================================================================

-- === PRODUCTION: compound_catalogue_map (INSERT only compound_id, catalogue_price_id, match_hash) ===

CREATE OR REPLACE FUNCTION designdb.trg_compound_catalogue_map_from_compound()
RETURNS trigger
LANGUAGE plpgsql
AS $$
BEGIN
    IF TG_OP = 'UPDATE' THEN
        IF OLD.compound_hash IS NOT DISTINCT FROM NEW.compound_hash THEN
            RETURN NEW;
        END IF;
        DELETE FROM designdb.compound_catalogue_map WHERE compound_id = NEW.id;
    END IF;

    INSERT INTO designdb.compound_catalogue_map (compound_id, catalogue_price_id, match_hash)
    SELECT NEW.id, cp.id, NEW.compound_hash
    FROM designdb.catalogue_prices cp
    JOIN designdb.catalogue_compounds cat ON cat.id = cp.catalogue_id
    WHERE cat.catalogue_hash = NEW.compound_hash
    ON CONFLICT (compound_id, catalogue_price_id) DO NOTHING;

    RETURN NEW;
END;
$$;

CREATE OR REPLACE FUNCTION designdb.trg_compound_catalogue_map_from_catalogue()
RETURNS trigger
LANGUAGE plpgsql
AS $$
BEGIN
    IF TG_OP = 'UPDATE' THEN
        IF OLD.catalogue_hash IS NOT DISTINCT FROM NEW.catalogue_hash THEN
            RETURN NEW;
        END IF;
        DELETE FROM designdb.compound_catalogue_map
        WHERE catalogue_price_id IN (
            SELECT id FROM designdb.catalogue_prices WHERE catalogue_id = NEW.id
        );
    END IF;

    INSERT INTO designdb.compound_catalogue_map (compound_id, catalogue_price_id, match_hash)
    SELECT c.id, cp.id, c.compound_hash
    FROM designdb.compounds c
    CROSS JOIN designdb.catalogue_prices cp
    WHERE cp.catalogue_id = NEW.id AND c.compound_hash = NEW.catalogue_hash
    ON CONFLICT (compound_id, catalogue_price_id) DO NOTHING;

    RETURN NEW;
END;
$$;

CREATE OR REPLACE FUNCTION designdb.trg_compound_catalogue_map_from_catalogue_price()
RETURNS trigger
LANGUAGE plpgsql
AS $$
BEGIN
    IF TG_OP = 'UPDATE' AND OLD.catalogue_id IS NOT DISTINCT FROM NEW.catalogue_id THEN
        RETURN NEW;
    END IF;
    IF TG_OP = 'UPDATE' THEN
        DELETE FROM designdb.compound_catalogue_map WHERE catalogue_price_id = NEW.id;
    END IF;

    INSERT INTO designdb.compound_catalogue_map (compound_id, catalogue_price_id, match_hash)
    SELECT c.id, NEW.id, c.compound_hash
    FROM designdb.compounds c
    JOIN designdb.catalogue_compounds cat ON cat.id = NEW.catalogue_id
    WHERE c.compound_hash = cat.catalogue_hash
    ON CONFLICT (compound_id, catalogue_price_id) DO NOTHING;

    RETURN NEW;
END;
$$;

DROP TRIGGER IF EXISTS trg_compound_catalogue_map_sync_compound ON designdb.compounds;
CREATE TRIGGER trg_compound_catalogue_map_sync_compound
    AFTER INSERT OR UPDATE OF compound_hash ON designdb.compounds
    FOR EACH ROW
    EXECUTE FUNCTION designdb.trg_compound_catalogue_map_from_compound();

DROP TRIGGER IF EXISTS trg_compound_catalogue_map_sync_catalogue ON designdb.catalogue_compounds;
CREATE TRIGGER trg_compound_catalogue_map_sync_catalogue
    AFTER INSERT OR UPDATE OF catalogue_hash ON designdb.catalogue_compounds
    FOR EACH ROW
    EXECUTE FUNCTION designdb.trg_compound_catalogue_map_from_catalogue();

DROP TRIGGER IF EXISTS trg_compound_catalogue_map_sync_catalogue_price ON designdb.catalogue_prices;
CREATE TRIGGER trg_compound_catalogue_map_sync_catalogue_price
    AFTER INSERT OR UPDATE OF catalogue_id ON designdb.catalogue_prices
    FOR EACH ROW
    EXECUTE FUNCTION designdb.trg_compound_catalogue_map_from_catalogue_price();
