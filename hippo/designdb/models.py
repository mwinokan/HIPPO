from pathlib import Path

import mrich

# from django.db.models import indexes
from django.conf import settings
from django.db import models
from django.db.models import Q
from django.utils import timezone
from rdkit import Chem

from .managers import CompoundManager

_MANAGE_MODELS = settings.MANAGE_MODELS


# Custom field type for text fields that store json. Once switchint to
# postgres, replace
class JSONTextField(models.TextField):
    def from_db_value(self, value, expression, connection):
        import json

        return json.loads(value) if value else {}

    def get_prep_value(self, value):
        import json

        if isinstance(value, dict):
            return json.dumps(value)
        return value


class RDKitMolField(models.TextField):
    """
    Stores RDKit molecules as MolBlock text (SDF format) in DB,
    but returns RDKit Mol objects in Python.
    """

    description = 'RDKit molecule stored as MolBlock text'

    # -------------------------
    # DB → Python (read path)
    # -------------------------
    def from_db_value(self, value, expression, connection):
        if not value:
            return None
        return Chem.MolFromMolBlock(value)

    # -------------------------
    # Python → DB (write path)
    # -------------------------
    def get_prep_value(self, value):
        if value is None:
            return None

        # Already serialized
        if isinstance(value, str):
            return value

        # RDKit Mol → MolBlock
        if isinstance(value, Chem.Mol):
            return Chem.MolToMolBlock(value)

        raise TypeError(
            'RDKitMolField only accepts RDKit Mol or MolBlock string, '
            f'got {type(value)}'
        )

    def deconstruct(self):
        name, path, args, kwargs = super().deconstruct()
        return name, path, args, kwargs


if settings.MANAGE_MODELS:
    # sqlite3: the RDKit cartridge field types are unavailable, fall back to the
    # plain-text RDKitMolField shim defined above.
    MolField = RDKitMolField
else:
    from django_rdkit.models import MolField


class BaseModel(models.Model):
    created_on = models.DateTimeField(null=True, blank=True, default=timezone.now)
    updated_on = models.DateTimeField(null=True, blank=True, default=timezone.now)

    class Meta:
        abstract = True
        managed = _MANAGE_MODELS
        app_label = 'designdb'
        default_related_name = '%(class)ss'


class Project(BaseModel):
    project_name = models.TextField(null=False, unique=True)
    open_to_public = models.BooleanField(default=False)

    class Meta(BaseModel.Meta):
        db_table = 'projects'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'project_name',
                ],
                name='uc_project',
            ),
        ]

    def __str__(self) -> str:
        return f'{self.project_name}'


class TargetModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    external_target_id = models.BigIntegerField(null=True, blank=True)
    target_name = models.TextField()
    target_metadata = models.TextField(null=True, blank=True)
    project = models.ForeignKey(
        Project,
        on_delete=models.RESTRICT,
        db_column='project_id',
    )

    class Meta(BaseModel.Meta):
        db_table = 'targets'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'target_name',
                ],
                name='uc_target',
            ),
        ]
        indexes = [
            models.Index(fields=['target_name'], name='idx_target_name'),
            models.Index(fields=['created_on'], name='idx_target_created'),
        ]


class CompoundModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    compound_inchikey = models.TextField(null=True, blank=True)
    compound_alias = models.TextField(null=True, blank=True)
    compound_smiles = models.TextField(null=True, blank=True)
    compound_hash = models.TextField(null=False, blank=True, default='a')

    base_compound = models.ForeignKey(
        'self',
        null=True,
        blank=True,
        on_delete=models.SET_NULL,
        db_column='base_compound_id',
        related_name='+',  # add if needed
    )

    # compound_mol = models.TextField(null=True, blank=True)
    # compound_pattern_bfp = models.TextField(null=True, blank=True)
    # compound_morgan_bfp = models.TextField(null=True, blank=True)
    # compound_mol = MolField(null=True)
    compound_mol = models.TextField(null=True, blank=True)
    # compound_pattern_bfp = BfpField(null=True)
    # compound_morgan_bfp = BfpField(null=True)
    compound_pattern_bfp = models.BinaryField(max_length=2048, null=True)
    compound_morgan_bfp = models.BinaryField(max_length=2048, null=True)

    compound_metadata = models.TextField(null=True, blank=True)
    note = models.TextField(null=True, blank=True)
    rdkit_version = models.TextField(null=True, blank=True)
    inchi_version = models.TextField(null=True, blank=True)

    tags = models.ManyToManyField(
        'CompoundTagModel',
        through='CompoundTagJunctionModel',
        related_name='compounds',
    )

    enumeration_methods = models.ManyToManyField(
        'EnumerationMethodModel',
        through='CompoundEnumerationMethodJunctionModel',
        related_name='compounds',
    )

    # unlike others, this wasn't clearly defined as m2m. may not want
    # to keep it
    scaffolds = models.ManyToManyField(
        'self',
        through='ScaffoldModel',
    )

    objects = models.Manager()
    compound_filter = CompoundManager()

    class Meta(BaseModel.Meta):
        db_table = 'compounds'
        constraints = [
            # I believe there were supposed to be changes to these
            # models.UniqueConstraint(
            #     fields=[
            #         'compound_alias',
            #     ],
            #     name='uc_compound_alias',
            # ),
            models.UniqueConstraint(
                fields=[
                    'compound_inchikey',
                ],
                name='uc_compound_inchikey',
            ),
            # tautomers mess this up
            # models.UniqueConstraint(
            #     fields=[
            #         'compound_smiles',
            #     ],
            #     name='uc_compound_smiles',
            # ),
        ]
        indexes = [
            # models.Index(fields=['base_compound'], name='idx_base_compound_id'),
            models.Index(fields=['compound_inchikey'], name='idx_compound_inchikey'),
            # models.Index(fields=['compound_smiles'], name='idx_compound_smiles'),
            models.Index(fields=['created_on'], name='idx_compound_created'),
        ]


class SubsiteModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    target = models.ForeignKey(
        TargetModel,
        on_delete=models.RESTRICT,
        db_column='target_id',
    )

    subsite_name = models.TextField()
    subsite_metadata = models.TextField(null=True, blank=True)

    class Meta(BaseModel.Meta):
        db_table = 'subsites'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'target',
                    'subsite_name',
                ],
                name='uc_subsite',
            ),
        ]
        indexes = [
            models.Index(fields=['target'], name='idx_subsite_target_id'),
            models.Index(fields=['created_on'], name='idx_subsite_created'),
        ]


class PoseModel(BaseModel):
    id = models.BigAutoField(primary_key=True)

    pose_inchikey = models.TextField(null=True, blank=True)
    pose_alias = models.TextField(null=True, blank=True)
    pose_smiles = models.TextField(null=True, blank=True)

    pose_reference = models.IntegerField(null=True, blank=True)
    protein_link = models.TextField(null=True, blank=True)

    compound = models.ForeignKey(
        CompoundModel,
        on_delete=models.RESTRICT,
        db_column='compound_id',
    )

    target = models.ForeignKey(
        TargetModel,
        on_delete=models.RESTRICT,
        db_column='target_id',
    )

    # pose_mol = models.TextField(null=True, blank=True)
    pose_mol = MolField(null=True)
    # this is integer in the db.. pretty sure this cannot be the case?
    pose_fingerprint = models.IntegerField(null=True, blank=True)

    # dicts dumped into that field, change to JSON?
    # pose_metadata = models.TextField(null=True, blank=True)
    pose_metadata = JSONTextField(null=True, blank=True)
    # pose_metadata = models.JSONField(null=True, blank=True)
    note = models.TextField(null=True, blank=True)

    rdkit_version = models.TextField(null=True, blank=True)
    inchi_version = models.TextField(null=True, blank=True)

    methods = models.ManyToManyField(
        'PoseMethodModel',
        through='PoseMethodJunctionModel',
        related_name='poses',
    )
    tags = models.ManyToManyField(
        'PoseTagModel',
        through='PoseTagJunctionModel',
        related_name='poses',
    )
    # unlike others, this wasn't clearly defined as m2m. may not want
    # to keep it
    # through_fields: the calling pose is the derivative, the added pose is the original
    inspirations = models.ManyToManyField(
        'self',
        through='InspirationModel',
        through_fields=('derivative_pose', 'original_pose'),
        symmetrical=False,
    )

    subsites = models.ManyToManyField(
        SubsiteModel,
        through='SubsiteTagModel',
    )

    class Meta(BaseModel.Meta):
        db_table = 'poses'
        # There are no constraints here, but they need to be unique,
        # verified in code (rdkit.align_pose coords from
        # file). Investigate adding coords to db and doing the search
        # there

        # although.. would alias-target combo work?
        indexes = [
            models.Index(fields=['compound'], name='idx_pose_compound_id'),
            models.Index(fields=['target'], name='idx_pose_target_id'),
            models.Index(fields=['protein_link'], name='idx_protein_link'),
            models.Index(fields=['created_on'], name='idx_pose_created'),
        ]

    @property
    def mol_path(self) -> Path | None:
        """Get Path to molecule file"""
        path = Path(self.protein_link)
        if path.name.endswith('.pdb'):
            mol_path = path.parent / path.name.replace('_hippo.pdb', '.pdb').replace(
                '.pdb', '_ligand.mol'
            )
            if not mol_path.exists():
                mol_path = path.parent / path.name.replace(
                    '_hippo.pdb', '.pdb'
                ).replace('.pdb', '_ligand.sdf')
                if not mol_path.exists():
                    mrich.error('Could not find ligand mol/sdf:', mol_path)
                    return None
            return mol_path
        elif path.name.endswith('.mol'):
            return path
        else:
            raise NotImplementedError

    @property
    def apo_path(self) -> Path | None:
        """Get path to the apo (de-liganded, desolvated) protein file"""
        path = Path(self.protein_link)
        if path.name.endswith('.pdb'):
            stem = path.name.replace('_hippo.pdb', '.pdb')
            # current Fragalysis protein-file naming
            apo_path = path.parent / stem.replace('.pdb', '_delig-desolv.pdb')
            if apo_path.exists():
                return apo_path
            # DEPRECATED(apo-naming): pre-'delig' Fragalysis naming, remove once
            # all data uses 'delig'
            legacy_path = path.parent / stem.replace('.pdb', '_apo-desolv.pdb')
            if legacy_path.exists():
                return legacy_path
            return None
        else:
            raise NotImplementedError


class SubsiteTagModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    pose = models.ForeignKey(
        PoseModel,
        on_delete=models.RESTRICT,
        db_column='pose_id',
    )
    subsite = models.ForeignKey(
        SubsiteModel,
        on_delete=models.RESTRICT,
        db_column='subsite_id',
    )

    subsite_tag_metadata = models.TextField(null=True, blank=True)

    class Meta(BaseModel.Meta):
        db_table = 'subsite_tags'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'subsite',
                    'pose',
                ],
                name='uc_subsite_tag',
            ),
        ]
        indexes = [
            models.Index(fields=['subsite'], name='idx_subsite_tag_subsite_id'),
            models.Index(fields=['pose'], name='idx_subsite_tag_pose_id'),
            models.Index(fields=['created_on'], name='idx_subsite_tag_created'),
        ]


class PoseMethodModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    pose_method_name = models.TextField(null=True, blank=True)
    pose_method_description = models.TextField(null=True, blank=True)
    pose_method_version = models.TextField(null=True, blank=True)
    pose_method_organization = models.TextField(null=True, blank=True)
    pose_method_link = models.TextField(null=True, blank=True)
    pose_method_note = models.TextField(null=True, blank=True)

    class Meta(BaseModel.Meta):
        db_table = 'pose_methods'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'pose_method_name',
                    'pose_method_version',
                ],
                name='uc_pose_method',
                nulls_distinct=False,
            )
        ]
        indexes = [
            models.Index(fields=['pose_method_name'], name='idx_pose_method_name'),
            models.Index(fields=['created_on'], name='idx_pose_method_created'),
        ]


class PoseMethodJunctionModel(BaseModel):
    pk = models.CompositePrimaryKey('pose_id', 'pose_method_id')
    pose = models.ForeignKey(
        'PoseModel',
        on_delete=models.CASCADE,
        db_column='pose_id',
    )

    pose_method = models.ForeignKey(
        'PoseMethodModel',
        on_delete=models.CASCADE,
        db_column='pose_method_id',
    )

    class Meta(BaseModel.Meta):
        db_table = 'has_pose_methods'
        indexes = [
            models.Index(
                fields=['pose_method'], name='idx_has_pose_methods_pose_method_id'
            ),
            models.Index(
                fields=['created_on'], name='idx_idx_has_pose_methods_created'
            ),
        ]


class PoseTagModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    pose_tag_name = models.TextField()
    pose_tag_description = models.TextField(null=True, blank=True)
    pose_tag_note = models.TextField(null=True, blank=True)

    class Meta(BaseModel.Meta):
        db_table = 'pose_tags'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'pose_tag_name',
                ],
                name='uc_pose_tag',
            )
        ]
        indexes = [
            models.Index(fields=['created_on'], name='idx_pose_tag_created'),
        ]


class PoseTagJunctionModel(BaseModel):
    pk = models.CompositePrimaryKey('pose_id', 'pose_tag_id')
    pose = models.ForeignKey(
        PoseModel,
        on_delete=models.CASCADE,
        db_column='pose_id',
    )
    pose_tag = models.ForeignKey(
        PoseTagModel,
        on_delete=models.CASCADE,
        db_column='pose_tag_id',
    )

    class Meta(BaseModel.Meta):
        db_table = 'has_pose_tags'
        indexes = [
            models.Index(fields=['pose_tag'], name='idx_has_pose_tag_pose_tag_id'),
            models.Index(fields=['created_on'], name='idx_has_pose_tag_created'),
        ]


# this was missing.. is this a m2m table as well? really looks like it
class InspirationModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    # original behaviour described in schema was SET_NULL but I don't
    # see how that makes sense. if either original or derivative is
    # deleted, you'll have orphaned entries
    original_pose = models.ForeignKey(
        PoseModel,
        # on_delete=models.SET_NULL,
        on_delete=models.CASCADE,
        db_column='original_pose_id',
        related_name='+',
    )
    derivative_pose = models.ForeignKey(
        PoseModel,
        # on_delete=models.SET_NULL,
        on_delete=models.CASCADE,
        db_column='derivative_pose_id',
        related_name='+',
    )

    class Meta(BaseModel.Meta):
        db_table = 'inspirations'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'original_pose',
                    'derivative_pose',
                ],
                name='uc_inspiration',
            )
        ]
        indexes = [
            models.Index(
                fields=['original_pose'], name='idx_inspiration_original_pose_id'
            ),
            models.Index(
                fields=['derivative_pose'], name='idx_inspiration_derivative_pose_id'
            ),
            models.Index(fields=['created_on'], name='idx_inspiration_created'),
        ]


class FeatureModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    feature_family = models.TextField(null=True, blank=True)
    target = models.ForeignKey(
        TargetModel,
        on_delete=models.RESTRICT,
        db_column='target_id',
    )

    feature_chain_name = models.TextField(null=True, blank=True)
    feature_residue_name = models.TextField(null=True, blank=True)
    feature_residue_number = models.IntegerField(null=True, blank=True)
    feature_atom_name = models.TextField(null=True, blank=True)

    class Meta(BaseModel.Meta):
        db_table = 'features'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'feature_family',
                    'target',
                    'feature_chain_name',
                    'feature_residue_name',
                    'feature_residue_number',
                    'feature_atom_name',
                ],
                name='uc_feature',
            )
        ]
        indexes = [
            models.Index(fields=['target'], name='idx_feature_target_id'),
            models.Index(fields=['created_on'], name='idx_feature_created'),
        ]


class InteractionModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    feature = models.ForeignKey(
        FeatureModel,
        on_delete=models.RESTRICT,
        db_column='feature_id',
    )
    pose = models.ForeignKey(
        PoseModel,
        on_delete=models.RESTRICT,
        db_column='pose_id',
    )

    interaction_type = models.TextField()
    interaction_family = models.TextField()
    interaction_atom_id = models.TextField()

    # could these be vectors?
    interaction_prot_coord = models.TextField()
    interaction_lig_coord = models.TextField()

    interaction_distance = models.FloatField()
    interaction_angle = models.FloatField(null=True, blank=True)
    interaction_energy = models.FloatField(null=True, blank=True)

    class Meta(BaseModel.Meta):
        db_table = 'interactions'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'feature',
                    'pose',
                    'interaction_type',
                    'interaction_family',
                    'interaction_atom_id',
                ],
                name='uc_interaction',
            )
        ]
        indexes = [
            models.Index(fields=['feature_id'], name='idx_interaction_feature_id'),
            models.Index(fields=['pose'], name='idx_interaction_pose_id'),
            models.Index(fields=['created_on'], name='idx_interaction_created'),
        ]


class CompoundTagModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    compound_tag_name = models.TextField()
    compound_tag_description = models.TextField(null=True, blank=True)
    compound_tag_note = models.TextField(null=True, blank=True)

    class Meta(BaseModel.Meta):
        db_table = 'compound_tags'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'compound_tag_name',
                ],
                name='uc_compound_tag_name',
            )
        ]
        indexes = [
            models.Index(fields=['created_on'], name='idx_compound_tag_created'),
        ]


class CompoundTagJunctionModel(BaseModel):
    pk = models.CompositePrimaryKey('compound_id', 'compound_tag_id')
    compound = models.ForeignKey(
        CompoundModel,
        on_delete=models.CASCADE,
        db_column='compound_id',
    )
    compound_tag = models.ForeignKey(
        CompoundTagModel,
        on_delete=models.CASCADE,
        db_column='compound_tag_id',
    )

    class Meta(BaseModel.Meta):
        db_table = 'has_compound_tags'
        indexes = [
            models.Index(
                fields=['compound_tag'], name='idx_has_compound_tag_compound_tag_id'
            ),
            models.Index(fields=['created_on'], name='idx_has_compound_tag_created'),
        ]


class EnumerationMethodModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    enum_name = models.TextField(null=True, blank=True)
    enum_description = models.TextField(null=True, blank=True)
    enum_version = models.TextField(null=True, blank=True)
    enum_organization = models.TextField(null=True, blank=True)
    enum_link = models.TextField(null=True, blank=True)
    enum_note = models.TextField(null=True, blank=True)

    class Meta(BaseModel.Meta):
        db_table = 'enumeration_methods'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'enum_name',
                    'enum_version',
                ],
                name='uc_enumeration_method',
                nulls_distinct=False,
            )
        ]
        indexes = [
            models.Index(fields=['enum_name'], name='idx_enumeration_method_name'),
            models.Index(fields=['created_on'], name='idx_enumeration_method_created'),
        ]


class CompoundEnumerationMethodJunctionModel(BaseModel):
    pk = models.CompositePrimaryKey('compound_id', 'enumeration_method_id')
    compound = models.ForeignKey(
        CompoundModel,
        on_delete=models.CASCADE,
        db_column='compound_id',
    )
    enumeration_method = models.ForeignKey(
        EnumerationMethodModel,
        on_delete=models.CASCADE,
        db_column='enumeration_method_id',
    )

    class Meta(BaseModel.Meta):
        db_table = 'has_enumeration_methods'
        indexes = [
            models.Index(
                fields=['enumeration_method'],
                name='idx_has_enumeration_methods_enumeration_method_id',
            ),
            models.Index(
                fields=['created_on'], name='idx_has_enumeration_methods_created'
            ),
        ]


class ScoringMethodModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    method_name = models.TextField(null=True, blank=True)
    method_description = models.TextField(null=True, blank=True)
    method_version = models.TextField(null=True, blank=True)
    method_organization = models.TextField(null=True, blank=True)
    method_link = models.TextField(null=True, blank=True)
    note = models.TextField(null=True, blank=True)

    class Meta(BaseModel.Meta):
        db_table = 'scoring_methods'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'method_name',
                    'method_version',
                ],
                name='uc_scoring_method',
                nulls_distinct=False,
            )
        ]
        indexes = [
            models.Index(fields=['method_name'], name='idx_scoring_method_name'),
            models.Index(fields=['created_on'], name='idx_scoring_method_created'),
        ]


class ScoreValueModel(BaseModel):
    pk = models.CompositePrimaryKey('pose_id', 'compound_id', 'scoring_method_id')
    pose = models.ForeignKey(
        PoseModel,
        on_delete=models.RESTRICT,
        db_column='pose_id',
        related_name='scores',
    )

    compound = models.ForeignKey(
        CompoundModel,
        on_delete=models.RESTRICT,
        db_column='compound_id',
        related_name='scores',
    )

    scoring_method = models.ForeignKey(
        ScoringMethodModel,
        on_delete=models.RESTRICT,
        db_column='scoring_method_id',
        related_name='scores',
    )

    score = models.JSONField()

    class Meta(BaseModel.Meta):
        db_table = 'score_values'
        indexes = [
            models.Index(fields=['pose'], name='idx_score_values_pose_id'),
            models.Index(fields=['compound'], name='idx_score_values_compound_id'),
            models.Index(
                fields=['scoring_method'], name='idx_score_values_scoring_method_id'
            ),
            models.Index(fields=['created_on'], name='idx_score_values_created'),
        ]


class ReactionModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    reaction_type = models.TextField(null=True, blank=True)
    product_compound = models.ForeignKey(
        CompoundModel,
        on_delete=models.RESTRICT,
        db_column='product_compound_id',
    )

    reaction_product_yield = models.FloatField(null=True, blank=True)
    reaction_metadata = models.TextField(null=True, blank=True)

    class Meta(BaseModel.Meta):
        db_table = 'reactions'
        indexes = [
            models.Index(
                fields=['product_compound'], name='idx_reaction_product_compound_id'
            ),
            models.Index(fields=['created_on'], name='idx_reaction_created'),
        ]


class ReactantModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    reactant_amount = models.FloatField(null=True, blank=True)
    reaction = models.ForeignKey(
        ReactionModel,
        on_delete=models.CASCADE,
        db_column='reaction_id',
        related_name='reactants',
    )
    compound = models.ForeignKey(
        CompoundModel,
        on_delete=models.RESTRICT,
        db_column='compound_id',
    )

    class Meta(BaseModel.Meta):
        db_table = 'reactants'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'reaction',
                    'compound',
                ],
                name='uc_reactant',
            )
        ]
        indexes = [
            models.Index(fields=['reaction'], name='idx_reactant_reaction_id'),
            models.Index(fields=['compound'], name='idx_reactant_compound_id'),
            models.Index(fields=['created_on'], name='idx_reactant_created'),
        ]


class CatalogueCompoundModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    catalogue_smiles = models.TextField(null=False, blank=True)
    catalogue_inchikey = models.TextField(null=False, blank=True)
    catalogue_hash = models.TextField(null=False, blank=True)
    rdkit_version = models.TextField(null=True, blank=True)
    inchi_version = models.TextField(null=True, blank=True)

    class Meta(BaseModel.Meta):
        db_table = 'catalogue_compounds'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'catalogue_smiles',
                ],
                name='uq_catalogue_compounds_smiles',
            ),
            models.CheckConstraint(
                condition=Q(catalogue_hash__isnull=False) & Q(catalogue_hash__gt=''),
                name='ck_catalogue_compounds_hash_nonempty',
            ),
        ]


class CataloguePriceModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    catalogue_compound = models.ForeignKey(
        CatalogueCompoundModel,
        null=True,
        on_delete=models.CASCADE,
        db_column='catalogue_id',
    )
    vendor = models.TextField(null=False, blank=True)
    supplier = models.TextField(null=True, blank=True)
    supplier_id = models.TextField(null=False, blank=True)
    amount = models.FloatField(null=True, blank=True)
    price = models.FloatField(null=True, blank=True)
    currency = models.TextField(null=True, blank=True)
    purity = models.FloatField(null=True, blank=True)
    lead_time = models.IntegerField(null=True, blank=True)

    compounds = models.ManyToManyField(
        CompoundModel,
        through='CataloguePriceCompoundJunctionModel',
        related_name='prices',
    )

    class Meta(BaseModel.Meta):
        db_table = 'catalogue_prices'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'catalogue_compound',
                    'vendor',
                    'supplier',
                    'supplier_id',
                    'amount',
                ],
                name='uc_catalogue_price',
            )
        ]


class CataloguePriceCompoundJunctionModel(BaseModel):
    ipk = models.CompositePrimaryKey('compound_id', 'catalogue_price_id')
    catalogue_price = models.ForeignKey(
        CataloguePriceModel,
        on_delete=models.CASCADE,
        db_column='catalogue_price_id',
    )
    compound = models.ForeignKey(
        CompoundModel,
        on_delete=models.CASCADE,
        db_column='compound_id',
    )

    match_hash = models.TextField(null=False, blank=True)

    # Not needed, remove
    # catalogue_inchikey = models.TextField(null=False, blank=True)
    # supplier = models.TextField(null=True, blank=True)
    # amount = models.FloatField(null=True, blank=True)
    # price = models.FloatField(null=True, blank=True)
    # lead_time = models.IntegerField(null=True, blank=True)

    class Meta(BaseModel.Meta):
        db_table = 'compound_catalogue_map'
        constraints = [
            models.CheckConstraint(
                condition=Q(match_hash__isnull=False) & Q(match_hash__gt=''),
                name='ck_compound_catalogue_map_match_hash_nonempty',
            )
        ]


class ScaffoldModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    # same comment as with inspiratons. original schema says SET_NULL
    # but doesn't seem right
    base_compound = models.ForeignKey(
        CompoundModel,
        # on_delete=models.SET_NULL,
        on_delete=models.CASCADE,
        db_column='base_compound_id',
        related_name='scaffold_bases',
    )
    superstructure_compound = models.ForeignKey(
        CompoundModel,
        # on_delete=models.SET_NULL,
        on_delete=models.CASCADE,
        db_column='superstructure_compound_id',
        related_name='scaffold_superstructures',
    )

    class Meta(BaseModel.Meta):
        db_table = 'scaffolds'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'base_compound',
                    'superstructure_compound',
                ],
                name='uc_scaffold',
            )
        ]
        indexes = [
            models.Index(
                fields=['base_compound'], name='idx_scaffold_base_compound_id'
            ),
            models.Index(
                fields=['superstructure_compound'],
                name='idx_scaffold_superstructure_compound_id',
            ),
            models.Index(fields=['created_on'], name='idx_scaffold_created'),
        ]


class RouteModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    product_compound = models.ForeignKey(
        CompoundModel,
        on_delete=models.RESTRICT,
        db_column='product_compound_id',
    )

    class Meta(BaseModel.Meta):
        db_table = 'routes'
        indexes = [
            models.Index(
                fields=['product_compound'], name='idx_route_product_compound_id'
            ),
            models.Index(fields=['created_on'], name='idx_route_created'),
        ]


class ComponentModel(BaseModel):
    id = models.BigAutoField(primary_key=True)
    route = models.ForeignKey(
        RouteModel,
        on_delete=models.RESTRICT,
        db_column='route_id',
    )
    component_type = models.IntegerField(null=True, blank=True)
    component_ref = models.IntegerField(null=True, blank=True)
    component_amount = models.FloatField(null=True, blank=True)

    class Meta(BaseModel.Meta):
        db_table = 'components'
        constraints = [
            models.UniqueConstraint(
                fields=[
                    'route',
                    'component_ref',
                    'component_type',
                ],
                name='uc_component',
            )
        ]
        indexes = [
            models.Index(fields=['route'], name='idx_component_route_id'),
            models.Index(fields=['created_on'], name='idx_component_created'),
        ]


# what follows is audit tables, indexes, materialised views (none),
# views, functions and triggers. I'm not sure I need them here, will create if do.


# these functions available in db
# CREATE OR REPLACE FUNCTION designdb.mol_from_smiles(smiles TEXT) RETURNS rdkit.mol
#   LANGUAGE SQL AS $$ SELECT rdkit.mol_from_smiles(smiles::cstring); $$;

# CREATE OR REPLACE FUNCTION designdb.mol_to_smiles(m rdkit.mol) RETURNS text
#   LANGUAGE SQL AS $$ SELECT rdkit.mol_to_smiles(m); $$;

# CREATE OR REPLACE FUNCTION designdb.mol_to_inchikey(m rdkit.mol) RETURNS text
#   LANGUAGE SQL AS $$ SELECT rdkit.mol_inchikey(m); $$;
