import mout

from .fingerprint import Fingerprint
from .tset import TagSet
from rdkit import Chem


class Pose:

    def __init__(
        self, name, compound, pdb_path, site_index=None, chain=None, tags=None
    ):

        # blanks
        self._fingerprint = None
        self._pdb_entry = None
        self._mol = None
        self._site_index = None
        self._chain = None

        self._stereo_smiles = None  # smiles from syndirella
        self._smiles = None

        # arguments
        self._name = name
        self._compound = compound
        self._pdb_path = pdb_path

        self._placement_mRMSD = None
        self._placement_ddG = None

        if site_index is not None:
            self.site_index = site_index
        if chain is not None:
            self.chain = chain

        self._tags = TagSet(tags or [])

    ### FACTORIES

    @classmethod
    def from_bound_pdb(
        cls, compound, pdb_path, metadata=None, tags=None, site_index=None, chain=None
    ):
        """Usually from crystalstructure / Fragalysis LHS"""

        self = cls.__new__(cls)

        tags = TagSet(tags or [])

        if metadata:

            if metadata["site_name"]:
                tags.add(metadata["site_name"])

            if metadata["pdb_entry"]:
                self._pdb_entry = metadata["pdb_entry"]
                self.tags.add("InPDB")

            pose_name = metadata["crystal_name"]

            if site_index is None:
                site_index = pose_name[-2]

            if chain is None:
                chain = pose_name[-1]

            name = f"{site_index}{chain}"

        else:

            if site_index is not None and chain is not None:
                name = f"{site_index}{chain}"
            else:
                name = "??"

        self.__init__(name, compound, pdb_path, site_index, chain, tags)

        return self

    @classmethod
    def from_mol_and_reference(
        cls,
        compound,
        mol,
        protein_reference,
        metadata=None,
        tags=None,
        site_index=None,
        chain=None,
    ):
        """Usually a virtual hit from Syndirella"""

        self = cls.__new__(cls)

        tags = TagSet(tags or [])

        name = compound.name

        self.__init__(name, compound, protein_reference, site_index, chain, tags)

        self._mol = mol

        return self

    ### PROPERTIES

    @property
    def smiles(self):
        if not self._smiles:
            self._smiles = Chem.MolToSmiles(self.mol)
        return self._smiles

    @property
    def tags(self):
        return self._tags

    @property
    def compound(self):
        return self._compound

    @property
    def name(self):
        return self._name

    @property
    def longname(self):
        return f"{self.compound.name}_{self.name}"

    @property
    def site_index(self):
        assert self._site_index is not None
        return self._site_index

    @site_index.setter
    def site_index(self, i):
        assert i is not None
        # mout.warning(f'{self.name} changing site index! ({i})')
        self._site_index = int(i)

    @property
    def chain(self):
        return self._chain

    @chain.setter
    def chain(self, i):
        s = str(i)
        assert len(s) == 1
        self._chain = s

    @property
    def pdb_path(self):
        return self._pdb_path

    @property
    def fingerprint(self):
        # if self._fingerprint is None:
        # 	self._fingerprint = self.calculate_fingerprint()
        return self._fingerprint

    @property
    def bound_system(self):
        import molparse as mp

        sys = mp.parsePDB(self.pdb_path, dry=True, verbosity=0)
        return sys

    @property
    def mol(self):

        # mout.debug(f'Pose.mol({self.longname})')

        if self._mol is None:

            lig_residues = self.bound_system["rLIG"]

            if self.chain:
                lig_residues = [l for l in lig_residues if l.chain == self.chain]

            if len(lig_residues) and self.chain is None:
                self.chain = lig_residues[0].atoms[0].chain

            split_lig_residues = []
            for lig in lig_residues:
                split_lig_residues += lig.split_by_site()

            if self.site_index is not None:
                ligand_group = split_lig_residues[self.site_index]
            elif len(split_lig_residues) == 1:
                ligand_group = split_lig_residues[0]
            else:
                mout.error(
                    f"Bound PDB has multiple ligands and {self.longname}.site_index == {self.site_index}"
                )
                return self.compound.mol

            from molparse.rdkit import mol_from_pdb_block

            self._mol = mol_from_pdb_block(ligand_group.pdb_block)

        return self._mol

    @property
    def dict(self):
        d = dict(
            name=self.name,
            # smiles = self.smiles,
            longname=self.longname,
            # orig_smiles = self.orig_smiles,
            site_index=self.site_index,
            chain=self.chain,
            pdb_path=self.pdb_path,
            # base = self.base,
        )

        if hasattr(self, "_placement_mRMSD"):
            d["_placement_mRMSD"] = self._placement_mRMSD

        if hasattr(self, "_placement_ddG"):
            d["_placement_ddG"] = self._placement_ddG

        # if self.reactions:
        #     d['reactions'] = [r.dict for r in self.reactions]

        # if self.base:
        #     d['base'] = self.base.dict

        # if self.inspirations:
        #     d['inspirations'] = [i.longname for i in self.inspirations]

        return d

    @property
    def placement_ddG(self):
        return self._placement_ddG

    @property
    def placement_RMSD(self):
        return self._placement_mRMSD

    ### METHODS

    def calculate_fingerprint(self):

        # need to load the bound PDB
        sys = self.bound_system

        fingerprint = Fingerprint(self, self.mol, sys.protein_system)

        # fingerprint = None
        self._fingerprint = fingerprint

        return fingerprint

    def summary(self):

        mout.header(self)
        mout.var("smiles", self.smiles)
        mout.var("fingerprint", self.fingerprint)
        # mout.var('pdb_entry', self.pdb_entry)
        # mout.var('mol', self.mol)
        mout.var("name", self.name)
        mout.var("compound", str(self.compound))
        mout.var("pdb_path", str(self.pdb_path))
        mout.var("_site_index", self._site_index)
        mout.var("_chain", self._chain)
        mout.var("tags", self.tags)
        mout.var("_stereo_smiles", self._stereo_smiles)

    ### DUNDERS

    def __repr__(self):
        return f'Pose("{self.name}", compound={self.compound.name}, tags={self.tags})'
