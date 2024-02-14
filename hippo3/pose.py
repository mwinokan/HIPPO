
import mcol

class Pose:

	def __init__(self, 
		db,
		id: int,
		name: str,
		longname: str,
		smiles: str,
		compound: int,
		target: str,
		path: str,
		mol: str,
		fingerprint: str,
	):

		self._db = db
		self._id = id
		self._name = name
		self._longname = longname
		self._smiles = smiles
		self._compound_id = compound
		self._target = target
		self._path = path
		self._mol = mol
		self._fingerprint = fingerprint
		
	### FACTORIES

	### PROPERTIES

	@property
	def db(self):
		return self._db

	@property
	def id(self):
		return self._id

	@property
	def name(self):
		return self._name

	@property
	def longname(self):
		return self._longname

	@property
	def smiles(self):
		return self._smiles

	@property
	def target(self):
		return self._target

	@property
	def compound_id(self):
		return self._compound_id

	@property
	def compound(self):
		return self.get_compound()

	@property
	def path(self):
		return self._path

	@property
	def mol(self):
		return self._mol

	@property
	def fingerprint(self):
		return self._fingerprint

	@property
	def tags(self):
		return self.get_tags()

	@property
	def inspirations(self):
		return self.get_inspirations()


	### METHODS

	def get_compound(self):
		return self.db.get_compound(id=self._compound_id)

	def get_tags(self):
		tags = self.db.select_where(query='tag_name', table='tag', key='pose', value=self.id, multiple=True)
		return {t[0] for t in tags}
	
	def get_inspirations(self):
		inspirations = self.db.select_where(query='inspiration_original', table='inspiration', key='derivative', value=self.id, multiple=True, none='quiet')
		
		if inspirations:
			inspirations = [self.db.get_pose(id=id[0]) for id in inspirations]

		return inspirations

	### DUNDERS

	def __str__(self):
		return f'P{self.id}'

	def __repr__(self):
		return f'{mcol.bold}{mcol.underline}{self.compound}->{self} "{self.longname}"{mcol.unbold}{mcol.ununderline}'
		