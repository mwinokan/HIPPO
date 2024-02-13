
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
		self._compound = compound
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
	def compound(self):
		return self._compound

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
	
	### METHODS

	def get_tags(self):
		tags = self.db.select_where(query='tag_name', table='tag', key='pose', value=self.id, multiple=True)
		return {t[0] for t in tags}

	### DUNDERS

	def __repr__(self):
		return f'Pose #{self.id} "{self.longname}"'
		