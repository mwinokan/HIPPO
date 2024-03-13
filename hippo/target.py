
import mcol

class Target:

	def __init__(self, db, id, name):

		self._db = db
		self._id = id
		self._name = name

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
	def feature_ids(self):
		feature_ids = self.db.select_where(
			query='feature_id', 
			table='feature', 
			key='target', 
			value=self.id, 
			none=False, 
			multiple=True,
			sort='feature_chain_name, feature_residue_number'
		)

		if feature_ids:
			feature_ids = [v for v, in feature_ids]

		return feature_ids

	@property
	def features(self):
		if feature_ids := self.feature_ids:
			return [self.db.get_feature(id=i) for i in feature_ids]
		return None

	def __str__(self):
		return f'T{self.id}'

	def __repr__(self):
		return f'{mcol.bold}{mcol.underline}{self} "{self.name}"{mcol.unbold}{mcol.ununderline}'

	def calculate_features(self, protein):

		features = protein.get_protein_features()

		for f in features:
			self.db.insert_feature(
				family=f.family,
				target=self.id,
				atom_names=[a.name for a in f.atoms],
				residue_name=f.res_name,
				residue_number=f.res_number,
				chain_name=f.res_chain,
			)

		return self.features
