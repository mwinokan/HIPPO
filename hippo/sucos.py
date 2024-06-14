
# implementation of SuCOS from https://github.com/susanhleung/SuCOS https://doi.org/10.26434/chemrxiv.8100203.v1

import os
from numpy import clip

from rdkit import RDConfig
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit.Chem.rdmolops import CombineMols

import logging
logger = logging.getLogger('HIPPO')

# feature setup

feature_factory = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))

feature_map_params = {}
for k in feature_factory.GetFeatureFamilies():
    feature_params = FeatMaps.FeatMapParams()
    feature_map_params[k] = feature_params

keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'ZnBinder',
        'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')

def feature_map_score(inspiration, derivative, score_mode=FeatMaps.FeatMapScoreMode.All, debug=False, draw=False):

	inspiration_features = [f for f in feature_factory.GetFeaturesForMol(inspiration) if f.GetFamily() in keep]
	derivative_features = [f for f in feature_factory.GetFeaturesForMol(derivative) if f.GetFamily() in keep]

	if draw:
		from molparse.rdkit import draw_mol
		draw_mol(inspiration, feats=inspiration_features)
		draw_mol(derivative, feats=derivative_features)
	
	inspiration_feature_map = FeatMaps.FeatMap(feats=inspiration_features, weights=[1] * len(inspiration_features), params=feature_map_params)
	derivative_feature_map = FeatMaps.FeatMap(feats=derivative_features, weights=[1] * len(derivative_features), params=feature_map_params)
	
	inspiration_feature_map.scoreMode = score_mode

	feature_score = inspiration_feature_map.ScoreFeats(derivative_features)

	score_vector = [0]*inspiration_feature_map.GetNumFeatures()
	
	inspiration_feature_map.ScoreFeats(derivative_features, mapScoreVect=score_vector)

	if debug:
		for feature, score in zip(inspiration_features, score_vector):
			logger.var(feature.GetFamily(), score)

	feature_score /= min(inspiration_feature_map.GetNumFeatures(), len(derivative_features))
	
	return feature_score

def multi_feature_map_score(inspirations, derivative, draw=False, debug=False, score_mode=FeatMaps.FeatMapScoreMode.All):

	inspiration_feature_lists = [[f for f in feature_factory.GetFeaturesForMol(x) if f.GetFamily() in keep] for x in inspirations]
	derivative_features = [f for f in feature_factory.GetFeaturesForMol(derivative) if f.GetFamily() in keep]

	combined_inspiration_features = sum(inspiration_feature_lists, [])

	if draw:
		from molparse.rdkit import draw_mol

		for inspiration, features in zip(inspirations, inspiration_feature_lists):
			draw_mol(inspiration, feats=features)
			
		draw_mol(derivative, feats=combined_inspiration_features)
		
		draw_mol(derivative, feats=derivative_features)

	derivative_feature_map = FeatMaps.FeatMap(feats=derivative_features, weights=[1] * len(derivative_features), params=feature_map_params)

	# loop over inspirations

	score_dict = {}

	for i,inspiration in enumerate(inspirations):

		score_dict[i] = {}

		if debug:
			logger.title(inspiration)

		inspiration_features = [f for f in feature_factory.GetFeaturesForMol(inspiration) if f.GetFamily() in keep]
		inspiration_feature_map = FeatMaps.FeatMap(feats=inspiration_features, weights=[1] * len(inspiration_features), params=feature_map_params)

		score_vector = [0]*derivative_feature_map.GetNumFeatures()

		feature_score = derivative_feature_map.ScoreFeats(inspiration_features, mapScoreVect=score_vector)

		for j, (feature, score) in enumerate(zip(derivative_features, score_vector)):
			if debug:
				logger.var(feature.GetFamily(), score)

			score_dict[i][j] = score

	i_list = list(score_dict.keys())

	combined_score_vector = []

	for j, feature in enumerate(derivative_features):

		scores = [score_dict[i][j] for i in i_list]

		best_score = max(scores)

		combined_score_vector.append(best_score)

		# print(feature.GetFamily(), scores, best_score)

	feature_score = sum(combined_score_vector)/len(derivative_features)
	
	return feature_score

def SuCOS_score(inspiration, derivative, return_all=False, print_scores=False, **kwargs):

	if isinstance(inspiration, list):
		multi=True
	else:
		multi=False

	if multi:
		feature_score = multi_feature_map_score(inspiration, derivative, **kwargs)
	else:
		feature_score = feature_map_score(inspiration, derivative, **kwargs)

	feature_score = clip(feature_score, 0, 1)
	
	if multi:

		mol = inspiration.pop()

		while inspiration:
			mol = CombineMols(mol, inspiration.pop())

		protrude_dist = rdShapeHelpers.ShapeProtrudeDist(mol, derivative, allowReordering=False)

	else:
		protrude_dist = rdShapeHelpers.ShapeProtrudeDist(inspiration, derivative, allowReordering=False)

	protrude_dist = clip(protrude_dist, 0, 1)
	volume_score = 1 - protrude_dist

	SuCOS_score = (feature_score + volume_score)*0.5

	if print_scores:
		logger.var("feature_score", feature_score)
		logger.var("volume_score", volume_score)
		logger.var("average_score", SuCOS_score)

	if return_all:
		return SuCOS_score, feature_score, volume_score
	else:
		return SuCOS_score
