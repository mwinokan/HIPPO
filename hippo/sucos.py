
# implementation of SuCOS from https://github.com/susanhleung/SuCOS https://doi.org/10.26434/chemrxiv.8100203.v1

import os
from numpy import clip

from rdkit import RDConfig
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps

import logging
logger = logging.getLogger('HIPPO')

# feature setup

fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))

fmParams = {}
for k in fdef.GetFeatureFamilies():
    fparams = FeatMaps.FeatMapParams()
    fmParams[k] = fparams

keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'ZnBinder',
        'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')

def get_FeatureMapScore(inspiration, derivative, score_mode=FeatMaps.FeatMapScoreMode.All):
	featLists = []
	for m in [inspiration, derivative]:
		rawFeats = fdef.GetFeaturesForMol(m)
		# filter that list down to only include the ones we're intereted in
		featLists.append([f for f in rawFeats if f.GetFamily() in keep])
	fms = [FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=fmParams) for x in featLists]
	fms[0].scoreMode = score_mode
	feature_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))
	return feature_score

def SuCOS_score(inspiration, derivative, return_all=False, print_scores=False):

	feature_score = get_FeatureMapScore(inspiration, derivative)
	feature_score = clip(feature_score, 0, 1)
	
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
