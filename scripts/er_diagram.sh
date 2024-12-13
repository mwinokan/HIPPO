#!/bin/bash

EXCLUDE="AbstractModel,TargetModel,CompoundModel,SupplierModel,PoseModel,QuoteModel,TagModel,ReactionModel,ProductModel,ReactantModel,FeatureModel,InteractionModel,SubsiteModel,ObservationModel,SolventModel,SupplierModel,PlacementModel,StructureModel"

echo "writing docs/images/hippo_models.dot"
python manage.py graph_models hippo --no-inheritance --exclude-models $EXCLUDE > docs/images/hippo_models.dot
echo "visualise on dreampuf.github.io/GraphvizOnline"

