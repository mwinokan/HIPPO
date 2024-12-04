#!/bin/bash

EXCLUDE="TargetModel,CompoundModel,SupplierModel,PoseModel,QuoteModel,TagModel,ReactionModel,ProductModel,ReactantModel,FeatureModel,InteractionModel,SubsiteModel,ObservationModel,SolventModel,SupplierModel"

python manage.py graph_models hippo --no-inheritance --exclude-models $EXCLUDE > hippo_models.dot
# python manage.py graph_models hippo --disable-abstract-fields > hippo_models.dot
# python manage.py graph_models hippo > hippo_models.dot
echo "visualise on dreampuf.github.io/GraphvizOnline"

