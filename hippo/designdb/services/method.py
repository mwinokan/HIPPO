import logging

from designdb.models import EnumerationMethodModel, PoseMethodModel, ScoringMethodModel

logger = logging.getLogger(__name__)


class MethodService:
    @classmethod
    def register_enumeration_method(
        cls, name: str, version: str, description: str = ''
    ):
        obj, created = EnumerationMethodModel.objects.get_or_create(
            enum_name=name,
            enum_version=version,
            defaults={'enum_description': description},
        )
        return obj, created

    @classmethod
    def register_pose_method(cls, name: str, version: str, description: str = ''):
        obj, created = PoseMethodModel.objects.get_or_create(
            pose_method_name=name,
            pose_method_version=version,
            defaults={'pose_method_description': description},
        )
        return obj, created

    @classmethod
    def register_scoring_method(cls, name: str, version: str, description: str = ''):
        obj, created = ScoringMethodModel.objects.get_or_create(
            method_name=name,
            method_version=version,
            defaults={'method_description': description},
        )
        return obj, created

    @classmethod
    def get_enumeration_methods(cls):
        return EnumerationMethodModel.objects.all()

    @classmethod
    def get_pose_methods(cls):
        return PoseMethodModel.objects.all()

    @classmethod
    def get_scoring_methods(cls):
        return ScoringMethodModel.objects.all()
