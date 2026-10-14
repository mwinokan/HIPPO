import logging

from designdb.models import EnumerationMethodModel, PoseMethodModel, ScoringMethodModel
from designdb.utils import admin_write

logger = logging.getLogger(__name__)


class MethodService:
    @classmethod
    def register_enumeration_method(
        cls, name: str, version: str, description: str = ''
    ):
        with admin_write(f'Registering enumeration method {name!r} v{version}'):
            obj, created = EnumerationMethodModel.objects.get_or_create(
                enum_name=name,
                enum_version=version,
                defaults={'enum_description': description},
            )
        return obj, created

    @classmethod
    def register_pose_method(cls, name: str, version: str, description: str = ''):
        with admin_write(f'Registering pose method {name!r} v{version}'):
            obj, created = PoseMethodModel.objects.get_or_create(
                pose_method_name=name,
                pose_method_version=version,
                defaults={'pose_method_description': description},
            )
        return obj, created

    @classmethod
    def register_scoring_method(cls, name: str, version: str, description: str = ''):
        with admin_write(f'Registering scoring method {name!r} v{version}'):
            obj, created = ScoringMethodModel.objects.get_or_create(
                method_name=name,
                method_version=version,
                defaults={'method_description': description},
            )
        return obj, created

    @classmethod
    def resolve_enumeration_method(
        cls, method: tuple[str, str] | None
    ) -> EnumerationMethodModel | None:
        """Resolve a ``(name, version)`` pair to a registered enumeration method.

        Lookup counterpart of :meth:`register_enumeration_method`, for callers that
        take a method as user input and need the row it refers to.

        :param method: ``(name, version)``, or ``None`` for no method
        :returns: the method, or ``None`` if ``method`` was ``None``
        :raises ValueError: if the method is not registered
        """
        if method is None:
            return None
        name, version = method
        try:
            return EnumerationMethodModel.objects.get(
                enum_name=name, enum_version=version
            )
        except EnumerationMethodModel.DoesNotExist:
            raise ValueError(
                f"Enumeration method '{name}' v{version} not found. "
                'Call register_enumeration_method() first.'
            ) from None

    @classmethod
    def resolve_pose_method(
        cls, method: tuple[str, str] | None
    ) -> PoseMethodModel | None:
        """Resolve a ``(name, version)`` pair to a registered pose method.

        :param method: ``(name, version)``, or ``None`` for no method
        :returns: the method, or ``None`` if ``method`` was ``None``
        :raises ValueError: if the method is not registered
        """
        if method is None:
            return None
        name, version = method
        try:
            return PoseMethodModel.objects.get(
                pose_method_name=name, pose_method_version=version
            )
        except PoseMethodModel.DoesNotExist:
            raise ValueError(
                f"Pose method '{name}' v{version} not found. "
                'Call register_pose_method() first.'
            ) from None

    @classmethod
    def resolve_scoring_method(cls, name: str, version: str) -> ScoringMethodModel:
        """Resolve a name and version to a registered scoring method.

        :param name: scoring method name
        :param version: scoring method version
        :returns: the scoring method
        :raises ValueError: if the method is not registered
        """
        try:
            return ScoringMethodModel.objects.get(
                method_name=name, method_version=version
            )
        except ScoringMethodModel.DoesNotExist:
            raise ValueError(
                f"Scoring method '{name}' v{version} not found. "
                'Call register_scoring_method() first.'
            ) from None

    @classmethod
    def get_enumeration_methods(cls):
        return EnumerationMethodModel.objects.all()

    @classmethod
    def get_pose_methods(cls):
        return PoseMethodModel.objects.all()

    @classmethod
    def get_scoring_methods(cls):
        return ScoringMethodModel.objects.all()
