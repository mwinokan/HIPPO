import logging
from collections import defaultdict

from designdb.models import InspirationModel, SubsiteModel, SubsiteTagModel

logger = logging.getLogger(__name__)


class SubsiteService:
    @classmethod
    def set_derivative_subsites(cls) -> None:
        """Propagate subsite assignments from inspiration poses to their derivatives."""
        subsite_tags = SubsiteTagModel.objects.filter(
            pose__in=InspirationModel.objects.values('original_pose')
        ).select_related('subsite')

        inspirations = InspirationModel.objects.filter(
            original_pose__in=subsite_tags.values('pose')
        ).select_related('derivative_pose')

        original_to_derivatives = defaultdict(list)
        for insp in inspirations:
            original_to_derivatives[insp.original_pose_id].append(insp.derivative_pose)

        new_tags = [
            SubsiteTagModel(pose=derivative, subsite=st.subsite)
            for st in subsite_tags
            for derivative in original_to_derivatives[st.pose_id]
        ]

        SubsiteTagModel.objects.bulk_create(new_tags, ignore_conflicts=True)

    @classmethod
    def set_subsites_from_metadata_field(
        cls,
        poses,
        field: str = 'CanonSites alias',
    ) -> None:
        """Create and assign subsite entries from a pose metadata field."""
        assignments = []

        for pose in poses.select_related('target'):
            metadata = pose.pose_metadata or {}
            name = metadata.get(field)
            if not name:
                logger.warning('Field "%s" not in metadata for pose pk=%s', field, pose.pk)
                continue
            subsite, _ = SubsiteModel.objects.get_or_create(
                target=pose.target,
                subsite_name=name,
            )
            assignments.append(SubsiteTagModel(pose=pose, subsite=subsite))

        SubsiteTagModel.objects.bulk_create(assignments, ignore_conflicts=True)
