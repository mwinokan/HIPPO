from ..models import CampaignModel
from .campaign_set import CampaignSet


class Campaign(CampaignModel):

    _objects = CampaignSet.as_manager()
    _parent_module = "projects"
