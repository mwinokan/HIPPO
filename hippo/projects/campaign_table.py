from .campaign import Campaign
from ..abstract import AbstractTable


class CampaignTable(AbstractTable):

    _name = "all campaigns"
    _model = Campaign
