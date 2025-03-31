from ..abstract import AbstractQuerySet
import mrich


class PoseSet(AbstractQuerySet):

    ### FACTORIES

    ### PROPERTIES

    ### METHODS

    def score_set(self, method="mocassin", *args, **kwargs):

        match method:
            case "mocassin":

                from mocassin.mocassin import calculate_mocassin

                query_pose = self.last()
                ref_poses = self.exclude(id=self.last().id)

                combined, shape, feature = calculate_mocassin(
                    [p.mol for p in ref_poses], query_pose.mol, *args, **kwargs
                )

                mrich.debug(combined, shape, feature)

                return dict(combined=combined, shape=shape, feature=feature)

            case _:
                raise NotImplementedError("Scoring method not supported '{method}'")

    ### DUNDERS
