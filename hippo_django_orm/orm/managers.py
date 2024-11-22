import mrich


class ManagerRouter:
    """Route attribute calls to _<field> via the custom manager stored in "_objects" on the model,
    except for those in REROUTED_ATTRIBUTES which go directly via the default RelatedManager
    """

    REROUTED_ATTRIBUTES = {"add", "remove", "clear", "set"}

    def __init__(self, model, field):
        self.model = model
        self.default_manager = getattr(self.model, f"_{field}")
        self.custom_manager = getattr(self.model, f"_{field}")(manager="_objects").all()

    def __getattr__(self, attr):
        if attr in self.REROUTED_ATTRIBUTES:
            descriptor = self.default_manager
        else:
            descriptor = self.custom_manager
        return getattr(descriptor, attr)

    def __getitem__(self, key):
        return self.custom_manager[key]

    def __str__(self):
        return self.custom_manager.__str__()

    def __repr__(self):
        return self.custom_manager.__repr__()

    def __rich__(self):
        return self.custom_manager.__rich__()
