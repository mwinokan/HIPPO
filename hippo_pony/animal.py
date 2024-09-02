# from .db import db

from .models import db


class HIPPO:

    def __init__(self, db_path):

        self._db = db

        db.bind(provider="sqlite", filename=db_path, create_db=True)

        db.generate_mapping(create_tables=True)

    @property
    def db(self):
        return self._db
