from pony import orm

db = orm.Database()


class Compound(db.Entity):

    name = orm.Required(str)
    smiles = orm.Required(str)
