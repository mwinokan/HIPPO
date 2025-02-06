import apsw


def executemany(path, sql, payload):
    connection = apsw.Connection(str(path.resolve()))
    result = list(connection.executemany(sql, payload))
    connection.close()
    return result
