"""Functions for interfacing with the apsw library"""

import apsw


def executemany(path: "Path", sql: str, payload: list[tuple]):
    """Bulk execution with apsw"""
    connection = apsw.Connection(str(path.resolve()))
    result = list(connection.executemany(sql, payload))
    connection.close()
    return result
