from pathlib import Path


def dict_formatter(payload) -> dict:
    if not payload:
        return {}

    if isinstance(payload, str):
        import json

        payload = json.loads(payload)

    elif not isinstance(payload, dict):
        payload = dict(payload)

    # payload could still be "None" if "null" is passed to json.loads
    payload = payload or {}

    assert isinstance(payload, dict)

    return payload


def str_formatter(payload) -> str:
    if not payload:
        return ""

    return str(payload)


def list_formatter(payload, separator: str | None = ",") -> list:
    if not payload:
        return []

    if isinstance(payload, str):
        import json

        try:
            payload = json.loads(payload)
        except json.JSONDecodeError:
            if separator:
                payload = payload.split(separator)
            else:
                raise

    elif not isinstance(payload, list):
        payload = list(payload)

    assert isinstance(payload, list)

    return payload


def path_formatter(payload, resolve: bool = True, exists: bool = False):

    if not payload:
        return None

    path = Path(payload)

    if resolve:
        path = path.resolve()

    if exists and not path.exists():
        raise FileNotFoundError(path)

    return path
