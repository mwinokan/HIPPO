"""Fragalysis download service tests (SQLite tier, no network).

``requests.Session`` is replaced with a stub, so these check what
``DownloadService`` sends and how it reports a refused download, not Fragalysis
itself.
"""

import json

import pytest
import requests

pytestmark = pytest.mark.sqlite

ANON_MESSAGE = 'User "|anon|" is not a member of lb32627-71'


def make_response(status: int, body: dict | None = None) -> requests.Response:
    response = requests.Response()
    response.status_code = status
    response.url = "https://fragalysis.example/api/download_structures/"
    response._content = json.dumps(body).encode() if body is not None else b""
    return response


class StubSession:
    """Records the cookies sent and answers every POST with ``post_response``."""

    def __init__(self, post_response: requests.Response):
        self.headers = {}
        self.cookies = requests.cookies.RequestsCookieJar()
        self.post_response = post_response
        self.posted_cookies = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def get(self, *args, **kwargs):
        return make_response(200)

    def post(self, *args, **kwargs):
        self.posted_cookies = dict(self.cookies)
        return self.post_response


@pytest.fixture
def stub_session(animal, monkeypatch):
    """Patch the download module's ``requests.Session`` to a refusing stub."""
    from designdb.services import download

    session = StubSession(make_response(404, {"message": ANON_MESSAGE}))
    monkeypatch.setattr(download.requests, "Session", lambda: session)
    monkeypatch.delenv(download.AUTH_TOKEN_ENV, raising=False)
    return session


def download(tmp_path, **kwargs):
    from designdb.services.download import DownloadService

    return DownloadService.download_target(
        target_name="D68EV3C",
        target_access_string="lb32627-71",
        url="https://fragalysis.example",
        destination=tmp_path,
        **kwargs,
    )


def test_auth_token_falls_back_to_env(stub_session, tmp_path, monkeypatch):
    monkeypatch.setenv("FRAGALYSIS_AUTH_TOKEN", "env-session")
    with pytest.raises(requests.HTTPError):
        download(tmp_path)
    assert stub_session.posted_cookies == {"sessionid": "env-session"}


def test_explicit_auth_token_beats_env(stub_session, tmp_path, monkeypatch):
    monkeypatch.setenv("FRAGALYSIS_AUTH_TOKEN", "env-session")
    with pytest.raises(requests.HTTPError):
        download(tmp_path, auth_token="explicit-session")
    assert stub_session.posted_cookies == {"sessionid": "explicit-session"}


def test_no_auth_token_sends_no_session(stub_session, tmp_path):
    with pytest.raises(requests.HTTPError):
        download(tmp_path)
    assert stub_session.posted_cookies == {}


def test_refusal_reports_fragalysis_message_and_login_hint(stub_session, tmp_path):
    with pytest.raises(requests.HTTPError) as exc_info:
        download(tmp_path)
    message = str(exc_info.value)
    assert "404" in message
    assert ANON_MESSAGE in message
    assert "No Fragalysis login was sent" in message
    assert "FRAGALYSIS_AUTH_TOKEN" in message
    assert exc_info.value.response.status_code == 404


def test_refusal_with_login_hints_at_access(stub_session, tmp_path):
    with pytest.raises(requests.HTTPError) as exc_info:
        download(tmp_path, auth_token="stale-session")
    message = str(exc_info.value)
    assert ANON_MESSAGE in message
    assert "session is still valid" in message
    assert "No Fragalysis login was sent" not in message


def test_server_error_without_json_body(stub_session, tmp_path):
    stub_session.post_response = make_response(500)
    with pytest.raises(requests.HTTPError) as exc_info:
        download(tmp_path)
    message = str(exc_info.value)
    assert "500" in message
    assert "Fragalysis says" not in message
    assert "login" not in message
