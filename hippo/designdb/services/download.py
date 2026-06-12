"""Service layer for downloading target data from Fragalysis.

Wraps the Fragalysis ``/api/download_structures/`` endpoint, which builds a
(g)zipped archive on demand and returns it in two steps:

1. ``POST`` the download specification -> the response JSON contains a
   ``file_url`` once the archive is ready (this is *not* async on the server, so
   the POST can block while the archive is assembled).
2. ``GET`` that ``file_url`` -> stream the archive to disk.

The endpoint's serializer exposes many boolean flags. For hippo we only need the
apo-desolvated protein structures, so :meth:`.DownloadService.download_target`
defaults ``apo_desolv_file=True`` and every other flag to ``False``. Individual
flags can still be overridden per call via keyword arguments.

Adapted and hardened from the ``downloader.py`` prototype.
"""

import os
import tarfile
import time
import zipfile
from pathlib import Path
from urllib.parse import urljoin

import mrich
import requests
from designdb.utils_frag import STACK_URLS
from requests.exceptions import JSONDecodeError

LOGIN_URL = '/accounts/login/'
DOWNLOAD_URL = '/api/download_structures/'
LANDING_PAGE_URL = '/viewer/react/landing/'

# Keep reasonably current; some endpoints reject obviously-bot user agents.
USER_AGENT = (
    'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) '
    'Chrome/120.0.0.0 Safari/537.36'
)

# Every BooleanField on DownloadStructuresSerializer except ``use_zip`` (which is
# handled as a dedicated parameter). This tuple is the single source of truth for
# "apo_desolv default, everything else False" and should track the serializer.
BOOLEAN_FLAGS = (
    'all_aligned_structures',
    'apo_file',
    'bound_file',
    'apo_solv_file',
    'apo_desolv_file',
    'ligand_pdb',
    'ligand_sdf',
    'ligand_smiles',
    'sdf_info',
    'smiles_info',
    'pdb_info',
    'cif_info',
    'mtz_info',
    'diff_file',
    'event_file',
    'sigmaa_file',
    'map_info',
    'single_sdf_file',
    'metadata_info',
    'trans_matrix_info',
    'compound_sets',
    'soakdb_files',
    'yaml_files',
    'extra_files',
    'pymol_scripts',
    'readme',
    'static_link',
)

# The flag enabled by default by :meth:`.DownloadService.download_target`.
DEFAULT_FLAG = 'apo_desolv_file'

# (connect timeout, read timeout) in seconds.
DEFAULT_TIMEOUT = (30, 1800)

# Short timeout for the best-effort CSRF-priming GET. It hits a normal page just
# to obtain the csrftoken cookie, so it must fail fast (e.g. on a backend-only
# deployment that doesn't serve the frontend landing page) rather than block on
# the long download read timeout.
CSRF_TIMEOUT = (10, 10)

# Task-status polling. The POST only *triggers* archive creation and returns a
# task status URL; we poll it until the task reaches a terminal state and a
# file_url becomes available.
FAILURE_STATUSES = ('FAILED', 'CANCELED', 'FATAL')
DEFAULT_POLL_INTERVAL = 2
DEFAULT_POLL_TIMEOUT = 1800


class DownloadService:
    """Download target data archives from Fragalysis."""

    @staticmethod
    def _build_payload(
        *,
        target_name: str,
        target_access_string: str = '',
        proteins: str = '',
        use_zip: bool = False,
        **flag_overrides: bool,
    ) -> dict:
        """Build the request payload.

        Every boolean flag defaults to ``False`` except :data:`.DEFAULT_FLAG`
        (``apo_desolv_file``). Pass any flag in :data:`.BOOLEAN_FLAGS` as a keyword
        argument to override it.

        :param target_name: Fragalysis target name
        :param target_access_string: target access string (proposal)
        :param proteins: comma-separated observation shortcodes (empty = all)
        :param use_zip: request a ``.zip`` instead of a ``.tar.gz`` archive
        :param flag_overrides: per-call overrides for any flag in
            :data:`.BOOLEAN_FLAGS`
        """

        unknown = set(flag_overrides) - set(BOOLEAN_FLAGS)
        if unknown:
            raise TypeError(
                f'Unknown download flag(s): {sorted(unknown)}. '
                f'Valid flags: {", ".join(BOOLEAN_FLAGS)}'
            )

        payload = {flag: False for flag in BOOLEAN_FLAGS}
        payload[DEFAULT_FLAG] = True
        payload.update(flag_overrides)

        payload.update(
            target_name=target_name,
            target_access_string=target_access_string or '',
            proteins=proteins or '',
            # NB: must be '' (not False) on the initial request, see prototype
            file_url='',
            use_zip=bool(use_zip),
        )

        return payload

    @staticmethod
    def download_target(
        *,
        target_name: str,
        target_access_string: str = '',
        proteins: str = '',
        stack: str = 'production',
        url: str | None = None,
        auth_token: str | None = None,
        destination: 'str | Path | None' = None,
        extract: bool = True,
        use_zip: bool = False,
        timeout: 'tuple[int, int] | int | None' = DEFAULT_TIMEOUT,
        poll_interval: float = DEFAULT_POLL_INTERVAL,
        poll_timeout: float = DEFAULT_POLL_TIMEOUT,
        **flag_overrides: bool,
    ) -> Path:
        """Download (and optionally extract) a target archive from Fragalysis.

        :param target_name: Fragalysis target name
        :param target_access_string: target access string (proposal)
        :param proteins: comma-separated observation shortcodes (empty = all)
        :param stack: key into :data:`.STACK_URLS` (``'production'``/``'staging'``)
        :param url: explicit base URL, overrides ``stack`` if given
        :param auth_token: Fragalysis ``sessionid`` cookie; falls back to the
            ``FRAGALYSIS_AUTH_TOKEN`` environment variable
        :param destination: directory to download into (default: current dir)
        :param extract: extract the archive and return the extracted directory
        :param use_zip: request a ``.zip`` instead of a ``.tar.gz`` archive
        :param timeout: requests timeout (connect, read) in seconds
        :param poll_interval: seconds between task-status polls
        :param poll_timeout: max seconds to wait for the archive task to finish
        :param flag_overrides: per-call overrides for any flag in
            :data:`.BOOLEAN_FLAGS`
        :returns: path to the extracted directory (``extract=True``) or the
            downloaded archive (``extract=False``)
        """

        base_url = url or STACK_URLS.get(stack)
        if not base_url:
            raise ValueError(
                f'Unknown stack {stack!r}; choose from {sorted(STACK_URLS)} '
                'or pass an explicit url'
            )


        destination = Path(destination) if destination else Path.cwd()
        destination.mkdir(parents=True, exist_ok=True)

        payload = DownloadService._build_payload(
            target_name=target_name,
            target_access_string=target_access_string,
            proteins=proteins,
            use_zip=use_zip,
            **flag_overrides,
        )

        download_api_url = urljoin(base_url, DOWNLOAD_URL)
        landing_page_url = urljoin(base_url, LANDING_PAGE_URL)

        mrich.var('target_name', target_name)
        mrich.var('target_access_string', target_access_string)
        mrich.var('stack url', base_url)

        with requests.Session() as session:
            session.headers.update(
                {
                    'User-Agent': USER_AGENT,
                    'Referer': landing_page_url,
                    'Referrer-policy': 'same-origin',
                }
            )

            # Best-effort: prime the csrftoken cookie by hitting a normal page.
            # Use a short timeout and tolerate failure so a backend-only stack
            # that doesn't serve the landing page fails fast instead of hanging
            # on the long download read timeout. If no token is obtained we still
            # proceed; a CSRF-enforcing server would then return a clear error.
            try:
                session.get(landing_page_url, timeout=CSRF_TIMEOUT)
            except requests.RequestException as exc:
                mrich.warning(
                    f'Could not reach {landing_page_url} to obtain a CSRF token '
                    f'({exc}); proceeding without it.'
                )
            csrftoken = session.cookies.get('csrftoken', None)
            if csrftoken:
                session.headers.update({'X-CSRFToken': csrftoken})

            if auth_token:
                session.cookies.update({'sessionid': auth_token})

            # Step 1: trigger archive creation. This returns a task status URL
            # (the archive is built asynchronously on the server).
            mrich.print('Requesting download from Fragalysis...')
            start_response = session.post(
                download_api_url, data=payload, timeout=timeout
            )
            start_response.raise_for_status()
            start_json = start_response.json()

            # Step 2: resolve the file_url. Newer servers return a task status
            # URL to poll; older/cached responses may return file_url directly.
            file_url = start_json.get('file_url')
            if not file_url:
                task_status_url = start_json.get('task_status_url')
                if not task_status_url:
                    raise RuntimeError(
                        'Fragalysis returned neither file_url nor task_status_url: '
                        f'{start_json}'
                    )
                task_status_url = urljoin(base_url, task_status_url)
                file_url = DownloadService._poll_task(
                    session,
                    task_status_url,
                    timeout=timeout,
                    poll_interval=poll_interval,
                    poll_timeout=poll_timeout,
                )

            # Step 3: stream the archive to disk
            archive_path = destination / Path(file_url).name
            mrich.writing(archive_path)
            with session.get(
                download_api_url,
                params={'file_url': file_url},
                stream=True,
                timeout=timeout,
            ) as r:
                r.raise_for_status()
                with open(archive_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)

        mrich.success('Downloaded', archive_path)

        if not extract:
            return archive_path

        return DownloadService._extract(archive_path, destination)

    @staticmethod
    def _poll_task(
        session: requests.Session,
        task_status_url: str,
        *,
        timeout: 'tuple[int, int] | int | None' = DEFAULT_TIMEOUT,
        poll_interval: float = DEFAULT_POLL_INTERVAL,
        poll_timeout: float = DEFAULT_POLL_TIMEOUT,
    ) -> str:
        """Poll a Fragalysis task-status URL until the archive is ready.

        :param session: the authenticated :class:`requests.Session`
        :param task_status_url: absolute task-status URL returned by the POST
        :param timeout: per-request timeout
        :param poll_interval: seconds between polls
        :param poll_timeout: max seconds to wait before giving up
        :returns: the ``file_url`` of the assembled archive
        """

        mrich.print('Waiting for Fragalysis to build the archive...')

        waited = 0.0
        seen_messages = 0

        while True:
            resp = session.get(task_status_url, timeout=timeout)
            try:
                data = resp.json()
            except JSONDecodeError:
                # task is too early in its lifecycle to have a JSON body yet
                data = {}

            if error := data.get('error'):
                raise RuntimeError(f'Fragalysis download task failed: {error}')

            status = data.get('status')
            if status in FAILURE_STATUSES:
                raise RuntimeError(f'Fragalysis download task {status}: {data}')

            finished = (
                status == 'SUCCESS'
                or data.get('finished') is True
                or data.get('ready') is True
            )

            if finished:
                # On success the download endpoint delivers the archive path in
                # `messages` (a string), not a dedicated `file_url` key.
                file_url = DownloadService._extract_file_url(data)
                if not file_url:
                    raise RuntimeError(
                        f'Download task finished but returned no file_url: {data}'
                    )
                return file_url

            # not finished yet: surface any new progress messages
            messages = data.get('messages')
            if messages is not None:
                messages = messages if isinstance(messages, list) else [messages]
                for m in messages[seen_messages:]:
                    mrich.print(m)
                seen_messages = len(messages)

            if waited >= poll_timeout:
                raise TimeoutError(
                    f'Download task did not finish within {poll_timeout}s '
                    f'({task_status_url})'
                )

            time.sleep(poll_interval)
            waited += poll_interval

    @staticmethod
    def _extract_file_url(data: dict) -> 'str | None':
        """Extract the archive path from a completed task-status response.

        The download endpoint returns the path in ``file_url`` on some servers and
        in ``messages`` (a string, occasionally a list) on others.
        """
        file_url = data.get('file_url')
        if file_url:
            return file_url

        messages = data.get('messages')
        if isinstance(messages, str):
            return messages
        if isinstance(messages, list) and messages:
            return messages[-1]
        return None

    @staticmethod
    def _extract(archive_path: Path, destination: Path) -> Path:
        """Extract a ``.zip`` or ``.tar(.gz)`` archive into a sibling directory.

        :param archive_path: path to the downloaded archive
        :param destination: directory to extract into
        :returns: path to the extracted directory
        """

        extract_dir = destination / _archive_stem(archive_path.name)
        extract_dir.mkdir(parents=True, exist_ok=True)

        mrich.print('Extracting', archive_path.name, '->', extract_dir)

        if zipfile.is_zipfile(archive_path):
            with zipfile.ZipFile(archive_path) as zf:
                zf.extractall(extract_dir)
        elif tarfile.is_tarfile(archive_path):
            with tarfile.open(archive_path) as tf:
                tf.extractall(extract_dir)
        else:
            raise ValueError(f'Unrecognised archive format: {archive_path}')

        mrich.success('Extracted to', extract_dir)
        return extract_dir


def _archive_stem(name: str) -> str:
    """Strip a ``.zip`` / ``.tar`` / ``.tar.gz`` / ``.tgz`` suffix from a filename."""
    for suffix in ('.tar.gz', '.tar.bz2', '.tgz', '.tar', '.zip', '.gz'):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem
