#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Standalone DB tunnel diagnostics for APERO RI.

This script mirrors the APERO RI database tunnel flow in a self-contained
way so you can debug SSH tunnel startup and MySQL connectivity without the
web app.
"""

import argparse
import importlib
import shlex
import socket
import subprocess
import sys
import time
import traceback
from pathlib import Path

# =============================================================================
# Define variables
# =============================================================================
DEFAULT_SSH_HOST = 'fir-login'
DEFAULT_REMOTE_HOST = 'cedar-mysql-vm.int.cedar.computecanada.ca'
DEFAULT_DB_USER = 'lmalo'
DEFAULT_DB_NAME = 'lmalo_spirou'
DEFAULT_REMOTE_PORT = 3306
DEFAULT_LOCAL_PORT = 3307
DEFAULT_PASSWORD_FILE = Path(__file__).with_name('pass.txt')


# =============================================================================
# Define functions
# =============================================================================
def _parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            'Standalone DB tunnel diagnostic for APERO RI. '
            'Use --auth, --verify, and --test to run individual steps.'
        )
    )
    parser.add_argument(
        '--auth',
        action='store_true',
        help=(
            'Run the interactive SSH tunnel command in the foreground '
            'so you can complete password/MFA prompts manually.'
        ),
    )
    parser.add_argument(
        '--verify',
        action='store_true',
        help=(
            'Run the background SSH tunnel command and verify that the '
            'local forwarded port is open.'
        ),
    )
    parser.add_argument(
        '--test',
        action='store_true',
        help=(
            'Run the APERO-style DB test using SQLAlchemy and a direct '
            'PyMySQL fallback.'
        ),
    )
    parser.add_argument(
        '--ssh-host',
        default=DEFAULT_SSH_HOST,
        help='SSH config host alias used for the tunnel.',
    )
    parser.add_argument(
        '--remote-host',
        default=DEFAULT_REMOTE_HOST,
        help='Remote MySQL host behind the SSH tunnel.',
    )
    parser.add_argument(
        '--remote-port',
        type=int,
        default=DEFAULT_REMOTE_PORT,
        help='Remote MySQL port.',
    )
    parser.add_argument(
        '--local-port',
        type=int,
        default=DEFAULT_LOCAL_PORT,
        help='Local forwarded port.',
    )
    parser.add_argument(
        '--db-user',
        default=DEFAULT_DB_USER,
        help='Database username.',
    )
    parser.add_argument(
        '--db-name',
        default=DEFAULT_DB_NAME,
        help='Database name to test.',
    )
    parser.add_argument(
        '--password-file',
        default=str(DEFAULT_PASSWORD_FILE),
        help='Path to a one-line password file (default: pass.txt).',
    )
    parser.add_argument(
        '--password',
        default='',
        help='Password override; if omitted the script reads password-file.',
    )
    parser.add_argument(
        '--no-wait',
        action='store_true',
        help='Do not wait for the forwarded port after verify.',
    )
    parser.add_argument(
        '--simple-ssh',
        action='store_true',
        help=(
            'Use docs-style SSH tunnel commands (no ControlMaster socket '
            'management).'
        ),
    )
    return parser.parse_args()


def _read_password(password_file: Path, override: str) -> str:
    """Read the password from a file or use an explicit override."""
    if override:
        return str(override)
    if password_file.exists():
        raw = password_file.read_text(encoding='utf-8', errors='replace')
        for line in raw.splitlines():
            value = line.strip()
            if value:
                return value
    return ''


def _mask_password(password: str) -> str:
    """Return a masked password for logs."""
    if not password:
        return '<empty>'
    return '*' * max(4, min(12, len(password)))


def _format_ssh_command(
    ssh_host: str,
    remote_host: str,
    remote_port: int,
    local_port: int,
    *,
    background: bool,
    interactive: bool,
    simple_ssh: bool,
) -> list[str]:
    """Build the SSH command used for the tunnel."""
    if simple_ssh:
        return _format_simple_ssh_command(
            ssh_host,
            remote_host,
            remote_port,
            local_port,
            background=background,
            interactive=interactive,
        )

    cmd = ['ssh']
    if interactive:
        cmd.extend(['-vvv', '-t', '-t'])
    if background:
        cmd.append('-f')
    cmd.extend(
        [
            '-N',
            '-M',
            '-S',
            str(
                _control_socket_path(
                    ssh_host,
                    remote_host,
                    local_port,
                    remote_port,
                )
            ),
            '-o',
            'StrictHostKeyChecking=accept-new',
            '-o',
            'ControlPersist=yes',
            '-o',
            'ExitOnForwardFailure=yes',
            '-o',
            'ConnectTimeout=15',
            '-L',
            f'{local_port}:{remote_host}:{remote_port}',
        ]
    )
    if background:
        cmd.extend(['-o', 'BatchMode=yes'])
    cmd.append(ssh_host)
    return cmd


def _format_simple_ssh_command(
    ssh_host: str,
    remote_host: str,
    remote_port: int,
    local_port: int,
    *,
    background: bool,
    interactive: bool,
) -> list[str]:
    """Build a docs-style SSH command with minimal options."""
    cmd = ['ssh']
    if interactive:
        cmd.append('-vvv')
    if background:
        cmd.extend(['-f', '-N', '-o', 'ExitOnForwardFailure=yes'])
    cmd.extend(['-L', f'{local_port}:{remote_host}:{remote_port}', ssh_host])
    return cmd


def _control_socket_path(
    ssh_host: str,
    remote_host: str,
    local_port: int,
    remote_port: int,
) -> Path:
    """Return the same style of tunnel control socket path APERO uses."""
    data_root = Path.home() / '.ari' / 'secret' / 'db_tunnels'
    data_root.mkdir(parents=True, exist_ok=True)
    import hashlib

    signature = hashlib.sha1(
        f'{ssh_host}|{remote_host}|{local_port}|{remote_port}'.encode('utf-8')
    ).hexdigest()[:16]
    return data_root / f'{signature}.sock'


def _show_command(cmd: list[str]) -> None:
    """Print a shell-friendly command line."""
    print(shlex.join(cmd))


def _run_command(cmd: list[str], label: str) -> subprocess.CompletedProcess:
    """Run a command and dump stdout/stderr."""
    print(f'\n[{label}] Running:')
    _show_command(cmd)
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    print(f'[{label}] Exit code: {result.returncode}')
    if result.stdout:
        print(f'[{label}] stdout:')
        print(result.stdout.rstrip())
    if result.stderr:
        print(f'[{label}] stderr:')
        print(result.stderr.rstrip())
    return result


def _port_open(host: str, port: int, timeout: float = 0.5) -> bool:
    """Return True when a local port accepts connections."""
    try:
        with socket.create_connection((host, port), timeout=timeout):
            return True
    except Exception:
        return False


def _wait_for_port(port: int, timeout_s: int = 20) -> bool:
    """Wait for localhost:port to become reachable."""
    print(f'[VERIFY] Waiting for 127.0.0.1:{port} ...')
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        if _port_open('127.0.0.1', port):
            print(f'[VERIFY] Port {port} is open.')
            return True
        time.sleep(0.5)
    print(f'[VERIFY] Timed out waiting for port {port}.')
    return False


def _show_listeners(port: int) -> None:
    """Print listener lines for the forwarded port."""
    try:
        result = subprocess.run(
            ['ss', '-lptn'], capture_output=True, text=True, check=False
        )
    except FileNotFoundError:
        print('[VERIFY] ss is not installed on this system.')
        return

    print(f'\n[VERIFY] Listener snapshot for port {port}:')
    found = False
    for line in result.stdout.splitlines():
        if f':{port} ' in line:
            print(line)
            found = True
    if not found:
        print('No matching listener lines found.')


def _run_auth_step(args: argparse.Namespace) -> bool:
    """Run the interactive SSH command in the foreground."""
    if _port_open('127.0.0.1', args.local_port):
        print(
            f'\n[AUTH] Port {args.local_port} is already in use. '
            'The SSH tunnel cannot bind to it.'
        )
        _show_listeners(args.local_port)
        print(
            '[AUTH] Close the existing tunnel/process, or rerun with a '
            'different --local-port.'
        )
        return False

    cmd = _format_ssh_command(
        args.ssh_host,
        args.remote_host,
        args.remote_port,
        args.local_port,
        background=False,
        interactive=True,
        simple_ssh=args.simple_ssh,
    )
    print('\n[AUTH] Interactive auth command:')
    _show_command(cmd)
    if args.simple_ssh:
        print('[AUTH] Running in docs-style simple SSH mode.')
    print(
        '[AUTH] This command stays attached to your terminal so you can '
        'complete password / MFA prompts.'
    )
    print('[AUTH] Stop it with Ctrl+C when you are done testing.')
    subprocess.run(cmd, check=False)
    return True


def _run_verify_step(args: argparse.Namespace) -> None:
    """Run the background tunnel command and verify the port."""
    if _port_open('127.0.0.1', args.local_port):
        print(
            f'\n[VERIFY] Port {args.local_port} is already open; '
            'leaving the tunnel alone.'
        )
        _show_listeners(args.local_port)
        return

    cmd = _format_ssh_command(
        args.ssh_host,
        args.remote_host,
        args.remote_port,
        args.local_port,
        background=True,
        interactive=False,
        simple_ssh=args.simple_ssh,
    )
    if args.simple_ssh:
        print('\n[VERIFY] Running in docs-style simple SSH mode.')
    result = _run_command(cmd, 'VERIFY')
    if result.returncode != 0:
        print('[VERIFY] Tunnel start command failed.')
        return
    if not args.no_wait:
        _wait_for_port(args.local_port)
    _show_listeners(args.local_port)


def _connect_url(
    host: str,
    user: str,
    password: str,
    db_name: str,
    port: int,
):
    """Build a SQLAlchemy URL."""
    url_module = importlib.import_module('sqlalchemy.engine')
    url_cls = getattr(url_module, 'URL')

    return url_cls.create(
        'mysql+pymysql',
        username=user,
        password=password,
        host=host,
        port=port,
        database=db_name or None,
    )


def _sqlalchemy_test(args: argparse.Namespace, password: str) -> bool:
    """Run the SQLAlchemy test used by APERO RI."""
    print('\n[TEST] SQLAlchemy connection test')
    url = _connect_url(
        '127.0.0.1',
        args.db_user,
        password,
        args.db_name,
        args.local_port,
    )
    masked_url = _connect_url(
        '127.0.0.1',
        args.db_user,
        _mask_password(password),
        args.db_name,
        args.local_port,
    )
    print(f'[TEST] URL: {masked_url}')
    try:
        sqlalchemy_module = importlib.import_module('sqlalchemy')
        create_engine = getattr(sqlalchemy_module, 'create_engine')
        text = getattr(sqlalchemy_module, 'text')
    except Exception:
        print('[TEST] SQLAlchemy is not installed.')
        return False

    engine = create_engine(url, future=True, pool_pre_ping=True)
    try:
        with engine.begin() as conn:
            result = conn.execute(text('SELECT 1 AS ok'))
            row = result.first()
            print(f'[TEST] SELECT 1 result: {row}')
        print('[TEST] SQLAlchemy test succeeded.')
        return True
    except Exception:
        print('[TEST] SQLAlchemy test failed:')
        print(traceback.format_exc().rstrip())
        return False
    finally:
        engine.dispose()


def _pymysql_probe(args: argparse.Namespace, password: str) -> bool:
    """Run a direct PyMySQL connection probe."""
    print('\n[TEST] Direct PyMySQL probe')
    try:
        pymysql = importlib.import_module('pymysql')
    except Exception:
        print('[TEST] PyMySQL is not installed.')
        return False

    print(
        '[TEST] Connection details: '
        f'host=127.0.0.1 port={args.local_port} user={args.db_user} '
        f'db={args.db_name or "<none>"} password={_mask_password(password)}'
    )
    try:
        with pymysql.connect(
            host='127.0.0.1',
            user=args.db_user,
            password=password,
            database=args.db_name or None,
            port=args.local_port,
            connect_timeout=10,
            read_timeout=10,
            write_timeout=10,
            autocommit=True,
            charset='utf8mb4',
        ) as conn:
            try:
                print(f'[TEST] Server info: {conn.get_server_info()}')
            except Exception:
                print('[TEST] Could not read server info.')
            try:
                with conn.cursor() as cursor:
                    cursor.execute('SELECT 1 AS ok')
                    row = cursor.fetchone()
                    print(f'[TEST] Cursor SELECT 1 result: {row}')
            except Exception:
                print('[TEST] Query failed during direct MySQL probe:')
                print(traceback.format_exc().rstrip())
                return False
        print('[TEST] Direct PyMySQL connection probe succeeded.')
        return True
    except Exception:
        print('[TEST] Direct PyMySQL connection probe failed:')
        print(traceback.format_exc().rstrip())
        return False


def _run_test_step(args: argparse.Namespace, password: str) -> None:
    """Run the APERO-style DB test flow."""
    if not _port_open('127.0.0.1', args.local_port):
        print(
            f'\n[TEST] Local port {args.local_port} is not open. '
            'Run --verify or --auth first.'
        )
        return

    sql_ok = _sqlalchemy_test(args, password)
    pymysql_ok = _pymysql_probe(args, password)

    if sql_ok:
        print('\n[TEST] Final result: connection test passed.')
        return
    if pymysql_ok:
        print(
            '\n[TEST] Final result: SQLAlchemy query failed, but the '
            'direct MySQL probe succeeded.'
        )
        return
    print(
        '\n[TEST] Final result: both the SQLAlchemy and direct MySQL '
        'probes failed.'
    )


def main() -> int:
    """Run the requested diagnostic steps."""
    args = _parse_args()
    password_file = Path(args.password_file).expanduser()
    password = args.password or _read_password(password_file, '')

    print('Standalone DB tunnel diagnostic')
    print(f'SSH host: {args.ssh_host}')
    print(f'Remote host: {args.remote_host}')
    print(f'Remote port: {args.remote_port}')
    print(f'Local port: {args.local_port}')
    print(f'DB user: {args.db_user}')
    print(f'DB name: {args.db_name or "<none>"}')
    print(f'Password file: {password_file}')
    print(f'Password loaded: {"yes" if password else "no"}')

    if not any((args.auth, args.verify, args.test)):
        print(
            '\nNo step flags were provided, so I will run --verify then '
            '--test.'
        )
        args.verify = True
        args.test = True

    if args.auth:
        if not _run_auth_step(args):
            return 1
    if args.verify:
        _run_verify_step(args)
    if args.test:
        _run_test_step(args, password)

    return 0


# =============================================================================
# Start of code
# =============================================================================
if __name__ == '__main__':
    raise SystemExit(main())

# =============================================================================
# End of code
# =============================================================================
