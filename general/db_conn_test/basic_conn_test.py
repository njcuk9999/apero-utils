#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Minimal SQLAlchemy MySQL connection test.

This does not require the mysql CLI. It only uses Python + SQLAlchemy.
"""

import importlib
from pathlib import Path

# =============================================================================
# Define variables
# =============================================================================
DEFAULT_HOST = '127.0.0.1'
DEFAULT_PORT = 3307
DEFAULT_DB_USER = 'lmalo'
DEFAULT_DB_NAME = 'lmalo_spirou'
try:
    DEFAULT_PASSWORD_FILE = Path(__file__).with_name('pass.txt')
except Exception:
    DEFAULT_PASSWORD_FILE = Path('pass.txt')


# =============================================================================
# Define functions
# =============================================================================
def read_password(password_file: Path) -> str:
	"""Return the first non-empty line from a password file."""
	if not password_file.exists():
		return ''
	raw = password_file.read_text(encoding='utf-8', errors='replace')
	for line in raw.splitlines():
		value = line.strip()
		if value:
			return value
	return ''


def main() -> int:
	"""Test a MySQL connection through SQLAlchemy."""
	try:
		sqlalchemy = importlib.import_module('sqlalchemy')
		sql_engine = importlib.import_module('sqlalchemy.engine')
		url_cls = getattr(sql_engine, 'URL')
	except Exception:
		print('Missing dependency: install sqlalchemy and pymysql first.')
		return 1

	password = read_password(DEFAULT_PASSWORD_FILE)
	if not password:
		print(f'No password found in {DEFAULT_PASSWORD_FILE}')
		return 1

	url = url_cls.create(
		'mysql+pymysql',
		username=DEFAULT_DB_USER,
		password=password,
		host=DEFAULT_HOST,
		port=DEFAULT_PORT,
		database=DEFAULT_DB_NAME,
	)
	engine = sqlalchemy.create_engine(url, future=True, pool_pre_ping=True)

	try:
		with engine.connect() as conn:
			result = conn.execute(sqlalchemy.text('SELECT 1 AS ok'))
			row = result.first()
			print('Connection OK')
			print(f'SELECT 1 result: {row}')
		return 0
	except Exception as exc:
		print('Connection FAILED')
		print(exc)
		return 1
	finally:
		engine.dispose()


# =============================================================================
# Start of code
# =============================================================================
if __name__ == '__main__':
	raise SystemExit(main())

# =============================================================================
# End of code
# =============================================================================

