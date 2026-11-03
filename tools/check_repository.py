#!/usr/bin/env python3
"""Check code, documentation and privacy without running analysis modules."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess


ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN = {'.bam', '.bai', '.cram', '.crai', '.fastq', '.fq', '.vcf', '.bcf',
             '.maf', '.tsv', '.csv', '.xlsx', '.xls', '.rds', '.rdata', '.h5',
             '.hdf5', '.mtx', '.h5ad', '.pdf', '.png', '.jpg', '.jpeg', '.svg',
             '.tif', '.tiff', '.doc', '.docx', '.pptx', '.zip', '.gz', '.pyc', '.html'}
PRIVATE_PATH = re.compile('/' + r'(?:Users|media|cbio/users|home)/[^/\s]+/')
SAMPLE_LITERAL = re.compile(r'\b\d{1,3}[TN](?:_(?:WES|RNA))?\b')
PATIENT_LIST = re.compile(r'\bpatients?\s+\d{1,3}\s*(?:,|and\b|;|\.)', re.I)
EDITORIAL_CONTEXT = re.compile(r'collaborat(?:or|ors)|reorganiz(?:ed|ation)|superseded|legacy\s+\d+|earlier STAR|previous manuscript', re.I)
LINK = re.compile(r'\[[^\]]+\]\(([^)]+)\)')


def files():
    listing = subprocess.check_output(
        ['git', 'ls-files', '--cached', '--others', '--exclude-standard'], cwd=ROOT, text=True
    )
    return [ROOT / name for name in listing.splitlines() if (ROOT / name).is_file()]


def main():
    problems = []
    paths = files()
    parsed = 0
    config = json.loads((ROOT / 'config/private.example.json').read_text())
    for name, values in config['cohorts'].items():
        if values:
            problems.append(f'Example contains real/nonempty cohort membership: {name}')
    for path in paths:
        relative = path.relative_to(ROOT)
        if path.is_symlink():
            problems.append(f'Symlink not allowed: {relative}')
        if set(s.lower() for s in path.suffixes) & FORBIDDEN:
            problems.append(f'Data/result/binary file not allowed: {relative}')
        try:
            text = path.read_text(encoding='utf-8')
        except UnicodeDecodeError:
            problems.append(f'Non-text file: {relative}')
            continue
        if PRIVATE_PATH.search(text):
            problems.append(f'Private absolute path: {relative}')
        if SAMPLE_LITERAL.search(text):
            problems.append(f'Embedded specimen identifier: {relative}')
        if PATIENT_LIST.search(text):
            problems.append(f'Embedded patient list: {relative}')
        if relative.as_posix() != 'tools/check_repository.py' and EDITORIAL_CONTEXT.search(text):
            problems.append(f'Editorial context in project documentation: {relative}')
        if path.suffix == '.py':
            try:
                tree = ast.parse(text, filename=str(relative))
                parsed += 1
                for node in ast.walk(tree):
                    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                        if node.func.id == 'configured_cohort' and node.args and isinstance(node.args[0], ast.Constant):
                            if node.args[0].value not in config['cohorts']:
                                problems.append(f'Unknown private cohort key in {relative}')
            except SyntaxError as error:
                problems.append(f'Syntax error in {relative}: {error}')
        if path.suffix == '.md':
            for target in LINK.findall(text):
                target = target.split('#', 1)[0]
                if not target or re.match(r'^(?:https?://|mailto:)', target):
                    continue
                if not (path.parent / target).exists():
                    problems.append(f'Broken documentation link in {relative}: {target}')
    manifest = json.loads((ROOT / 'docs/source_manifest.json').read_text())
    for entry in manifest['files']:
        target = ROOT / entry['path']
        if not target.is_file():
            problems.append(f'Missing code file: {entry["path"]}')
        elif hashlib.sha256(target.read_bytes()).hexdigest() != entry['sha256']:
            problems.append(f'Update published source checksum: {entry["path"]}')
    if problems:
        for problem in problems:
            print('FAIL:', problem)
        return 1
    print(f'PASS: {len(paths)} text files; {parsed} Python files parsed; {len(manifest["files"])} code hashes checked.')
    print('No analysis code was imported or executed.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
