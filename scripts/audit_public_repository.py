#!/usr/bin/env python3
"""Fail if the public repository contains private data or unsafe references."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_SUFFIXES = {
    ".bam", ".bai", ".cram", ".crai", ".fastq", ".fq", ".vcf", ".bcf",
    ".maf", ".xlsx", ".xls", ".rds", ".rdata", ".h5", ".h5ad", ".mtx",
    ".pdf", ".png", ".jpg", ".jpeg", ".tif", ".tiff", ".html", ".pyc",
}
PRIVATE_PATH_PATTERNS = (
    re.compile("/" + r"Users/[^/]+/"),
    re.compile("/" + r"media/[^/]+/"),
    re.compile("/" + r"home/(?!EXAMPLE_USER\b)[^/]+/"),
    re.compile("New_" + "Volume3|New " + "Volume3"),
)
REALISTIC_SAMPLE_PATTERNS = (
    re.compile(r"\bPATIENT\d+\b", re.I),
    re.compile(r"\b\d{1,3}[TN]_(?:WES|RNA)\b", re.I),
    re.compile(r"\bTumou?r\s+\d{2}\b", re.I),
)
MARKDOWN_LINK = re.compile(r"(?<!!)\[[^]]+\]\(([^)]+)\)")


def repository_files() -> list[Path]:
    try:
        output = subprocess.check_output(
            ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
            cwd=ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return sorted(
            path for path in ROOT.rglob("*")
            if path.is_file() and ".git" not in path.relative_to(ROOT).parts
        )
    return [ROOT / line for line in output.splitlines() if line and (ROOT / line).is_file()]


def audit_file(path: Path) -> list[str]:
    relative = path.relative_to(ROOT)
    problems: list[str] = []
    if path.is_symlink():
        return [f"unexpected symbolic link: {relative}"]
    lower_name = str(relative).lower()
    if any(word in lower_name for word in ("manuscript", "reviewer_comment", "review_comments")):
        problems.append(f"disallowed publication document name: {relative}")
    suffixes = {suffix.lower() for suffix in path.suffixes}
    if suffixes.intersection(FORBIDDEN_SUFFIXES):
        problems.append(f"forbidden patient-data/result extension: {relative}")
    if path.stat().st_size > 5 * 1024 * 1024:
        problems.append(f"file exceeds 5 MB code/documentation limit: {relative}")

    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        problems.append(f"unexpected binary file: {relative}")
        return problems

    for pattern in PRIVATE_PATH_PATTERNS:
        if pattern.search(text):
            problems.append(f"private or retired absolute path in {relative}: {pattern.pattern}")
    for pattern in REALISTIC_SAMPLE_PATTERNS:
        if pattern.search(text):
            problems.append(f"realistic patient/sample identifier in {relative}: {pattern.pattern}")

    if path.suffix.lower() == ".md":
        for target in MARKDOWN_LINK.findall(text):
            target = target.strip().split("#", 1)[0]
            if not target or re.match(r"^(?:https?://|mailto:)", target):
                continue
            resolved = (path.parent / target).resolve()
            if not resolved.exists():
                problems.append(f"broken local Markdown link in {relative}: {target}")
    return problems


def main() -> int:
    problems = []
    files = repository_files()
    for path in files:
        problems.extend(audit_file(path))
    if problems:
        print("PUBLIC REPOSITORY AUDIT FAILED")
        for problem in problems:
            print(f"- {problem}")
        return 1
    print(f"PUBLIC REPOSITORY AUDIT PASSED: {len(files)} files checked")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
