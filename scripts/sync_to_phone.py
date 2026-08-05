#!/usr/bin/env python3
"""
Sync modified project files to iCloud Drive for ChatGPT on the phone.

DESIGN (phone-friendly, one file at a time):
  - DEFAULT = BUNDLE ONLY: one merged markdown `latest_bundle.md` with the
    full content of every file changed since the last run — a single file to
    open on the phone and feed to ChatGPT.
  - Regenerable data (e.g. *_results.json) is never synced — it carries no
    context for ChatGPT.
  - `--files` opt-in: additionally copy the individual changed files into
    files/ (incremental, sha256 manifest).
  - `--images` opt-in: include images/binaries in the individual copies.

Usage
-----
  python3 scripts/sync_to_phone.py              # bundle of changed files
  python3 scripts/sync_to_phone.py --dry-run    # show what would change
  python3 scripts/sync_to_phone.py --force      # treat all files as changed
  python3 scripts/sync_to_phone.py --files      # also copy individual files
  python3 scripts/sync_to_phone.py --images     # ... including images

iPhone: Files app → iCloud Drive → UED_Sync/ → latest_bundle.md
"""

import argparse
import hashlib
import json
import os
import shutil
import sys
from datetime import datetime

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ICLOUD_ROOT = os.path.expanduser(
    "~/Library/Mobile Documents/com~apple~CloudDocs")
SYNC_DIR_NAME = "UED_Sync"
MANIFEST_NAME = "manifest.json"
SUMMARY_NAME = "sync_summary.md"
BUNDLE_NAME = "latest_bundle.md"
FILES_SUBDIR = "files"

# Everything outside these is treated as live project content.
DEFAULT_EXCLUDE_DIRS = {
    ".git", "__pycache__", ".vscode", ".DS_Store", "scripts",
    "Okemos1.22b_imagSystem",     # legacy MFC C++ control code
    "AG-prime",                   # legacy AG version
    "notes-codes", "codes-result", "code-update-log", "log",  # historical
    "文献",                        # reference literature (PDFs, not code)
}
DEFAULT_EXCLUDE_FILES = {".DS_Store"}
TEXT_EXTS = {".py", ".yaml", ".yml", ".md", ".tex", ".json", ".csv",
             ".sh", ".txt", ".cfg", ".toml"}
MAX_FILE_MB = 20.0


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def is_text(name, include_images):
    if include_images:
        return True
    return os.path.splitext(name)[1].lower() in TEXT_EXTS


def is_regenerable(rel):
    """Data artifacts with no context value for ChatGPT — never synced."""
    return rel.endswith("_results.json") or "/results/" in rel


def walk_project(exclude_dirs, exclude_files, max_bytes, include_images):
    out = []
    for root, dirs, files in os.walk(PROJECT_ROOT):
        dirs[:] = [d for d in dirs if d not in exclude_dirs]
        for name in files:
            if name in exclude_files:
                continue
            if not is_text(name, include_images):
                continue
            rel = os.path.relpath(os.path.join(root, name), PROJECT_ROOT)
            if is_regenerable(rel):
                continue
            full = os.path.join(root, name)
            try:
                size = os.path.getsize(full)
            except OSError:
                continue
            if size > max_bytes:
                continue
            out.append((rel, full, size))
    return out


def load_manifest(path):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def save_manifest(path, manifest):
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)


def build_bundle(rels, changed_paths):
    """One markdown file with the full content of the changed files."""
    lines = ["# UED project — changed files (bundle)",
             "",
             f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
             f"Files in this bundle: {len(rels)}",
             ""]
    for rel in rels:
        full = changed_paths[rel]
        lines += [f"## FILE: {rel}", ""]
        try:
            with open(full, encoding="utf-8", errors="replace") as f:
                content = f.read()
        except OSError as e:
            content = f"<unreadable: {e}>"
        lines.append(content)
        lines.append("")
        lines.append("---")
        lines.append("")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--files", action="store_true",
                    help="also copy individual changed files into files/")
    ap.add_argument("--images", action="store_true",
                    help="include images/binaries in the individual copies")
    ap.add_argument("--no-bundle", action="store_true",
                    help="skip the merged latest_bundle.md")
    ap.add_argument("--max-mb", type=float, default=MAX_FILE_MB)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.out:
        sync_dir = os.path.abspath(args.out)
    else:
        if not os.path.isdir(ICLOUD_ROOT):
            print(f"iCloud Drive not found at {ICLOUD_ROOT}")
            sys.exit(1)
        sync_dir = os.path.join(ICLOUD_ROOT, SYNC_DIR_NAME)
    os.makedirs(os.path.join(sync_dir, FILES_SUBDIR), exist_ok=True)

    manifest_path = os.path.join(sync_dir, MANIFEST_NAME)
    manifest = {} if args.force else load_manifest(manifest_path)

    files = walk_project(DEFAULT_EXCLUDE_DIRS, DEFAULT_EXCLUDE_FILES,
                         args.max_mb * 1e6, args.images)
    hashes, changed_paths = {}, {}
    added = changed = 0
    for rel, full, size in files:
        h = sha256(full)
        hashes[rel] = h
        if rel not in manifest:
            added += 1
            changed_paths[rel] = full
        elif manifest[rel] != h:
            changed += 1
            changed_paths[rel] = full
    removed = [rel for rel in manifest if rel not in hashes]

    print(f"iCloud sync dir : {sync_dir}")
    print(f"scanned (text)  : {len(files)}  "
          f"(images included: {args.images}, results-JSON excluded)")
    print(f"added           : {added}")
    print(f"changed         : {changed}")
    print(f"removed         : {len(removed)}")

    if args.dry_run:
        for rel in sorted(changed_paths):
            print(f"  [bundle] {rel}")
        return 0

    if args.files:
        for rel in sorted(changed_paths):
            src = os.path.join(PROJECT_ROOT, rel)
            dst = os.path.join(sync_dir, FILES_SUBDIR, rel)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(src, dst)
        for rel in removed:
            try:
                os.remove(os.path.join(sync_dir, FILES_SUBDIR, rel))
            except OSError:
                pass

    save_manifest(manifest_path, hashes)

    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    summary = [f"# UED project sync — {now}",
               "",
               f"- Added: {added}, changed: {changed}, removed: {len(removed)}.",
               f"- Files synced to iCloud Drive → {SYNC_DIR_NAME}/files/",
               "",
               "## Changed in this run", ""]
    for rel in sorted(changed_paths):
        summary.append(f"- `{rel}`")
    if not changed_paths:
        summary.append("- (nothing changed since last sync)")
    with open(os.path.join(sync_dir, SUMMARY_NAME), "w") as f:
        f.write("\n".join(summary) + "\n")

    if not args.no_bundle and changed_paths:
        bundle = build_bundle(sorted(changed_paths), changed_paths)
        bundle_path = os.path.join(sync_dir, BUNDLE_NAME)
        with open(bundle_path, "w", encoding="utf-8") as f:
            f.write(bundle)
        bsize = os.path.getsize(bundle_path)
        if bsize > 2e6:
            print(f"WARNING: bundle is {bsize/1e6:.1f} MB — consider --dry-run "
                  f"or editing which files change")

    print(f"\nsummary -> {os.path.join(sync_dir, SUMMARY_NAME)}")
    if not args.no_bundle:
        print(f"bundle  -> {os.path.join(sync_dir, BUNDLE_NAME)}"
              + ("" if changed_paths else " (no changes, kept previous)"))
    print("iPhone: Files app → iCloud Drive → UED_Sync/ → latest_bundle.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
