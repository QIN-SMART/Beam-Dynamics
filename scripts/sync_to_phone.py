#!/usr/bin/env python3
"""
Sync modified project files to iCloud Drive for ChatGPT on the phone.

DESIGN (phone-friendly, one file at a time):
  - DEFAULT = BUNDLE ONLY: one merged markdown `latest_bundle.txt` with the
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
  python3 scripts/sync_to_phone.py --gpt-review # upload the A-level GPT review
                                                # package (core code + latest
                                                # report/figure/changelog) to
                                                # UED_Sync/gpt_review/

iPhone: Files app → iCloud Drive → UED_Sync/ → latest_bundle.txt
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
def _txt(name):
    """iCloud-side files use .txt (content unchanged) — local .md untouched."""
    return name[:-3] + ".txt" if name.lower().endswith(".md") else name


SYNC_DIR_NAME = "UED_Sync"
MANIFEST_NAME = "manifest.json"
SUMMARY_NAME = "sync_summary.md"
BUNDLE_NAME = "latest_bundle.txt"
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


# ══════════════════════════════════════════════════════════════════════
#  GPT-review package — the A-level files most worth sending to GPT
# ══════════════════════════════════════════════════════════════════════

# Always-shipped core files (highest data-flow / physics risk, per the
# project review spec).  Flattened names used in the iCloud folder.
GPT_REVIEW_CORE = [
    "validation/backend.py",
    "shared/beamline_config.yaml",
    "shared/params.py",
    # v0.13 contract layer (single-source derivation & coordinate access)
    "shared/beam_physics.py",
    "shared/ocelot_coords.py",
    "shared/constants.py",
    "validation/config_check.py",
    "validation/test_rf.py",
    "validation/test_full_beamline.py",
    "validation/test_gpt_route_equivalence.py",
    # v0.14.1 task 1: SC scheduler characterization + production acceptance
    "validation/test_sc_scheduler_equivalence.py",
    "validation/run_all.py",
    "GPT模拟/ued_beamline_v2.py",
    "AG/run_shared.py",
    "AGENTS.md",
    "validation/CHECKPOINTS.md",
]

# A-level curated files (flattened names) archived per run — the "most
# worth sending to GPT" set, kept lean to avoid duplication.
GPT_REVIEW_A_LEVEL = [
    "core_validation_backend.py",
    "core_GPT模拟_ued_beamline_v2.py",
    "core_shared_beam_physics.py",
    "core_shared_ocelot_coords.py",
    "core_validation_test_rf.py",
    "core_validation_test_full_beamline.py",
    "core_validation_test_gpt_route_equivalence.py",
    "core_validation_test_sc_scheduler_equivalence.py",
    "latest_report.txt",
    "latest_review.png",
    "latest_changelog.txt",
    "latest_manifest.txt",
    "sc_ocelot_source_audit.txt",
    "sc_ag_model_audit.txt",
    "sc_diagnostics.png",
]

# Auto-picked latest artifacts (one per category, by mtime).
GPT_REVIEW_LATEST = [
    ("validation/reports/*.md", "latest_report.txt"),
    ("validation/reports/review_summary*.png", "latest_review.png"),
    ("validation/reports/SC_ocelot_source_audit.md", "sc_ocelot_source_audit.txt"),
    ("validation/reports/SC_AG_model_audit.md", "sc_ag_model_audit.txt"),
    ("validation/reports/SC_diagnostics.png", "sc_diagnostics.png"),
    ("CHANGELOG_*.md", "latest_changelog.txt"),
    ("validation/baselines/*/BASELINE_MANIFEST.md", "latest_manifest.txt"),
]


def glob_md(directory):
    import glob
    return glob.glob(os.path.join(directory, "*.md"))


def _newest(glob_pattern):
    import glob
    cands = [p for p in glob.glob(os.path.join(PROJECT_ROOT, glob_pattern))
             if os.path.isfile(p)]
    if not cands:
        return None
    return max(cands, key=os.path.getmtime)


def sync_gpt_review(sync_dir):
    """Copy the A-level GPT-review package into iCloud/UED_Sync/gpt_review/."""
    dest = os.path.join(sync_dir, "gpt_review")
    os.makedirs(dest, exist_ok=True)
    # remove stale .md leftovers (iCloud side is .txt only)
    for old in glob_md(dest):
        os.remove(old)
    copied = []          # (flattened_name, original_relpath, why)

    for rel in GPT_REVIEW_CORE:
        src = os.path.join(PROJECT_ROOT, rel)
        if not os.path.isfile(src):
            continue
        flat = _txt("core_" + rel.replace("/", "_"))
        shutil.copy2(src, os.path.join(dest, flat))
        copied.append((flat, rel, "core"))

    for glob_pattern, out_name in GPT_REVIEW_LATEST:
        latest = _newest(glob_pattern)
        if latest is None:
            continue
        rel = os.path.relpath(latest, PROJECT_ROOT)
        shutil.copy2(latest, os.path.join(dest, _txt(out_name)))
        copied.append((_txt(out_name), rel, "latest"))

    # A/B/C level map (matches the project review spec)
    readme = [
        "# GPT-review package — what to send and why",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "Folder: iCloud Drive → UED_Sync → gpt_review/",
        "",
        "## A-level (strongly recommended — high-risk, hard to verify by text)",
        "",
        "- core_validation_backend.py  : dual-backend adapter, RF/R56 conversions,",
        "                                element assembly, switches.",
        "- core_GPT模拟_ued_beamline_v2.py : main route (lattice single source).",
        "- core_shared_beam_physics.py : single γ/β/p0 derivation (v0.13).",
        "- core_shared_ocelot_coords.py: semantic rparticles access (v0.13).",
        "- core_shared_constants.py    : single physical-constants source.",
        "- core_validation_test_rf.py / test_full_beamline.py /",
        "  test_gpt_route_equivalence.py : acceptance tests.",
        "- latest_review.png           : one merged figure with the key curves.",
        "",
        "## B-level (send when something looks wrong)",
        "",
        "- core_shared_beamline_config.yaml : when suspecting parameter source.",
        "- core_shared_params.py / AG_run_shared.py / core_validation_config_check.py.",
        "- latest_report.txt / latest_changelog.txt.",
        "",
        "## C-level (archive only)",
        "",
        "- core_AGENTS.txt, core_validation_CHECKPOINTS.txt : provenance/history,",
        "  not direct evidence.",
        "- latest_manifest.txt (baseline snapshot).",
        "",
        "## Copied in this run",
        "",
    ]
    for flat, rel, why in copied:
        readme.append(f"- {flat}   ← {rel}   [{why}]")
    with open(os.path.join(dest, "REVIEW_LIST.txt"), "w") as f:
        f.write("\n".join(readme) + "\n")

    # ── version snapshot: A-level files only (no duplicates) ──
    ts = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    vdir = os.path.join(dest, "versions", ts)
    os.makedirs(vdir, exist_ok=True)
    a_copied = []
    for name in os.listdir(dest):
        if name == "versions" or not os.path.isfile(os.path.join(dest, name)):
            continue
        if name in GPT_REVIEW_A_LEVEL:
            shutil.copy2(os.path.join(dest, name), os.path.join(vdir, name))
            a_copied.append(name)
    with open(os.path.join(vdir, "info.txt"), "w") as f:
        f.write(f"A-level snapshot at {datetime.now():%Y-%m-%d %H:%M:%S}\n"
                f"A-level files: {len(a_copied)}\n"
                f"This is the curated 'most worth sending to GPT' set; "
                f"gpt_review/ root holds the full latest package.\n")
    with open(os.path.join(dest, "VERSIONS.txt"), "a") as f:
        f.write(f"{ts}  ({len(a_copied)} A-level files)\n")

    print(f"GPT-review package -> {dest}")
    print(f"  files copied: {len(copied)}  (list: gpt_review/REVIEW_LIST.txt)")
    print(f"  snapshot     -> {vdir}")
    return dest


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--files", action="store_true",
                    help="also copy individual changed files into files/")
    ap.add_argument("--images", action="store_true",
                    help="include images/binaries in the individual copies")
    ap.add_argument("--no-bundle", action="store_true",
                    help="skip the merged latest_bundle.txt")
    ap.add_argument("--gpt-review", action="store_true",
                    help="upload the A-level GPT-review package to gpt_review/")
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

    if args.gpt_review:
        sync_gpt_review(sync_dir)
        return 0
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
            dst = os.path.join(sync_dir, FILES_SUBDIR, _txt(rel))
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
    print("iPhone: Files app → iCloud Drive → UED_Sync/ → latest_bundle.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
