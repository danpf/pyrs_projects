#!/usr/bin/env python
"""
For now we'll run most with shell=True so that we can move away from this quickly if necessary
"""

import argparse
import os
import platform
from typing import List
import sys
import logging
import subprocess

LOGLEVEL = os.environ.get("LOGLEVEL", "INFO").upper()
logging.basicConfig(level=LOGLEVEL)
log = logging.getLogger(__name__)

CIBUILDWHEEL_VERSION = "2.16.2"
CIBUILDWHEEL_PLATFORM = "auto"


def parseargs(args: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=["install_ci_deps", "cicd_dpf_ssw_aligner_rs", "cicd_dpf_mrcfile"],
        help="mode to run",
        required=True,
    )
    parser.add_argument(
        "--cibw_platform",
        choices=["auto", "linux", "macos", "windows"],
        help="cibuildwheel platform to build",
        default=CIBUILDWHEEL_PLATFORM,
    )
    return parser.parse_args(args)


def run_l(cmd: list[str]) -> None:
    str_cmd = " ".join(cmd)
    log.debug(f"Running: {str_cmd}")
    ret = subprocess.run(cmd)
    if ret.returncode:
        raise RuntimeError(f"Failure with command: {str_cmd}")


def run_s(cmd: str, env_add: dict[str, str]) -> None:
    log.debug(f"Running: {cmd}")
    env = dict(os.environ)
    env.update(env_add)
    ret = subprocess.run(cmd, env=env, shell=True)
    if ret.returncode:
        raise RuntimeError(f"Failure with command: {cmd}")


def install_rust_cibw_command():
    return (
        "curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs"
        " | sh -s -- --default-toolchain none -y"
        " && rustup toolchain install nightly --allow-downgrade --profile minimal --component clippy"
        " && rustup show"
    )


def install_ci_deps():
    cmds = f"""
python -m pip install cibuildwheel=={CIBUILDWHEEL_VERSION}
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- --default-toolchain none -y
rustup toolchain install nightly --allow-downgrade --profile minimal --component clippy
rustup show
pip install pipx
"""
    match platform.system():
        case "Darwin":
            cmds += "rustup target add aarch64-apple-darwin x86_64-apple-darwin"
        case "Windows":
            cmds += "rustup target add x86_64-pc-windows-msvc i686-pc-windows-msvc"
    for cmd in cmds.split("\n"):
        if cmd:
            run_s(cmd, {})


def cicd_dpf_ssw_aligner_rs():
    env = {
        "CIBW_PLATFORM": CIBUILDWHEEL_PLATFORM,
        "CIBW_ARCHS": "auto",
        "CIBW_ENVIRONMENT": 'PATH="$HOME/.cargo/bin:$PATH" CARGO_TERM_COLOR="always"',
        "CIBW_ENVIRONMENT_WINDOWS": 'PATH="$UserProfile\\.cargo\\bin;$PATH"',
        "CIBW_BEFORE_BUILD": install_rust_cibw_command(),
        "CIBW_BUILD_VERBOSITY": "1",
    }

    cmds = """
python -m cibuildwheel ./src/python/dpf-ssw-aligner-rs --output-dir wheelhouse
ls wheelhouse
pipx run build --sdist --outdir wheelhouse ./src/python/dpf-ssw-aligner-rs
ls wheelhouse
    """
    for cmd in cmds.split("\n"):
        if cmd:
            run_s(cmd, env)


def cicd_dpf_mrcfile():
    env = {
        "CIBW_PLATFORM": CIBUILDWHEEL_PLATFORM,
        "CIBW_ARCHS": "auto",
        "CIBW_ENVIRONMENT": 'PATH="$HOME/.cargo/bin:$PATH" CARGO_TERM_COLOR="always"',
        "CIBW_ENVIRONMENT_WINDOWS": 'PATH="$UserProfile\\.cargo\\bin;$PATH"',
        "CIBW_BEFORE_BUILD": install_rust_cibw_command(),
        "CIBW_BUILD_VERBOSITY": "1",
    }

    cmds = """
python -m cibuildwheel ./src/python/dpf-mrcfile --output-dir wheelhouse
ls wheelhouse
pipx run build --sdist --outdir wheelhouse ./src/python/dpf-mrcfile
ls wheelhouse
    """
    for cmd in cmds.split("\n"):
        if cmd:
            run_s(cmd, env)


def main(_args: List[str]) -> None:
    args = parseargs(_args)
    global CIBUILDWHEEL_PLATFORM
    CIBUILDWHEEL_PLATFORM = args.cibw_platform
    match args.mode:
        case "install_ci_deps":
            install_ci_deps()
        case "cicd_dpf_ssw_aligner_rs":
            install_ci_deps()
            cicd_dpf_ssw_aligner_rs()
        case "cicd_dpf_mrcfile":
            install_ci_deps()
            cicd_dpf_mrcfile()


if __name__ == "__main__":
    main(sys.argv[1:])
