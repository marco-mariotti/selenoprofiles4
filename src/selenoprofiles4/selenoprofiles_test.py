#!/usr/bin/python -u
from string import *
import os, sys, traceback, subprocess
from easyterm import command_line_options, write
import importlib.resources
import tempfile
import shutil

help_msg="""selenoprofiles test: utility to verify that the installation works correctly.

Usage:   selenoprofiles test [-b path_to_selenoprofiles] [-i input_fasta] [-o output_folder]

Description:
  This utility runs a series of quick tests to verify that the selenoprofiles
  installation is functional. It runs several families (MSRB, DIO, SEPHS, GPX)
  against a small test FASTA dataset and checks that expected outputs are produced.

Options:
    -b       Path to the Selenoprofiles executable (optional). If not provided,
           it is auto-detected using 'which selenoprofiles' or 'which selenoprofiles4.py'.
    -i       Input FASTA file containing test sequences (required).
    -o       Output folder where results will be stored (required).
    -h, --help     Show this help message and exit
"""

def get_default_test_fasta():
    try:
        # Use importlib to safely find the installed data file
        return importlib.resources.files("selenoprofiles4.tests") / "test_sequences.fa"
    except Exception:
        return "test_sequences.fa"

# Default options
def_opt = {
    "cmd": "test"
}

class NotracebackException(Exception):
    """Cleaner output exception"""
    pass

def bash(cmd, print_it=False):
    """Run a bash command and return (exit_status, stdout+stderr)."""
    if print_it:
        write(f"$ {cmd}", end="\n")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.returncode, result.stdout + result.stderr

#########################################################
###### start main program function

def run_single_test(family, selenoprofiles_bin, results_folder, sequence_file, logname, success_marker):
    """Run a single selenoprofiles family test and verify output."""
    write(f"\n## Test: {family}", end="\n")
    cmd = (
        f"{selenoprofiles_bin} -o {results_folder} -t {sequence_file} "
        f'-s "Homo sapiens" -P {family} -B -output_p2g > {logname} 2>&1'
    )
    status, output = bash(cmd, print_it=True)
    if status != 0:
        raise NotracebackException(f"ERROR: test {family} failed. See {logname}\n\n{output}")

    # Check for errors or expected outputs
    e, _ = bash(f"grep ERROR {logname}")
    if e == 0:
        raise NotracebackException(f"ERROR found in log: {logname}")
    _, warn = bash(f"grep WARNING {logname}")
    if warn.strip():
        write(f"-- WARNING:\n{warn}", end="\n")
    s, _ = bash(f'grep P2G {logname} | grep "{success_marker}"')
    if s != 0:
        raise NotracebackException(f"-- ERROR! Expected '{success_marker}' not found in output.")
    write("##### OK!\n", end="\n")

def main(opt=None):

    # If no options were passed, use defaults
    if not opt:
        opt = def_opt.copy()

    write("Looking for Selenoprofiles executable...", end="\n")

    # Try autodetect
    for candidate in ["selenoprofiles", "selenoprofiles4.py"]:
        status, path = bash(f"which {candidate}")
        if status == 0 and os.path.isfile(path.strip()):
            selenoprofiles_bin = path.strip()
            break
    if not selenoprofiles_bin:
        raise NotracebackException("ERROR: selenoprofiles executable not found. Use -b to specify it manually.")

    # Paths
    sequence_file = str(get_default_test_fasta())

    if not os.path.isfile(sequence_file):
        raise NotracebackException(
            "ERROR: cannot run tests. test_sequences.fa not found! Run this in the selenoprofiles installation directory."
        )
    
    # Create temporary folder if no output folder is given
    temp_folder = None
    temp_folder = tempfile.mkdtemp(prefix="selenoprofiles_test_")
    results_folder = temp_folder

    write("Starting tests...\n", end="\n")

    try:
        # Run all four tests
        run_single_test("MSRB", selenoprofiles_bin, results_folder, sequence_file, "log_test1", "cysteine")
        run_single_test("DIO", selenoprofiles_bin, results_folder, sequence_file, "log_test2", "selenocysteine")
        run_single_test("SEPHS2", selenoprofiles_bin, results_folder, sequence_file, "log_test3", "threonine")
        run_single_test("GPX", selenoprofiles_bin, results_folder, sequence_file, "log_test4", "cysteine")

        write("\nAll tests completed successfully.\n", end="\n")
    finally:
        # Clean up temporary folder if created
        if temp_folder:
            shutil.rmtree(temp_folder)
            write(f"Temporary folder {temp_folder} removed.\n", end="\n")


"""
if __name__ == "__main__":
    try:
        main()
    except NotracebackException as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(1)
    except Exception:
        traceback.print_exc()
        sys.exit(1)
"""