#!/usr/bin/env python3

import glob, io, json, os, platform, stat, sys, zipfile
from urllib.parse import quote_plus
from urllib.request import urlopen
import datetime, logging

def main():
    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    logging.basicConfig(level = logging.DEBUG, format = "%(asctime)s - %(message)s", handlers = [logging.FileHandler(f"{timestamp}_fragpipe_install.log"), logging.StreamHandler()])

    system_platform = platform.platform(aliased = True, terse = True)
    if system_platform.lower().find("linux") == -1:
        logging.error(f"Error: Unsupported system platform '{system_platform}'.")
        if system_platform.lower().find("windows") != -1:
            logging.error("Although FragPipe supports Windows, GenomeProt has not been tested on it yet.")
        else:
            logging.error("Neither FragPipe nor GenomeProt support this platform.")
        logging.error("Stopping the installation of FragPipe and its tools...")
        sys.exit(1)

    # Default: Extract FragPipe and its tools in the current directory
    # Users can specify a different extraction directory by providing an argument
    if not 1 <= len(sys.argv) <= 2 or (len(sys.argv) == 2 and sys.argv[1].strip().lower() in ["-h", "--help"]):
        logging.info("Usage: python fragpipe_installer.py (extraction_dir)")
        sys.exit(1)

    extraction_dir = os.getcwd()
    if len(sys.argv) == 2:
        try:
            extraction_dir = os.path.realpath(sys.argv[1])
            if not os.path.exists(extraction_dir):
                os.makedirs(extraction_dir)
                os.chdir(extraction_dir)
        except KeyboardInterrupt:
            logging.error("User has quit the installation script. Stopping the installation of FragPipe and its tools...")
            sys.exit(1)
        except Exception as e:
            logging.error(f"Error: An exception has occurred when attempting to create and navigate to the user-specified extraction directory '{sys.argv[1]}'.")
            logging.error(e)
            logging.error("Stopping the installation of FragPipe and its tools...")
            sys.exit(1)
    else:
        logging.info("No arguments provided. Installing FragPipe and its tools in the current directory...")

    logging.info(f"Extraction directory: {extraction_dir}")

    while True:
        try:
            resp = input("Would you like to install FragPipe ([Y]/n)? ").strip().lower()
            if resp.startswith('n'):
                logging.error("Stopping the installation of FragPipe and its tools...")
                sys.exit(1)
            elif resp.startswith('y') or resp == '':
                logging.info("Beginning installation of FragPipe...")
                break
            logging.error(f"Error: Invalid response '{resp}'. Please enter 'y' or 'n' (case-insensitive).")
            continue
        except KeyboardInterrupt:
            logging.error("User has quit the installation script. Stopping the installation of FragPipe and its tools...")
            sys.exit(1)

    fragpipe_download_url = "https://github.com/Nesvilab/FragPipe/releases/download/23.1/FragPipe-23.1-linux.zip"
    logging.info(f"Downloading FragPipe v23.1 from {fragpipe_download_url}...")

    zip_magic_number = b"PK\x03\x04"

    while True:
        try:
            with urlopen(fragpipe_download_url) as f:
                fragpipe_zip = f.read()
            if not fragpipe_zip.startswith(zip_magic_number):
                logging.error("Error: The FragPipe zip file does not appear to be a zip file. Retrying...")
                continue
            logging.info("FragPipe downloaded!")
            break
        except KeyboardInterrupt:
            logging.error("User has quit the installation script. Stopping the installation of FragPipe and its tools...")
            sys.exit(1)
        except Exception as e:
            logging.error("Error: An exception has occurred when downloading FragPipe.")
            logging.error(e)
            logging.error("Retrying...")
            continue

    # Extract FragPipe
    logging.info("Unzipping FragPipe...")

    with io.BytesIO(fragpipe_zip) as fragpipe_zip_stream:
        with zipfile.ZipFile(fragpipe_zip_stream, 'r') as fragpipe_zip_file:
            fragpipe_zip_file.extractall(extraction_dir)

    del fragpipe_zip

    logging.info("Done!")

    # Set execute permissions for FragPipe
    logging.info("Setting execute permissions for FragPipe...")

    execute_permissions_mask = stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH

    fragpipe_dir_glob = os.path.join(extraction_dir, "fragpipe-*")
    fragpipe_dirs = glob.glob(fragpipe_dir_glob)
    if len(fragpipe_dirs) == 0:
        logging.info(f"Error: There is no FragPipe installation directory fitting the glob pattern {fragpipe_dir_glob}.")
        logging.info(f"Stopping the installation of FragPipe and its tools...")
        sys.exit(1)

    fragpipe_dir = fragpipe_dirs[0]

    globs = ["bin/fragpipe", "tools/diann/*/linux/diann-*", "tools/diann/*/linux/*.so*", "tools/diann_so/*.so*", "tools/percolator_*/linux/percolator",
             "tools/Philosopher/philosopher-*", "tools/PTMProphet/PTMProphetParser-*", "tools/SAINTexpress/SAINTexpress-*"]
    for executable_or_shared_lib_glob in globs:
        files = glob.glob(os.path.join(fragpipe_dir, executable_or_shared_lib_glob))
        for file in files:
            try:
                file_permission_mask = os.stat(file).st_mode & 0o777
                os.chmod(file, file_permission_mask | execute_permissions_mask)
            except Exception as e:
                logging.warning(f"Warning: An exception has occurred when setting execute permissions for '{file}'.")
                logging.warning(e)
                logging.warning("Continuing...")
                continue

    logging.info("Done!")

    while True:
        try:
            resp = input("Would you like to install MSFragger, IonQuant and diaTracer (required additional tools to run the FragPipe analysis) ([Y]/n)? ").strip().lower()
            if resp.startswith('n'):
                logging.error("Stopping the installation of the FragPipe tools...")
                sys.exit(1)
            elif resp.startswith('y') or resp == '':
                logging.info("Beginning installation of the FragPipe tools...")
                break
            logging.error(f"Error: Invalid response '{resp}'. Please enter 'y' or 'n' (case-insensitive).")
            continue
        except KeyboardInterrupt:
            logging.error("User has quit the installation script. Stopping the installation of the FragPipe tools...")
            sys.exit(1)

    tools_extraction_dir_glob = os.path.join(extraction_dir, "fragpipe-*", "tools")
    tools_extraction_dirs = glob.glob(tools_extraction_dir_glob)
    if len(tools_extraction_dirs) == 0:
        logging.error(f"Error: There is no FragPipe tools directory fitting the glob pattern {tools_extraction_dir_glob}.")
        logging.error(f"Stopping the installation of the FragPipe tools...")
        sys.exit(1)

    tools_extraction_dir = tools_extraction_dirs[0]

    msfragger_latest_version_url = "https://msfragger-upgrader.nesvilab.org/upgrader/latest_version.php"
    ionquant_latest_version_url = "https://msfragger-upgrader.nesvilab.org/ionquant/latest_version.php"
    diatracer_latest_version_url = "https://msfragger-upgrader.nesvilab.org/diatracer/latest_version.php"

    logging.info("Fetching the latest versions of MSFragger, IonQuant and diaTracer...")

    while True:
        try:
            with urlopen(msfragger_latest_version_url) as f:
                msfragger_version = f.read().decode().strip()
            logging.info(f"MSFragger: {msfragger_version}")
            with urlopen(ionquant_latest_version_url) as f:
                ionquant_version = f.read().decode().strip()
            logging.info(f"IonQuant: {ionquant_version}")
            with urlopen(diatracer_latest_version_url) as f:
                diatracer_version = f.read().decode().strip()
            logging.info(f"diaTracer: {diatracer_version}")
            break
        except KeyboardInterrupt:
            logging.error("User has quit the installation script. Stopping the installation of the FragPipe tools...")
            sys.exit(1)
        except Exception as e:
            logging.error("Error: An exception has occurred when fetching the latest versions of the FragPipe tools.")
            logging.error(e)
            logging.error("Retrying...")
            continue

    logging.info("To install the FragPipe tools, first head to https://msfragger-upgrader.nesvilab.org/upgrader/.")
    logging.info("Then, enter your first name, last name, academic email address and academic institution, check all tickboxes for the academic license, license agreement and SDK library distribution conditions, and click on the 'Download' button.")
    logging.info("Next, wait for an email from no-reply@fragpipe.info at the email address you have specified. It should contain a download link that has the following format:")
    logging.info("https://msfragger-upgrader.nesvilab.org/upgrader/download.php?token=<6-digit token>&download=<version>%24zip")
    logging.info("The 6-digit token can be used to download MSFragger, IonQuant and diaTracer.")

    while True:
        try:
            token = input("Please enter your token (6 digits): ").strip()
            if len(token) != 6:
                logging.error("Error: Please enter exactly 6 characters for the token.")
                continue
            if not token.isdigit():
                logging.error("Error: Please ensure all 6 characters entered for the token are digits.")
                continue
            msfragger_download_url = f"https://msfragger-upgrader.nesvilab.org/upgrader/download.php?token={token}&download={quote_plus(msfragger_version)}$zip"
            ionquant_download_url = f"https://msfragger-upgrader.nesvilab.org/ionquant/download.php?token={token}&download={quote_plus(ionquant_version)}$zip"
            diatracer_download_url = f"https://msfragger-upgrader.nesvilab.org/diatracer/download.php?token={token}&download={quote_plus(diatracer_version)}$zip"
            logging.info(f"Using token '{token}'...")
            logging.info(f"Downloading the latest version of MSFragger from {msfragger_download_url}...")
            with urlopen(msfragger_download_url) as f:
                msfragger_zip = f.read()
            if not msfragger_zip.startswith(zip_magic_number):
                logging.error("Error: The token you have entered may be invalid or have expired. Please enter a valid token or request a new token by re-submitting your details at https://msfragger-upgrader.nesvilab.org/upgrader/ and waiting for an email containing the new token.")
                continue
            logging.info("MSFragger downloaded!")
            logging.info(f"Downloading the latest version of IonQuant from {ionquant_download_url}...")
            with urlopen(ionquant_download_url) as f:
                ionquant_zip = f.read()
            if not ionquant_zip.startswith(zip_magic_number):
                logging.error("Error: The token you have entered may be invalid or have expired. Please enter a valid token or request a new token by re-submitting your details at https://msfragger-upgrader.nesvilab.org/upgrader/ and waiting for an email containing the new token.")
                continue
            logging.info("IonQuant downloaded!")
            logging.info(f"Downloading the latest version of diaTracer from {diatracer_download_url}...")
            with urlopen(diatracer_download_url) as f:
                diatracer_zip = f.read()
            if not diatracer_zip.startswith(zip_magic_number):
                logging.error("Error: The token you have entered may be invalid or have expired. Please enter a valid token or request a new token by re-submitting your details at https://msfragger-upgrader.nesvilab.org/upgrader/ and waiting for an email containing the new token.")
                continue
            logging.info("diaTracer downloaded!")
            break
        except KeyboardInterrupt:
            logging.error("User has quit the installation script. Stopping the installation of the FragPipe tools...")
            sys.exit(1)
        except Exception as e:
            logging.error("Error: An exception has occurred when downloading the FragPipe tools.")
            logging.error(e)
            logging.error("Retrying...")
            continue

    # Extract MSFragger
    logging.info("Unzipping MSFragger...")

    with io.BytesIO(msfragger_zip) as msfragger_zip_stream:
        with zipfile.ZipFile(msfragger_zip_stream, 'r') as msfragger_zip_file:
            msfragger_zip_file.extractall(tools_extraction_dir)

    del msfragger_zip

    logging.info("Done!")

    # Extract IonQuant
    logging.info("Unzipping IonQuant...")

    with io.BytesIO(ionquant_zip) as ionquant_zip_stream:
        with zipfile.ZipFile(ionquant_zip_stream, 'r') as ionquant_zip_file:
            ionquant_zip_file.extractall(tools_extraction_dir)

    del ionquant_zip

    logging.info("Done!")

    # Extract diaTracer
    logging.info("Unzipping diaTracer...")

    with io.BytesIO(diatracer_zip) as diatracer_zip_stream:
        with zipfile.ZipFile(diatracer_zip_stream, 'r') as diatracer_zip_file:
            diatracer_zip_file.extractall(tools_extraction_dir)

    logging.info("Done!")

    # Set execute permissions for MSFragger
    logging.info("Setting execute permissions for MSFragger...")

    msfragger_shared_lib_glob = os.path.join(fragpipe_dir, "tools/MSFragger-*/ext/*/*.so")
    msfragger_shared_libs = glob.glob(msfragger_shared_lib_glob)
    for file in msfragger_shared_libs:
        try:
            file_permission_mask = os.stat(file).st_mode & 0o777
            os.chmod(file, file_permission_mask | execute_permissions_mask)
        except Exception as e:
            logging.warning(f"Warning: An exception has occurred when setting execute permissions for '{file}'.")
            logging.warning(e)
            logging.warning("Continuing...")
        continue

    logging.info("Done!")

    logging.info("FragPipe installation complete.")
    logging.info(f"FragPipe extracted into {extraction_dir}.")
    logging.info(f"FragPipe tools extracted into {tools_extraction_dir}.")

if __name__ == "__main__":
    main()
